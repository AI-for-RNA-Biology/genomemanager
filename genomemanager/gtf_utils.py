#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileCopyrightText: Copyright © 2025-2026 University of Bern
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@unibe.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
"""
GTF file utilities for pandas-based GTF parsing.

This module provides functions to read GTF files using pandas instead of gtfparse,
avoiding polars compatibility issues.
"""

import pandas as pd
from typing import Optional
import re


def load_gtf_file(gtf_path: str, expand_attributes: bool = True) -> pd.DataFrame:
    """
    Load a GTF file using pandas and optionally expand attribute column.
    
    This function provides an alternative to gtfparse that avoids polars compatibility issues.
    
    Parameters:
    -----------
    gtf_path : str
        Path to the GTF file
    expand_attributes : bool, default True
        Whether to parse the attribute column and create separate columns for each attribute
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing the GTF data with standard GTF columns and optionally
        expanded attribute columns
        
    Examples:
    ---------
    >>> # Basic usage (expands attributes by default)
    >>> annotation = load_gtf_file("path/to/file.gtf")
    
    >>> # Load without expanding attributes
    >>> annotation = load_gtf_file("path/to/file.gtf", expand_attributes=False)
    """
    # Define standard GTF column names
    gtf_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    
    # Read GTF file with pandas
    df = pd.read_csv(
        gtf_path,
        sep='\t',
        header=None,
        names=gtf_columns,
        comment='#',
        dtype={
            'seqname': 'str',
            'source': 'str', 
            'feature': 'str',
            'start': 'int64',
            'end': 'int64',
            'score': 'str',  # Keep as string to handle '.' values
            'strand': 'str',
            'frame': 'str',  # Keep as string to handle '.' values
            'attribute': 'str'
        }
    )
    
    if expand_attributes:
        def parse_attributes(attr_string):
            """Parse GTF attribute string into a dictionary"""
            if pd.isna(attr_string):
                return {}
            
            attributes = {}
            # Split by semicolon and process each key-value pair
            for item in attr_string.split(';'):
                item = item.strip()
                if item and ' ' in item:
                    # Split on first space to separate key and value
                    key, value = item.split(' ', 1)
                    # Remove quotes if present
                    value = value.strip('"')
                    attributes[key] = value
            return attributes
        
        # Extract attributes and create new columns
        attr_dicts = df['attribute'].apply(parse_attributes)
        
        # Get all unique attribute keys
        all_keys = set()
        for attr_dict in attr_dicts:
            all_keys.update(attr_dict.keys())
        
        # Create columns for each attribute
        for key in sorted(all_keys):  # Sort for consistent column order
            df[key] = attr_dicts.apply(lambda x: x.get(key, None))
    
    return df


def read_gtf(filepath_or_buffer: str, expand_attribute_column: bool = True, **kwargs) -> pd.DataFrame:
    """
    Drop-in replacement for gtfparse.read_gtf() function.
    
    This function mimics the gtfparse.read_gtf() API for easier migration.
    
    Parameters:
    -----------
    filepath_or_buffer : str
        Path to the GTF file
    expand_attribute_column : bool, default True
        Whether to parse the attribute column and create separate columns for each attribute
    **kwargs
        Additional keyword arguments (ignored for compatibility)
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing the GTF data
    """
    return load_gtf_file(filepath_or_buffer, expand_attributes=expand_attribute_column)

def from_gtf_to_bed(gtf, bed):
    """Write the exons of a GTF file as a BED6 file (0-based start, name = transcript_id or ".", score ".").

    Args:
        gtf: input GTF file.
        bed: output BED file.
    """
    with open(gtf) as fin, open(bed, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            chrom = fields[0]
            start = int(fields[3]) - 1  # 0-based for BED
            end = fields[4]
            strand = fields[6]
            attributes = fields[8]

            # Try to extract transcript_id
            match = re.search(r'transcript_id "([^"]+)"', attributes)
            if match:
                name = match.group(1)
            else:
                name = "."

            fout.write(f"{chrom}\t{start}\t{end}\t{name}\t.\t{strand}\n")
