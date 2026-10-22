#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileCopyrightText: Copyright © 2025-2026 University of Bern
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@unibe.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import numpy as np

def select_region_around_3UTR(row, start_region_around=300, end_region_around=400):
    """
    Selects a region around the 3' UTR of a genomic feature based on the strand orientation.
    Parameters:
        row (pandas.Series): A row containing at least 'start', 'end', and 'strand' fields.
        start_region_around (int, optional): The offset from the start position (for '+' strand) or end position (for '-' strand) to define the region's start. Default is 300.
        end_region_around (int, optional): The offset from the start position (for '+' strand) or end position (for '-' strand) to define the region's end. Default is 400.
    Returns:
        pandas.Series: The input row with updated 'start' and 'end' fields representing the selected region around the 3' UTR.
    """
    

    old_start = row['start']
    old_end = row['end']

    if row['strand'] == '+':
        row['start'] = old_start + start_region_around
        row['end'] = old_start + end_region_around
    else:
        row['start'] = old_end - end_region_around
        row['end'] = old_end - start_region_around
    return row

import re

def is_suspected_to_not_be_5utr(start, end, start_start_codon, end_start_codon, strand):
    """True if a 5' UTR interval is not entirely upstream of the start codon of its gene (False when the gene has no
    start codon).
    """
    if strand == "+":
        if np.isnan(end_start_codon):
            return False
        if start <= end_start_codon:
            if end <= end_start_codon:
                return False
            else:
                return True
        else:
            return True
    else:
        if np.isnan(start_start_codon):
            return False
        if end >= start_start_codon:
            if start >= start_start_codon:
                return False
            else:
                return True
        else:
            return True


def is_suspected_to_not_be_3utr(start, end, start_stop_codon, end_stop_codon, strand):
    """True if a 3' UTR interval is not entirely downstream of the stop codon of its gene (False when the gene has no
    stop codon).
    """
    if strand == "+":
        if np.isnan(start_stop_codon):
            return False
        if end >= start_stop_codon:
            if start >= start_stop_codon:
                return False
            else:
                return True
        else:
            return True
    else:
        if np.isnan(end_stop_codon):
            return False
        if start <= end_stop_codon:
            if end <= end_stop_codon:
                return False
            else:
                return True
        else:
            return True


def remove_suspicious_5_utrs(whole_annotation, utr_annotation):
    """Remove the 5' UTRs that are not upstream of the start codon of their gene.

    Args:
        whole_annotation: GTF table of the genes, with their start_codon features.
        utr_annotation: 5' UTR table whose `feature` is "<gene_id>_utr..."; a `gene_id` column is added to it.

    Returns:
        pd.DataFrame: the 5' UTRs kept, with the start codon coordinates of their gene.
    """
    utr_annotation["gene_id"] = utr_annotation["feature"].apply(
        lambda x: x.split("_utr")[0]
    )
    annot_advanced = (
        whole_annotation.groupby(["gene_id", "feature", "strand"])
        .agg({"start": min, "end": max})
        .reset_index()
    )
    utr_codon_info = annot_advanced[annot_advanced["feature"] == "start_codon"][
        ["strand", "gene_id", "start", "end"]
    ].merge(
        utr_annotation,
        left_on=["gene_id", "strand"],
        how="right",
        right_on=["gene_id", "strand"],
        suffixes=("_start_codon", ""),
    )
    utr_codon_info["is_suspicious"] = utr_codon_info.apply(
        lambda row: is_suspected_to_not_be_5utr(
            row["start"],
            row["end"],
            row["start_start_codon"],
            row["end_start_codon"],
            row["strand"],
        ),
        axis=1,
    )
    suspicious_utrs = utr_codon_info[utr_codon_info["is_suspicious"]].index
    return utr_codon_info[~(utr_codon_info.index.isin(suspicious_utrs))]


def remove_suspicious_3_utrs(whole_annotation, utr_annotation):
    """Remove the 3' UTRs that are not downstream of the stop codon of their gene.

    Args:
        whole_annotation: GTF table of the genes, with their stop_codon features.
        utr_annotation: 3' UTR table whose `feature` is "<gene_id>_utr..."; a `gene_id` column is added to it.

    Returns:
        pd.DataFrame: the 3' UTRs kept, with the stop codon coordinates of their gene.
    """
    utr_annotation["gene_id"] = utr_annotation["feature"].apply(
        lambda x: x.split("_utr")[0]
    )
    annot_advanced = (
        whole_annotation.groupby(["gene_id", "feature", "strand"])
        .agg({"start": min, "end": max})
        .reset_index()
    )
    utr_codon_info = annot_advanced[annot_advanced["feature"] == "stop_codon"][
        ["strand", "gene_id", "start", "end"]
    ].merge(
        utr_annotation,
        left_on=["gene_id", "strand"],
        how="right",
        right_on=["gene_id", "strand"],
        suffixes=("_stop_codon", ""),
    )
    utr_codon_info["is_suspicious"] = utr_codon_info.apply(
        lambda row: is_suspected_to_not_be_3utr(
            row["start"],
            row["end"],
            row["start_stop_codon"],
            row["end_stop_codon"],
            row["strand"],
        ),
        axis=1,
    )
    suspicious_utrs = utr_codon_info[utr_codon_info["is_suspicious"]].index

    return utr_codon_info[~(utr_codon_info.index.isin(suspicious_utrs))]
