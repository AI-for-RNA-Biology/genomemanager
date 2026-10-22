#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileCopyrightText: Copyright © 2026 University of Bern
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@unibe.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import pyranges as pr
from gtfparse import read_gtf


def adjust_start_position(
    strand, bed_start_col, bed_end_col, annotation_start_col, annotation_end_col
):
    """Genomic start of an interval given in transcript coordinates (BED start / end relative to the transcript), for a
    transcript annotated at [annotation_start_col, annotation_end_col] on `strand`.
    """
    if strand == "+":
        return annotation_start_col + bed_start_col
    else:
        return annotation_end_col - bed_end_col


def adjust_end_position(
    strand, bed_start_col, bed_end_col, annotation_start_col, annotation_end_col
):
    """Genomic end of an interval given in transcript coordinates (BED start / end relative to the transcript), for a
    transcript annotated at [annotation_start_col, annotation_end_col] on `strand`.
    """
    if strand == "+":
        return annotation_start_col + bed_end_col
    else:
        return annotation_end_col - bed_start_col


def from_transcript_to_chromosome_coverage(
    annotation_path: str,
    bed_transcripts_path: str,
    bed_transcript_col: str = "Chromosome",
    annotation_chromosome_col: str = "seqname",
    annotation_start_col: str = "start",
    annotation_end_col: str = "end",
    bed_start_col: str = "Start",
    bed_end_col: str = "End",
    annotation_transcript_col: str = "transcript_id",
    annotation_strand_col: str = "strand",
    split: bool = True,
):
    """Convert a .bed coverage file, containing the coverage along the transcripts, to a .bed coverage containing
    the coverage along the chromosomes.

    Args:
        annotation_path: path to the .gtf file containing the transcripts coordinates along the chromosomes.
        bed_transcripts_path: path to the .bed transcripts coverage file to be converted.
        bed_transcript_col: Name of the column containing the transcripts name in the original bedfile.
        annotation_chromosome_col: Name of the column in the annotation file containing the chromosome names.
        annotation_start_col: Name of the column in the annotation file containing the start position of the
        transcripts along the chromosomes.
        bed_start_col: Name of the column in the original bedfile containing the start positions along the transcripts.
        bed_end_col: Name of the column in the original bedfile containing the end positions along the transcripts.
        annotation_transcript_col: Name of the column in the annotation file containing the name of the transcripts.
        These should match with the names in *bed_transcript_col*.
        split: whether to split or not the name of the transcripts in the bedfile,
        by taking the first part before the first '|'

    Returns:
        pd.DataFrame: pandas dataframe representing the converted bedfile. Can be saved into a .bed.
    """

    # Read
    annotation = read_gtf(annotation_path)
    bedfile = pr.read_bed(bed_transcripts_path)
    bedfile_chr = bedfile.df.copy()

    if split:
        bedfile_chr[bed_transcript_col] = bedfile_chr[bed_transcript_col].apply(
            lambda x: x.split("|")[0]
        )
        print("After split:")
        print(bedfile_chr.head())

    # Adjust the positions
    chr_df = bedfile_chr.merge(
        annotation[
            [
                annotation_chromosome_col,
                annotation_end_col,
                annotation_start_col,
                annotation_strand_col,
                annotation_transcript_col,
            ]
        ],
        left_on=bed_transcript_col,
        right_on=annotation_transcript_col,
    )
    # chr_df[bed_start_col] = chr_df.apply(
    #     lambda row: row[bed_start_col] + row[annotation_start_col], axis=1
    # )
    # chr_df[bed_end_col] = chr_df.apply(
    #     lambda row: row[bed_end_col] + row[annotation_start_col], axis=1
    # )
    print("After merge:")
    print(chr_df.head())
    chr_df[f"{bed_start_col}_adj"] = chr_df.apply(
        lambda row: adjust_start_position(
            strand=row[annotation_strand_col],
            bed_start_col=row[bed_start_col],
            bed_end_col=row[bed_end_col],
            annotation_start_col=row[annotation_start_col],
            annotation_end_col=row[annotation_end_col],
        ),
        axis=1,
    )
    chr_df[f"{bed_end_col}_adj"] = chr_df.apply(
        lambda row: adjust_end_position(
            strand=row[annotation_strand_col],
            bed_start_col=row[bed_start_col],
            bed_end_col=row[bed_end_col],
            annotation_start_col=row[annotation_start_col],
            annotation_end_col=row[annotation_end_col],
        ),
        axis=1,
    )
    print("After adjustement:")
    print(chr_df.head())
    chr_df["Chromosome"] = chr_df[annotation_chromosome_col].values
    chr_df.drop([bed_start_col, bed_end_col], axis=1, inplace=True)
    chr_df.rename(
        columns={
            f"{bed_start_col}_adj": bed_start_col,
            f"{bed_end_col}_adj": bed_end_col,
        },
        inplace=True,
    )

    if len(chr_df[chr_df[bed_start_col] > chr_df[bed_end_col]]) > 0:
        raise Exception("Start position is bigger than end position")

    # chr_df = chr_df[list(bedfile.df.columns) + [annotation_strand]]
    chr_df = chr_df[list(bedfile.df.columns)]

    print("Before return:")
    print(chr_df.head())

    return chr_df
