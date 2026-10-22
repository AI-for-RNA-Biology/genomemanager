#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileCopyrightText: Copyright © 2026 University of Bern
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@unibe.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import pickle  # nosec

import numpy as np
import pandas as pd
from gtfparse import read_gtf


def get_overlapping_genes(interval, others_intervals, other_gene_list):
    """Genes (from `other_gene_list`) whose interval in `others_intervals` overlaps `interval`, or NaN if none."""
    if len(other_gene_list[others_intervals.overlaps(interval)]) > 0:
        return list(other_gene_list[others_intervals.overlaps(interval)].values)
    else:
        return np.nan


def get_overlap_file(annotation_path, save_path=None):
    """Genes that overlap a gene on the opposite strand of the same chromosome (gene extent = min start to max end of
    its annotation lines).

    Args:
        annotation_path: GTF annotation.
        save_path: if given, the result is pickled there and None is returned.

    Returns:
        pd.DataFrame: one row per overlapping gene (index gene_id) with its `start`, `end`, `interval` and the list of
        genes it overlaps (`overlaps_with`); None if `save_path` is given.
    """
    annotations = read_gtf(annotation_path)
    possible_chr = list(set(list(annotations["seqname"])))
    overlaps = []
    for chromosome in possible_chr:
        annotations_chr = annotations[annotations["seqname"] == chromosome]
        positive_strand_genes = (
            annotations_chr[annotations_chr["strand"] == "+"]
            .groupby("gene_id")
            .agg({"start": min, "end": max})
        )
        negative_strand_genes = (
            annotations_chr[annotations_chr["strand"] == "-"]
            .groupby("gene_id")
            .agg({"start": min, "end": max})
        )
        positive_strand_genes["interval"] = positive_strand_genes.apply(
            lambda gene: pd.Interval(
                left=gene["start"], right=gene["end"], closed="both"
            ),
            axis=1,
        )
        negative_strand_genes["interval"] = negative_strand_genes.apply(
            lambda gene: pd.Interval(
                left=gene["start"], right=gene["end"], closed="both"
            ),
            axis=1,
        )
        positive_intervals = pd.arrays.IntervalArray(positive_strand_genes["interval"])
        negative_intervals = pd.arrays.IntervalArray(negative_strand_genes["interval"])
        positive_strand_genes["overlaps_with"] = positive_strand_genes[
            "interval"
        ].apply(
            lambda pos_interval: get_overlapping_genes(
                pos_interval, negative_intervals, negative_strand_genes.index
            )
        )
        negative_strand_genes["overlaps_with"] = negative_strand_genes[
            "interval"
        ].apply(
            lambda neg_interval: get_overlapping_genes(
                neg_interval, positive_intervals, positive_strand_genes.index
            )
        )
        overlaps.append(
            pd.concat([positive_strand_genes.dropna(), negative_strand_genes.dropna()])
        )
    overlaps = pd.concat(overlaps)

    if save_path is not None:
        with open(save_path, "wb") as f:
            pickle.dump(overlaps, f)
    else:
        return overlaps
