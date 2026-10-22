#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileCopyrightText: Copyright © 2025-2026 University of Bern
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@unibe.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import csv
import itertools

import pandas as pd
from gtfparse import read_gtf
import re


def reorder_gtf(df_gtf):
    """Columns of a GTF table in GTF order (seqname, source, feature, start, end, score, strand, frame, attributes)."""
    return df_gtf[
        [
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attributes",
        ]
    ]


def focus_annotation_focus(
    annotation, focus_type: str = "end", nt: int = 300, save_path: str = None
):
    """Restrict every isoform of an annotation to its last `nt` nucleotides (3' end, along the strand), e.g. to count
    reads of 3' end sequencing.

    Args:
        annotation: GTF table (one row per isoform, with `start`, `end` and `strand`); the other non-GTF columns become
            attributes.
        focus_type: "end" focuses on the 3' end. Defaults to "end".
        nt: length of the focus. Defaults to 300. Isoforms shorter than `nt` are kept whole.
        save_path: if given, the focused annotation is written there (GTF).

    Returns:
        pd.DataFrame: the focused annotation, in GTF format.
    """
    def end_focus_coordinates(strand, start, end, nt):
        """(start, end) of the last `nt` nucleotides of an interval along the strand (the whole interval if it is
        shorter).
        """
        if (end - start) <= nt:
            start_focus = start
            end_focus = end
        else:
            if strand == "+":
                start_focus = end - nt
                end_focus = end
            else:
                start_focus = start
                end_focus = start + nt
        return start_focus, end_focus

    annotation_focus = annotation.copy()
    if focus_type == "end":
        annotation_focus["start_focus"] = annotation_focus.apply(
            lambda row: end_focus_coordinates(
                row["strand"], row["start"], row["end"], nt
            )[0],
            axis=1,
        )
        annotation_focus["end_focus"] = annotation_focus.apply(
            lambda row: end_focus_coordinates(
                row["strand"], row["start"], row["end"], nt
            )[1],
            axis=1,
        )
    # TODO: should be modified: only focus_type="end" works. This branch swaps the end-focus coordinates (start > end
    # in the output) instead of focusing on the 5' end; implement a start focus or reject other values.
    else:
        annotation_focus["end_focus"] = annotation_focus.apply(
            lambda row: end_focus_coordinates(
                row["strand"], row["start"], row["end"], nt
            )[0],
            axis=1,
        )
        annotation_focus["start_focus"] = annotation_focus.apply(
            lambda row: end_focus_coordinates(
                row["strand"], row["start"], row["end"], nt
            )[1],
            axis=1,
        )
    annotation_focus.drop(["start", "end"], axis=1, inplace=True)
    annotation_focus.rename(
        columns={"start_focus": "start", "end_focus": "end"}, inplace=True
    )
    annotation_focus = format_dataframe_as_gtf(annotation_focus)

    if save_path is not None:
        annotation_focus.to_csv(
            save_path,
            index=False,
            sep="\t",
            quoting=csv.QUOTE_NONE,
            escapechar="\\",
            header=False,
        )
    return annotation_focus


def adjust_positions_prox(strand, start_annot, end_annot, predicted_proximal_APA):
    """(start, end) of the proximal isoform predicted by DaPars: from the annotated start to the predicted proximal PAS
    on the + strand, from the predicted proximal PAS to the annotated end on the - strand.
    """
    if strand == "+":
        start = start_annot
        end = predicted_proximal_APA
    else:
        start = predicted_proximal_APA
        end = end_annot
    return start, end


def reannotate_transcripts(
    df,
    gene_id: str = "gene_id",
    strand_col: str = "strand",
    start_col: str = "start",
    end_col: str = "end",
):
    """Number the isoforms of every gene with several PAS from proximal to distal: transcript_id = "<gene_id>.<rank>"
    (ranked by end on the + strand, by decreasing start on the - strand).

    Args:
        df: annotation with `gene_id`, strand, start and end columns (index is reset).
        gene_id, strand_col, start_col, end_col: column names.

    Returns:
        pd.DataFrame: `df` with the new `transcript_id` values.
    """
    df.reset_index(drop=True, inplace=True)
    multiple_pas_genes = list(
        df.groupby("gene_id")
        .count()[df.groupby("gene_id").count()["seqname"] > 1]
        .index
    )
    for gene in multiple_pas_genes:
        strand = list(df[df[gene_id] == gene][strand_col].values)[0]
        if strand == "-":
            for j, (i, row) in enumerate(
                df[df[gene_id] == gene]
                .sort_values(start_col, ascending=False)
                .iterrows()
            ):
                df.loc[i, "transcript_id"] = f"{df.loc[i, 'gene_id']}.{str(j+1)}"
        else:
            for j, (i, row) in enumerate(
                df[df[gene_id] == gene].sort_values(end_col, ascending=True).iterrows()
            ):
                df.loc[i, "transcript_id"] = f"{df.loc[i, 'gene_id']}.{str(j+1)}"
    return df


def get_gtf_annotation_from_DaPars_output(df: pd.DataFrame, save_path=None):
    """GTF annotation of the proximal and distal 3' UTR isoforms found by DaPars.

    Each DaPars row (gene, locus, predicted proximal APA site) gives a distal isoform (the whole locus) and a proximal
    isoform (up to the predicted site); isoforms are numbered per gene from proximal to distal.

    Args:
        df: DaPars results (`Gene`, `Loci`, `Predicted_Proximal_APA` columns).
        save_path: if given, the annotation is written there (GTF).

    Returns:
        pd.DataFrame: the annotation, in GTF format (duplicate intervals removed).
    """
    df["gene_id"] = df["Gene"].apply(lambda x: x.split("|")[1])
    df["strand"] = df["Gene"].apply(lambda x: x.split("|")[-1])
    df["start"] = df["Loci"].apply(lambda x: int(x.split("-")[0].split(":")[1]))
    df["end"] = df["Loci"].apply(lambda x: int(x.split("-")[-1]))
    df_dapars_proxs = df.copy()
    df_dapars_proxs["start_"] = df_dapars_proxs.apply(
        lambda row: adjust_positions_prox(
            row["strand"], row["start"], row["end"], row["Predicted_Proximal_APA"]
        )[0],
        axis=1,
    )
    df_dapars_proxs["end_"] = df_dapars_proxs.apply(
        lambda row: adjust_positions_prox(
            row["strand"], row["start"], row["end"], row["Predicted_Proximal_APA"]
        )[1],
        axis=1,
    )
    df_dapars_proxs.drop(["start", "end"], axis=1, inplace=True)
    df_dapars_proxs.rename(columns={"start_": "start", "end_": "end"}, inplace=True)

    df_dapars_dist = df.copy()

    dapars_annotation = pd.concat(
        [df_dapars_proxs.reset_index(drop=True), df_dapars_dist.reset_index(drop=True)]
    )

    dapars_annotation["seqname"] = dapars_annotation["Loci"].apply(
        lambda x: x.split(":")[0]
    )
    dapars_annotation["feature"] = "exon"
    dapars_annotation["source"] = "rtracklayer"
    dapars_annotation["frame"] = "."
    dapars_annotation["score"] = "."

    dapars_annotation["transcript_id"] = dapars_annotation["gene_id"] + ".1"
    dapars_annotation = reannotate_transcripts(dapars_annotation)

    dapars_annotation["attributes"] = dapars_annotation.apply(
        lambda row: f"uniqueID \"{row['transcript_id']}\"; \
        iso {str(row['transcript_id'].split('.')[-1])}; \
        transcript_id \"{row['transcript_id']}\"; \
        gene_id \"{row['gene_id']}\"; \
        gene_name \"{row['gene_id']}\"; \
        gID \"{row['gene_id']}\"; \
        ord {str(row['transcript_id'].split('.')[-1])}",
        axis=1,
    )
    dapars_annotation = dapars_annotation[
        [
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attributes",
        ]
    ]
    dapars_annotation.drop_duplicates(["start", "end"], keep="first", inplace=True)

    if save_path is not None:
        dapars_annotation.to_csv(
            save_path,
            index=False,
            sep="\t",
            quoting=csv.QUOTE_NONE,
            escapechar="\\",
            header=False,
        )
    return dapars_annotation


def format_dataframe_as_gtf(df):
    """Format a table as a GTF: every column that is not a GTF field becomes a `key "value";` attribute, score and frame
    are set to ".", and the columns are put in GTF order.

    Args:
        df: table with the GTF fields and any attribute columns.

    Returns:
        pd.DataFrame: the GTF table.
    """
    def create_attribute_string(row, cols):
        """GTF attribute string (`key "value";` for each column of `cols`) of one row."""
        for i, col in enumerate(cols):
            if i == 0:
                attribute = f'{col} "{row[col]}";'
            else:
                attribute = f'{attribute} {col} "{row[col]}";'
        return attribute

    cols = df.columns
    attributes_cols = [
        col
        for col in cols
        if col
        not in [
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
        ]
    ]

    df_formatted = df.copy()
    df_formatted["attributes"] = df_formatted.apply(
        lambda row: create_attribute_string(row, attributes_cols), axis=1
    )
    df_formatted.drop(
        attributes_cols,
        axis=1,
        inplace=True,
    )

    df_formatted["score"] = "."
    df_formatted["frame"] = "."

    df_formatted = reorder_gtf(df_formatted)
    return df_formatted


def remove_overlaps_in_annotations(annotation_path, save_path: str = None):
    """Remove overlapping isoforms of the same gene: going along the strand, an isoform that starts before the end of
    the previous kept isoform is dropped.

    Args:
        annotation_path: GTF annotation.
        save_path: if given, the filtered annotation is written there (GTF).

    Returns:
        pd.DataFrame: the filtered annotation, in GTF format.
    """
    annotation = read_gtf(annotation_path)
    annotation_pos = annotation[annotation["strand"] == "+"]
    annotation_neg = annotation[annotation["strand"] == "-"]

    for gene in list(set(list(annotation_pos["gene_id"].values))):
        # positive strand
        df_gene = annotation_pos[annotation_pos["gene_id"] == gene].sort_values(
            "start", ascending=True
        )
        for j, (i, row) in enumerate(df_gene.iterrows()):
            if j == 0:
                first_end = row["end"]
            else:
                if row["start"] < first_end:
                    annotation_pos.drop(i, inplace=True)
                else:
                    first_end = row["end"]

    # negative strand
    for gene in list(set(list(annotation_neg["gene_id"].values))):
        df_gene = annotation_neg[annotation_neg["gene_id"] == gene].sort_values(
            "start", ascending=False
        )
        for j, (i, row) in enumerate(df_gene.iterrows()):
            if j == 0:
                first_end = row["start"]
            else:
                if row["end"] > first_end:
                    annotation_neg.drop(i, inplace=True)
                else:
                    first_end = row["start"]

    annotation_filtered = format_dataframe_as_gtf(
        pd.concat([annotation_pos, annotation_neg])
    )

    if save_path is not None:
        annotation_filtered.to_csv(
            save_path,
            index=False,
            sep="\t",
            quoting=csv.QUOTE_NONE,
            escapechar="\\",
            header=False,
        )

    return annotation_filtered


# Get the isoforms who overlap with ground truth
def check_overlap(df_, row_to_test):
    """IDs (`uniqueID`) of the isoforms of `df_` (other than `row_to_test`) with a start or an end inside `row_to_test`.

    Returns:
        tuple: (uniqueID of `row_to_test`, list of the overlapping uniqueIDs).
    """
    df = df_.drop(row_to_test.name)
    overlap = df.apply(
        lambda x: (x["start"] >= row_to_test["start"])
        & (x["start"] <= row_to_test["end"])
        | (x["end"] >= row_to_test["start"]) & (x["end"] <= row_to_test["end"]),
        axis=1,
    )
    if len(overlap) == 0:
        list_overlaps = []
    else:
        list_overlaps = list(df[overlap]["uniqueID"].values)
    return row_to_test["uniqueID"], list_overlaps


def keep_annotated_isoforms(myUTRs, ends_to_check, end_threshold=50, strand="+"):
    """Isoforms of `myUTRs` whose 3' end lies within `end_threshold` nt of one of the annotated ends `ends_to_check`.

    Args:
        myUTRs: isoforms of one gene (`start`, `end`, `uniqueID`).
        ends_to_check: annotated 3' end positions (RefSeq / Ensembl).
        end_threshold: tolerance in nt. Defaults to 50.
        strand: "+" (3' end = `end`) or "-" (3' end = `start`). Defaults to "+".

    Returns:
        pd.DataFrame: the matching isoforms, sorted from distal to proximal.
    """
    annotated_isoforms = []
    if strand == "+":
        for end in ends_to_check:
            is_close_to_end = (myUTRs["end"] >= end - end_threshold) & (
                myUTRs["end"] <= end + end_threshold
            )
            row_close_to_end = myUTRs[is_close_to_end]

            if not row_close_to_end.empty:
                annotated_isoforms.append(row_close_to_end["uniqueID"].values[0])

        return myUTRs[myUTRs["uniqueID"].isin(annotated_isoforms)].sort_values(
            "start", ascending=False
        )
    else:
        for end in ends_to_check:
            is_close_to_end = (myUTRs["start"] >= end - end_threshold) & (
                myUTRs["start"] <= end + end_threshold
            )
            row_close_to_end = myUTRs[is_close_to_end]

            if not row_close_to_end.empty:
                annotated_isoforms.append(row_close_to_end["uniqueID"].values[0])

        return myUTRs[myUTRs["uniqueID"].isin(annotated_isoforms)].sort_values(
            "start", ascending=True
        )


def find_3UTR_Ensembl(ensembl, gene_name, strand):
    """
    Find the start and end coordinates of the 3' UTR for a given gene in the Ensembl annotation.

    Parameters:
    ensembl (DataFrame): The Ensembl annotation containing gene information.
    gene_name (str): The name of the gene.
    strand (str): The strand of the gene ("+" for forward strand, "-" for reverse strand).

    Returns:
    tuple: A tuple containing the start and end coordinates of the 3' UTR.

    """
    if strand == "+":
        start_3UTR = ensembl[
            (ensembl["gene_name"] == gene_name) & (ensembl["feature"] == "UTR")
        ][
            "start"
        ].max()  # get the 3' UTR start from Ensembl
        end_3UTR = ensembl[
            (ensembl["gene_name"] == gene_name) & (ensembl["feature"] == "UTR")
        ][
            "end"
        ].max()  # get the 3' UTR end from Ensembl
    else:
        start_3UTR = ensembl[
            (ensembl["gene_name"] == gene_name) & (ensembl["feature"] == "UTR")
        ][
            "start"
        ].min()  # this is technically the longest end; but the start in terms of coordinates
        end_3UTR = ensembl[
            (ensembl["gene_name"] == gene_name) & (ensembl["feature"] == "UTR")
        ][
            "end"
        ].min()  # this is technically the first start; but the end in terms of coordinates
    return start_3UTR, end_3UTR


def find_3UTR_RefSeq(refSeq, gene_name):
    """
    Find the start and end coordinates of the 3' UTR for a given gene in the RefSeq annotation.

    Parameters:
    refSeq (DataFrame): The RefSeq annotation containing gene information.
    gene_name (str): The name of the gene.
    strand (str): The strand of the gene ("+" for forward strand, "-" for reverse strand).

    Returns:
    tuple: A tuple containing the start and end coordinates of the 3' UTR.

    """

    # For RefSeq, no need to distinguish the strands because the 3UTR is already annotated as such
    start_3UTR = refSeq[
        (refSeq["gene_name"] == gene_name) & refSeq["feature"].str.contains("3UTR")
    ][
        "start"
    ].min()  # get the 3' UTR end from refSeq
    end_3UTR = refSeq[
        (refSeq["gene_name"] == gene_name) & refSeq["feature"].str.contains("3UTR")
    ][
        "end"
    ].max()  # get the 3' UTR start from refSeq
    return start_3UTR, end_3UTR


def find_largest_3UTR(ensembl, refSeq, gene_name, strand):
    """
    Finds the largest 3' UTR (untranslated region) for a given gene.

    Args:
        ensembl (DataFrame): The Ensembl annotation containing gene information.
        refSeq (DataFrame): The RefSeq annotation containing gene information.
        gene_name (str): The name of the gene.
        strand (str): The strand of the gene.

    Returns:
        tuple: A tuple containing the start and end positions of the largest 3' UTR.

    """
    start_3UTR_ensembl, end_3UTR_ensembl = find_3UTR_Ensembl(ensembl, gene_name, strand)
    start_3UTR_refseq, end_3UTR_refseq = find_3UTR_RefSeq(refSeq, gene_name)

    return min(start_3UTR_ensembl, start_3UTR_refseq), max(
        end_3UTR_ensembl, end_3UTR_refseq
    )


def get_ground_truth_one_gene(
    gene_name, refSeq, ensembl, myUTRs_focus, end_threshold=50, strand="+", margin=300
):
    """Ground-truth 3' UTR isoforms of one gene: isoforms whose 3' end matches a RefSeq or Ensembl 3' UTR end, completed
    by the non-overlapping isoforms, overlaps being removed going from the distal end; isoforms more than `margin` nt
    outside the largest annotated 3' UTR are dropped.

    Args:
        gene_name: gene.
        refSeq: RefSeq GTF table (`gene_id`, feature "3UTR").
        ensembl: Ensembl GTF table (`gene_name`, feature "UTR").
        myUTRs_focus: isoform annotation to filter (`gene_id`, `uniqueID`, `start`, `end`).
        end_threshold: tolerance in nt to match annotated ends. Defaults to 50.
        strand: strand of the gene. Defaults to "+".
        margin: tolerance in nt around the largest annotated 3' UTR. Defaults to 300.

    Returns:
        tuple: (ground-truth isoforms, the other isoforms of the gene).
    """
    if strand == "+":
        end_name = "end"
        ascending = False
    else:
        end_name = "start"
        ascending = True
    ref_seq_end = refSeq[
        (refSeq["gene_id"] == gene_name) & (refSeq["feature"] == "3UTR")
    ][end_name]
    ensembl_end = ensembl[
        (ensembl["gene_name"] == gene_name) & (ensembl["feature"] == "UTR")
    ][end_name]
    myUTRS = myUTRs_focus[myUTRs_focus["gene_id"] == gene_name]

    # Concatenate the two series
    concatenated_series = pd.concat([ref_seq_end, ensembl_end])

    # sort the values in descending order and get the unique values
    ends_to_check = concatenated_series.sort_values(ascending=ascending).unique()

    # Create ground truth -- keep annotated isoforms
    annots_all = keep_annotated_isoforms(
        myUTRS, ends_to_check, end_threshold, strand=strand
    )

    annots_ensembl = keep_annotated_isoforms(
        myUTRS,
        ensembl_end.sort_values(ascending=ascending).unique(),
        end_threshold,
        strand=strand,
    )

    # Get the isoforms who overlap with ground truth
    overlapping_isoforms = {}
    for _, row in annots_all.iterrows():
        key, vals = check_overlap(myUTRS, row)
        overlapping_isoforms[key] = vals
    # to_remove = list(itertools.chain.from_iterable(to_remove))
    # Do not throw away the annotated isoforms even if they overlap
    overlapping_not_annots = [
        to_rem
        for to_rem in list(
            itertools.chain.from_iterable(list(overlapping_isoforms.values()))
        )
        if to_rem not in annots_all["uniqueID"].values
    ]
    overlapping_annots = [
        to_rem
        for to_rem in list(
            itertools.chain.from_iterable(list(overlapping_isoforms.values()))
        )
        if to_rem in annots_all["uniqueID"].values
    ]

    # If there are overlapping annotations; keep the Ensembl one;
    # or the longest one if they belong to the same annotation
    isoforms_to_keep = []
    isoforms_to_remove = []
    for overlap_iso in overlapping_annots:
        list_to_test = [overlap_iso] + overlapping_isoforms[overlap_iso]
        for i, iso in enumerate(list_to_test):
            if i == 0:
                iso_to_keep = iso
            else:
                is_ensembl_keep = iso_to_keep in annots_ensembl["uniqueID"].values
                is_ensembl_to_test = iso in annots_ensembl["uniqueID"].values
                if strand == "+":
                    is_longer_to_test = (
                        myUTRS[myUTRS["uniqueID"] == iso]["end"].values[0]
                        > myUTRS[myUTRS["uniqueID"] == iso_to_keep]["end"].values[0]
                    )
                else:
                    is_longer_to_test = (
                        myUTRS[myUTRS["uniqueID"] == iso]["start"].values[0]
                        < myUTRS[myUTRS["uniqueID"] == iso_to_keep]["start"].values[0]
                    )
                if is_ensembl_keep:
                    if is_ensembl_to_test:
                        if is_longer_to_test:
                            iso_to_keep = iso
                else:
                    if is_ensembl_to_test:
                        iso_to_keep = iso
        isoforms_to_keep.append(iso_to_keep)
        isoforms_to_remove += [iso for iso in list_to_test if iso != iso_to_keep]
    isoforms_to_keep = list(set(isoforms_to_keep))
    isoforms_to_remove = list(set(isoforms_to_remove))

    # Remove overlapping isoforms not annotated
    ground_truth = myUTRS[~myUTRS["uniqueID"].isin(overlapping_not_annots)].sort_values(
        "start", ascending=ascending
    )
    # Remove overlapping isoforms annotated
    ground_truth = ground_truth[
        ~ground_truth["uniqueID"].isin(isoforms_to_remove)
    ].sort_values("start", ascending=ascending)

    # Now remove the remaining overlapping isoforms

    if strand == "+":
        for j, (i, row) in enumerate(ground_truth.iterrows()):
            if j == 0:
                last_start = row["start"]
            else:
                if row["end"] > last_start:
                    ground_truth.drop(i, inplace=True)
                else:
                    last_start = row["start"]
    else:
        for j, (i, row) in enumerate(ground_truth.iterrows()):
            if j == 0:
                last_start = row["end"]
            else:
                if row["start"] < last_start:
                    ground_truth.drop(i, inplace=True)
                else:
                    last_start = row["end"]

    # Finally, remove the isoforms that fall outside of the annotated 3' UTRs,
    # with a margin of 300 nt
    start_3UTR, end_3UTR = find_largest_3UTR(
        ensembl, refSeq, gene_name=gene_name, strand=strand
    )
    ground_truth = ground_truth[
        (ground_truth["start"] <= end_3UTR + margin)
        & (ground_truth["end"] >= start_3UTR - margin)
    ]

    return (
        ground_truth,
        myUTRS[~myUTRS["uniqueID"].isin(ground_truth["uniqueID"].values)],
    )


def get_ground_truth_one_gene_side_prox(
    gene_name,
    refSeq,
    ensembl,
    myUTRs_focus,
    end_threshold=50,
    strand="+",
    margin=300,
):
    """Same as `get_ground_truth_one_gene`, but overlaps are removed going from the proximal end."""
    if strand == "+":
        end_name = "end"
        ascending = True
    else:
        end_name = "start"
        ascending = False
    ref_seq_end = refSeq[
        (refSeq["gene_id"] == gene_name) & (refSeq["feature"] == "3UTR")
    ][end_name]
    ensembl_end = ensembl[
        (ensembl["gene_name"] == gene_name) & (ensembl["feature"] == "UTR")
    ][end_name]
    myUTRS = myUTRs_focus[myUTRs_focus["gene_id"] == gene_name]

    # Concatenate the two series
    concatenated_series = pd.concat([ref_seq_end, ensembl_end])

    # sort the values in descending order and get the unique values
    ends_to_check = concatenated_series.sort_values(ascending=ascending).unique()

    # Create ground truth -- keep annotated isoforms
    annots_all = keep_annotated_isoforms(
        myUTRS, ends_to_check, end_threshold, strand=strand
    )

    annots_ensembl = keep_annotated_isoforms(
        myUTRS,
        ensembl_end.sort_values(ascending=ascending).unique(),
        end_threshold,
        strand=strand,
    )

    # Get the isoforms who overlap with ground truth
    overlapping_isoforms = {}
    for _, row in annots_all.iterrows():
        key, vals = check_overlap(myUTRS, row)
        overlapping_isoforms[key] = vals
    # to_remove = list(itertools.chain.from_iterable(to_remove))
    # Do not throw away the annotated isoforms even if they overlap
    overlapping_not_annots = [
        to_rem
        for to_rem in list(
            itertools.chain.from_iterable(list(overlapping_isoforms.values()))
        )
        if to_rem not in annots_all["uniqueID"].values
    ]
    overlapping_annots = [
        to_rem
        for to_rem in list(
            itertools.chain.from_iterable(list(overlapping_isoforms.values()))
        )
        if to_rem in annots_all["uniqueID"].values
    ]

    # If there are overlapping annotations; keep the Ensembl one;
    # or the longest one if they belong to the same annotation
    isoforms_to_keep = []
    isoforms_to_remove = []
    for overlap_iso in overlapping_annots:
        list_to_test = [overlap_iso] + overlapping_isoforms[overlap_iso]
        for i, iso in enumerate(list_to_test):
            if i == 0:
                iso_to_keep = iso
            else:
                is_ensembl_keep = iso_to_keep in annots_ensembl["uniqueID"].values
                is_ensembl_to_test = iso in annots_ensembl["uniqueID"].values
                if strand == "+":
                    is_longer_to_test = (
                        myUTRS[myUTRS["uniqueID"] == iso]["end"].values[0]
                        > myUTRS[myUTRS["uniqueID"] == iso_to_keep]["end"].values[0]
                    )
                else:
                    is_longer_to_test = (
                        myUTRS[myUTRS["uniqueID"] == iso]["start"].values[0]
                        < myUTRS[myUTRS["uniqueID"] == iso_to_keep]["start"].values[0]
                    )
                if is_ensembl_keep:
                    if is_ensembl_to_test:
                        if is_longer_to_test:
                            iso_to_keep = iso
                else:
                    if is_ensembl_to_test:
                        iso_to_keep = iso
        isoforms_to_keep.append(iso_to_keep)
        isoforms_to_remove += [iso for iso in list_to_test if iso != iso_to_keep]
    isoforms_to_keep = list(set(isoforms_to_keep))
    isoforms_to_remove = list(set(isoforms_to_remove))

    # Remove overlapping isoforms not annotated
    ground_truth = myUTRS[~myUTRS["uniqueID"].isin(overlapping_not_annots)].sort_values(
        "start", ascending=ascending
    )
    # Remove overlapping isoforms annotated
    ground_truth = ground_truth[
        ~ground_truth["uniqueID"].isin(isoforms_to_remove)
    ].sort_values("start", ascending=ascending)

    # Now remove the remaining overlapping isoforms

    if strand == "+":
        for j, (i, row) in enumerate(ground_truth.iterrows()):
            if j == 0:
                first_end = row["end"]
            else:
                if row["start"] < first_end:
                    ground_truth.drop(i, inplace=True)
                else:
                    first_end = row["end"]
    else:
        for j, (i, row) in enumerate(ground_truth.iterrows()):
            if j == 0:
                first_end = row["start"]
            else:
                if row["end"] > first_end:
                    ground_truth.drop(i, inplace=True)
                else:
                    first_end = row["start"]

    # Finally, remove the isoforms that fall outside of the annotated 3' UTRs,
    # with a margin of 300 nt
    start_3UTR, end_3UTR = find_largest_3UTR(
        ensembl, refSeq, gene_name=gene_name, strand=strand
    )
    ground_truth = ground_truth[
        (ground_truth["start"] <= end_3UTR + margin)
        & (ground_truth["end"] >= start_3UTR - margin)
    ]

    return (
        ground_truth,
        myUTRS[~myUTRS["uniqueID"].isin(ground_truth["uniqueID"].values)],
    )


def keep_annotated_and_remove_overlaps(
    annotation_path,
    refSeq_path,
    ensembl_path,
    remove_overlaps_starting_from_dist=True,
    save_path: str = None,
):
    """Filter an isoform annotation to its ground-truth isoforms (`get_ground_truth_one_gene`, or
    `get_ground_truth_one_gene_side_prox`), gene by gene, with an end tolerance of 30 nt.

    Args:
        annotation_path: isoform annotation (GTF) to filter.
        refSeq_path: RefSeq GTF.
        ensembl_path: Ensembl GTF.
        remove_overlaps_starting_from_dist: remove overlaps going from the distal (True) or the proximal (False) end.
            Defaults to True.
        save_path: if given, the filtered annotation is written there (GTF).

    Returns:
        pd.DataFrame: the filtered annotation, in GTF format.
    """
    # Read the annotations
    annotation = read_gtf(annotation_path)
    refSeq = read_gtf(refSeq_path)
    ensembl = read_gtf(ensembl_path)

    # Separate into positive and negative strand
    annotation_pos = annotation[annotation["strand"] == "+"]
    annotation_neg = annotation[annotation["strand"] == "-"]

    if remove_overlaps_starting_from_dist:
        # positive strand
        for gene in list(set(list(annotation_pos["gene_id"].values))):
            _, inverse_ground_truth = get_ground_truth_one_gene(
                gene, refSeq, ensembl, annotation_pos, end_threshold=30, strand="+"
            )
            annotation_pos.drop(inverse_ground_truth.index, inplace=True)

        # negative strand
        for gene in list(set(list(annotation_neg["gene_id"].values))):
            _, inverse_ground_truth = get_ground_truth_one_gene(
                gene, refSeq, ensembl, annotation_neg, end_threshold=30, strand="-"
            )
            annotation_neg.drop(inverse_ground_truth.index, inplace=True)

    else:
        # positive strand
        for gene in list(set(list(annotation_pos["gene_id"].values))):
            _, inverse_ground_truth = get_ground_truth_one_gene_side_prox(
                gene, refSeq, ensembl, annotation_pos, end_threshold=30, strand="+"
            )
            annotation_pos.drop(inverse_ground_truth.index, inplace=True)

        # negative strand
        for gene in list(set(list(annotation_neg["gene_id"].values))):
            _, inverse_ground_truth = get_ground_truth_one_gene_side_prox(
                gene, refSeq, ensembl, annotation_neg, end_threshold=30, strand="-"
            )
            annotation_neg.drop(inverse_ground_truth.index, inplace=True)

    annotation_filtered = format_dataframe_as_gtf(
        pd.concat([annotation_pos, annotation_neg])
    )

    if save_path is not None:
        annotation_filtered.to_csv(
            save_path,
            index=False,
            sep="\t",
            quoting=csv.QUOTE_NONE,
            escapechar="\\",
            header=False,
        )

    return annotation_filtered

