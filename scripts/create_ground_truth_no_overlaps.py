#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileCopyrightText: Copyright © 2026 University of Bern
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@unibe.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import argparse

from genomemanager.annotation import keep_annotated_and_remove_overlaps


def main():
    """Filter an isoform annotation (GTF) to its ground-truth isoforms (matching RefSeq / Ensembl 3' UTR ends), removing
    overlaps from the distal end.
    """
    FLAGS = argparse.ArgumentParser(
        description="Create ground truth annotation focused without overlaps"
    )
    FLAGS.add_argument(
        "-a",
        "--annotation_gtf",
        help="<Required> Path to initial annotation GTF file",
        required=True,
    )

    FLAGS.add_argument(
        "-r",
        "--refseq_gtf",
        help="<Required> Path to RefSeq GTF file",
        required=True,
    )
    FLAGS.add_argument(
        "-e",
        "--ensembl_gtf",
        help="<Required> Path to Ensembl GTF file",
        required=True,
    )
    FLAGS.add_argument(
        "-o", "--output", help="Path where to save the final GTF", required=True
    )
    args = FLAGS.parse_args()

    gtf_corrected = keep_annotated_and_remove_overlaps(
        annotation_path=args.annotation_gtf,
        refSeq_path=args.refseq_gtf,
        ensembl_path=args.ensembl_gtf,
        save_path=args.output,
    )

    print(
        f"Annotation with with annotated ground truth and without overlaps \
            created! Contains {len(gtf_corrected)} isoforms."
    )


if __name__ == "__main__":
    main()
