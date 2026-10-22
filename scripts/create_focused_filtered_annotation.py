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
import os

from gtfparse import read_gtf

from genomemanager.annotation import focus_annotation_focus
from genomemanager.annotation import remove_overlaps_in_annotations


def main():
    """Create, from an isoform annotation (GTF), its 3' end focused version and the same without overlapping isoforms
    (`<name>_<nt>_focus_<suffix>.gtf` and `<name>_<nt>_focus_no_overlaps_<suffix>.gtf` in the output folder).
    """
    FLAGS = argparse.ArgumentParser(
        description="Create annotation focused without overlaps."
    )
    FLAGS.add_argument(
        "-a",
        "--annotation_gtf",
        help="<Required> Path to original GTF file",
        required=True,
    )
    FLAGS.add_argument("-f", "--focus_type", help="Type of focus", default="end")
    FLAGS.add_argument(
        "-n", "--nucleotides", help="Number of nucleotides for the focus", default=300
    )
    FLAGS.add_argument(
        "-o", "--output", help="Directory where to save the final GTFs", required=True
    )
    FLAGS.add_argument("-s", "--suffix", help="Suffix for the new GTFs", default="")
    args = FLAGS.parse_args()

    file_name = os.path.basename(args.annotation_gtf).split(".")[0]
    gtf_no_focus = read_gtf(args.annotation_gtf)
    # TODO: should be modified: --focus_type and --nucleotides are not passed to focus_annotation_focus (the defaults,
    # "end" and 300 nt, are always used, while --nucleotides is used in the output file names).
    gtf_focus = focus_annotation_focus(
        annotation=gtf_no_focus,
        save_path=os.path.join(
            args.output, f"{file_name}_{args.nucleotides}_focus_{args.suffix}.gtf"
        ),
    )
    print(
        f"Annotation with {args.nucleotides} focus created! Contains {len(gtf_focus)} isoforms."
    )

    gtf_focus_no_overlaps = remove_overlaps_in_annotations(
        annotation_path=os.path.join(
            args.output, f"{file_name}_{args.nucleotides}_focus_{args.suffix}.gtf"
        ),
        save_path=os.path.join(
            args.output,
            f"{file_name}_{args.nucleotides}_focus_no_overlaps_{args.suffix}.gtf",
        ),
    )

    print(
        f"Annotation with focus and no overlaps created! Contains {len(gtf_focus_no_overlaps)} isoforms."
    )


if __name__ == "__main__":
    main()
