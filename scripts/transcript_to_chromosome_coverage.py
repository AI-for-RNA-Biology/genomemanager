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

import yaml

from genomemanager.coverage import from_transcript_to_chromosome_coverage


def main():
    """Convert a coverage BED file in transcript coordinates to chromosome coordinates, with the annotation of the
    transcripts; arguments from a YAML config (`--config_file`) or from the command line.
    """
    FLAGS = argparse.ArgumentParser(
        description="Convert transcript coverage to chromosome coverage."
    )
    FLAGS.add_argument("--config_file", help="Location of config file")
    FLAGS.add_argument("--annotation_path", help="Annotation path")
    FLAGS.add_argument("--bed_transcripts_path", help="Bedfile path")
    FLAGS.add_argument("--output", help="output file path")

    args = FLAGS.parse_args()

    if args.config_file is not None:
        with open(args.config_file, "r") as stream:
            config = yaml.safe_load(stream)

        bedfile_chr = from_transcript_to_chromosome_coverage(
            config["annotation_path"],
            config["bed_transcripts_path"],
            config["bed_transcript_col"],
            config["annotation_chromosome_col"],
            config["annotation_start_col"],
            config["bed_start_col"],
            config["bed_end_col"],
            config["annotation_transcript_col"],
            config["annotation_strand"],
        )
        bedfile_chr.to_csv(
            config["path_to_processed_bedfile"], sep="\t", header=False, index=False
        )
    elif (
        args.annotation_path is not None
        and args.bed_transcripts_path is not None
        and args.output
    ):
        bedfile_chr = from_transcript_to_chromosome_coverage(
            annotation_path=args.annotation_path,
            bed_transcripts_path=args.bed_transcripts_path,
        )

        bedfile_chr.to_csv(args.output, sep="\t", header=False, index=False)
    else:
        raise ValueError(
            "You need to enter a config file or the annotation path, the bed path and the output path."
        )


if __name__ == "__main__":
    main()
