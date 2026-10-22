#!/bin/bash
#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
PATH=$PATH:/path/to/samtools
# Get the bam file from the command line
DATA=$1
PREFIX=$2

# Forward strand.
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse  strand
#
samtools view -b -f 128 -F 16 $DATA > ${PREFIX}_fwd1.bam
samtools index ${PREFIX}_fwd1.bam

samtools view -b -f 80 $DATA > ${PREFIX}_fwd2.bam
samtools index ${PREFIX}_fwd2.bam

#
# Combine alignments that originate on the forward strand.
#
samtools merge -f ${PREFIX}_fwd.bam ${PREFIX}_fwd1.bam ${PREFIX}_fwd2.bam
rm ${PREFIX}_fwd1.bam ${PREFIX}_fwd2.bam
samtools index ${PREFIX}_fwd.bam
samtools sort ${PREFIX}_fwd.bam -o ${PREFIX}_fwd_sorted.bam
rm ${PREFIX}_fwd.bam

# Reverse strand
#
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
#
samtools view -b -f 144 $DATA > ${PREFIX}_rev1.bam
samtools index ${PREFIX}_rev1.bam

samtools view -b -f 64 -F 16 $DATA > ${PREFIX}_rev2.bam
samtools index ${PREFIX}_rev2.bam

#
# Combine alignments that originate on the reverse strand.
#
samtools merge -f ${PREFIX}_rev.bam ${PREFIX}_rev1.bam ${PREFIX}_rev2.bam
rm ${PREFIX}_rev1.bam ${PREFIX}_rev2.bam
samtools index ${PREFIX}_rev.bam
samtools sort ${PREFIX}_rev.bam -o ${PREFIX}_rev_sorted.bam
rm ${PREFIX}_rev.bam
