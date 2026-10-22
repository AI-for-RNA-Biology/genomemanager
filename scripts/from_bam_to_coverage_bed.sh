#!/bin/bash
#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
PATH=$PATH:/path/to/IGV:/path/to/bedtools2

DATA=$1
PREFIX=$2
CHROM_SIZES=$3

igvtools count ${DATA} ${PREFIX}.tdf $CHROM_SIZES
igvtools tdftobedgraph ${PREFIX}.tdf ${PREFIX}.bed
bedtools sort -i ${PREFIX}.bed -g $CHROM_SIZES > ${PREFIX}_sorted.bedgraph
