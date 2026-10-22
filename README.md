<!--
#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileCopyrightText: Copyright © 2025-2026 University of Bern
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@unibe.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
-->

# Project

This project aims at creating a package with useful functions to deal with annotation files.

## Description of the scripts
The script `create_focused_filtered_annotation.py` creates a focused annotation, without overlaps from a NON-FOCUSED annotation. For example, `annotation.gtf` gives `annotation_300_focus_<suffix>.gtf` (last 300 nt of every isoform) and `annotation_300_focus_no_overlaps_<suffix>.gtf`. The overlaps are removed by starting from the shortest one.

The script `create_ground_truth_no_overlaps.py` creates a annotation from an annotation already focused but with overlaps. The overlaps will be removed based on RefSeq and Ensembl annotation; and in case of overlapping known annotation, the Ensembl one will be kept. In case of overlapping Ensembl annotations, the longest one will be kept. The remaining overlaps are removed by starting from the longest one.

## Origin

genomemanager was first developed at the [Idiap Research Institute](https://www.idiap.ch), where [version 0.1](https://github.com/idiap/genomemanager) was released. It has been developed since at the University of Bern: this repository (version 0.2.0 and later) continues it.

## Installation of this package in another project

To install this library in another Python project, execute simply:

```bash
pip install git+https://github.com/AI-for-RNA-Biology/genomemanager.git
```
