#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
import sys

from genomemanager.overlaps import get_overlap_file

if __name__ == "__main__":
    # Compute overlaps
    get_overlap_file(sys.argv[1], sys.argv[2])
