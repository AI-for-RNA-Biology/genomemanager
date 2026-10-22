<!--
#
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
#
# SPDX-FileContributor: Lisa Fournier <lisa.fournier@idiap.ch>
#
# SPDX-License-Identifier: GPL-3.0-only
#
-->

# SETUP

Procedure to be followed on any new project:

1. Go to <https://gitlab.idiap.ch/cverzat/template>
2. Fork the project to your namespace
3. Rename your fork (Settings/General) to your project name
4. Change URL (Settings/General/Advanced/Change path) to the same as your project name
5. Remove the fork relationship (Settings/General/Advanced/Remove fork relationship)
6. Clone your new project

   ```bash
   git clone git@gitlab.idiap.ch:{your_namespace}/{project_name}.git
   ```

7. Rename the [package](package) folder with your package name
8. Change the package name and other info (author, license, dependencies, etc) in [pyproject.toml](pyproject.toml). You could also start a new `pyproject.toml` from scratch using `poetry init`.
9. Update the [README.md](README.md)
10. Commit/push your changes

The proposed structure is always changing and comes from my experience with various projects at Idiap.

Don't hesitate to contact me on Mattermost (@cverzat) or by email (<cverzat@idiap.ch>) if you have questions.
