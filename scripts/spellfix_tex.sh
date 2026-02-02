#!/usr/bin/env bash
set -euo pipefail

REPO="/Volumes/Seagate/Lattie-Boltzmann-Method-Note_of_PaperReView"
cd "$REPO"

# Auto-fix common English typos in .tex files.
# codespell is conservative and only fixes known misspellings.
find . -type f -name "*.tex" -print0 | \
  xargs -0 python3 -m codespell \
    --write-changes \
    --quiet-level=2 \
    --check-filenames
