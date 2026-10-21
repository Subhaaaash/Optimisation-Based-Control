#!/usr/bin/env bash
# Rebuild the PDF copies of every exercise's Markdown notes.
# Requires: pandoc and a TeX distribution with xelatex.
# Usage (from the repository root):  bash tools/build_pdfs.sh
set -euo pipefail
for dir in Exercise-*/; do
  mkdir -p "$dir/pdf"
  for md in "$dir"0*_*.md; do
    [ -e "$md" ] || continue
    out="$dir/pdf/$(basename "${md%.md}").pdf"
    # drop the first H1 (the PDF title block already shows it)
    sed '0,/^# /{/^# /d}' "$md" | pandoc -f markdown -o "$out" \
      --resource-path="$dir" --pdf-engine=xelatex \
      -V geometry:margin=2.2cm -V mainfont="DejaVu Sans" -V monofont="DejaVu Sans Mono" \
      -V fontsize=10pt -V colorlinks=true -V linkcolor=blue -V urlcolor=blue
    echo "built $out"
  done
done
