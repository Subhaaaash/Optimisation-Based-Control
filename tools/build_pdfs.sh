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
    # drop the first H1 (the title block shows it) and replace GIF animations, which
    # LaTeX cannot embed, by a pointer to the file in the repository
    sed -e '0,/^# /{/^# /d}' \
        -e 's|^!\[\([^]]*\)\](\([^)]*\.gif\))$|*Animation "\1": open `\2` in the repository.*|' "$md" | pandoc -f markdown-implicit_figures -o "$out" \
      --resource-path="$dir" --pdf-engine=xelatex \
      -V geometry:margin=2.2cm -V mainfont="DejaVu Sans" -V monofont="DejaVu Sans Mono" \
      -V fontsize=10pt -V colorlinks=true -V linkcolor=blue -V urlcolor=blue
    echo "built $out"
  done
done
