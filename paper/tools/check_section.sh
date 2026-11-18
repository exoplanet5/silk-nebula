#!/bin/bash
# Compile one section file standalone (AASTeX 6.31, twocolumn) to catch LaTeX errors.
# Usage: paper/tools/check_section.sh sections/s2_velocity.tex
# Undefined citations/references are expected and ignored; errors are printed.
set -u
PAPER=/Users/mickey/astronomy/silk-nebula/paper
SEC="$1"
TMP=$(mktemp -d "${TMPDIR:-/tmp}/silkchk.XXXXXX")
cat > "$TMP/t.tex" <<TEX
\documentclass[twocolumn]{aastex631}
\input{$PAPER/macros.tex}
\graphicspath{{$PAPER/}}
\begin{document}
\input{$PAPER/$SEC}
\end{document}
TEX
cd "$TMP" && pdflatex -interaction=nonstopmode -halt-on-error t.tex > log.txt 2>&1
rc=$?
grep -n -A3 '^!' t.log | head -40
grep -E 'Overfull \\hbox \(([2-9][0-9]|[0-9]{3,})' t.log | head -10
echo "exit=$rc  pages=$(grep -o 'Output written on t.pdf ([0-9]* page' t.log | grep -o '[0-9]* page')  pdf=$TMP/t.pdf"
