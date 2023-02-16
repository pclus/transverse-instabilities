#!/bin/bash
DIR="$(pwd)"
FIGTEX=$1.tex
FIGPDF=$1.pdf
cp $FIGTEX ~/TeX_Generator/temporalTeX_file.tex
cd ~/TeX_Generator/
pdflatex -synctex=1 -interaction=nonstopmode  --shell-escape temporalTeX_file.tex
cd $DIR
cp ~/TeX_Generator/temporalTeX_file.pdf $FIGPDF

