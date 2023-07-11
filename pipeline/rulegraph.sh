#!/bin/bash

snakemake --rulegraph Snakefile > rulegraph.tmp
tail -n +10 rulegraph.tmp > rulegraph.txt
dot -Tpdf rulegraph.txt > rulegraph.pdf
rm rulegraph.tmp rulegraph.txt
echo "File saved to rulegraph.pdf"