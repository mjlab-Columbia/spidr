#!/bin/bash

snakemake --dag > dag.tmp
tail -n +10 dag.tmp > dag.txt
dot -Tpdf dag.txt > dag.pdf
echo "File saved to dag.pdf"