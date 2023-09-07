#!/bin/bash

echo "Generating rulegraph for fastq files --> target bam files"
snakemake \
    --rulegraph \
    --snakefile "generate-targets.smk" \
    --configfile "config.generate-targets.yaml" \
    > rulegraph-generate-targets.smk.tmp
tail -n +10 rulegraph-generate-targets.smk.tmp > rulegraph-generate-targets.smk.txt
dot -Gsize=9,16\! -Gdpi=100 -Tpdf rulegraph-generate-targets.smk.txt -o rulegraph-generate-targets.pdf
rm rulegraph-generate-targets.smk.tmp rulegraph-generate-targets.smk.txt
echo "File saved to rulegraph-generate-targets.pdf"

echo "Generating rulegraph for bam targets --> bedgraph, bigwig, and merged bam file"
snakemake \
    --rulegraph \
    --snakefile "postprocess-targets.smk" \
    --configfile "config.postprocess-targets.yaml" \
    | dot -Tpdf > rulegraph-postprocess-targets.pdf
echo "File saved to rulegraph-postprocess-targets.pdf"

