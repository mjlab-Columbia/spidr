#!/bin/bash

set -e # Exit on first error

# Default value is "both" which runs both pipeline sections sequentially
pipeline_to_run="both"

# Check if a command-line argument is provided
if [ $# -gt 0  ]; then
    pipeline_to_run=$1
fi

# Run the selected pipeline section(s)
if [ "$pipeline_to_run" = "generate"  ] || [ "$pipeline_to_run" = "both"  ]; then
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
fi

if [ "$pipeline_to_run" = "postprocess"  ] || [ "$pipeline_to_run" = "both"  ]; then
    if [ -f "workup/splitbams" ]; then
        echo "Generating rulegraph for bam targets --> bedgraph, bigwig, and merged bam file"
        snakemake \
            --rulegraph \
            --snakefile "postprocess-targets.smk" \
            --configfile "config.postprocess-targets.yaml" \
            | dot -Tpdf > rulegraph-postprocess-targets.pdf
        echo "File saved to rulegraph-postprocess-targets.pdf"
    else
        echo "workup/splitbams does not exist, skipping postprocess rulegraph generation."
    fi
fi

# Inform the user if an invalid argument is provided
if [ "$pipeline_to_run" != "generate"  ] && [ "$pipeline_to_run" != "postprocess"  ] && [ "$pipeline_to_run" != "both"  ]; then
   echo "Invalid argument. Usage: $0 [generate|postprocess|both]"
fi

