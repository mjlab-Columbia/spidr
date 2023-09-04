#!/bin/sh

# TODO: Add logic to allow switching between running each pipeline section or both sections sequentially

# Run pipeline for fastq --> target bam files
snakemake \
    --snakefile generate-targets.smk \
    --use-conda \
    --rerun-incomplete \
    -j 128 \
    --cluster-config cluster.generate-targets.yaml \
    --configfile config.generate-targets.yaml \
    --cluster "sbatch -A {cluster.account} -c {cluster.cpus} -t {cluster.time} -N {cluster.nodes} --mem {cluster.mem} --output {cluster.output} --error {cluster.error}" \

# TODO: Write run command for target bam --> postprocessing outputs
