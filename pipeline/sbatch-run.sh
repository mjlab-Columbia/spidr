#!/bin/sh
#
# Script for running SPIDR pipeline as a SLURM batch job
#
#
#SBATCH --account=mjlab
#SBATCH --job-name=spidr
#SBATCH -c 4
#SBATCH --time=0-12:00           # Time format: D-HH:MM
#SBATCH --mem-per-cpu=5G

# TODO: Add logic to allow switching between running each pipeline section or both sections sequentially
# Run the whole pipeline w/ <= 128 jobs concurrently 
snakemake \
    --snakefile generate-targets.smk \
    --use-conda \
    --rerun-incomplete \
    -j 128 \
    --cluster-config cluster.yaml \
    --configfile config.yaml \
    --cluster "sbatch -A {cluster.account} -c {cluster.cpus} -t {cluster.time} -N {cluster.nodes} --mem {cluster.mem} --output {cluster.output} --error {cluster.error}" \

# TODO: Write run command for target bam --> postprocessing outputs
