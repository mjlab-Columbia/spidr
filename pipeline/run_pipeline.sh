#!/bin/sh
#
# Script for running SPIDR pipeline as a batch job
#
#
#SBATCH --account=mjlab
#SBATCH --job-name=spidr
#SBATCH -c 4
#SBATCH --time=0-12:00           # Time format: D-HH:MM
#SBATCH --mem-per-cpu=16G

# Run the whole pipeline w/ <= 128 jobs concurrently and retry multiqc if it fails
snakemake \
--snakefile Snakefile \
--use-conda \
--rerun-incomplete \
-j 128 \
--cluster-config cluster.yaml \
--configfile config.yaml \
--cluster "sbatch -A {cluster.account} -c {cluster.cpus} -t {cluster.time} -N {cluster.nodes} --mem {cluster.mem} --output {cluster.output} --error {cluster.error}" \
