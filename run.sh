#!/bin/sh
#
# Script for running SPIDR pipeline as a SLURM batch job
#
#
#SBATCH --account=mjlab
#SBATCH --job-name=spidr
#SBATCH --output=spidr.log
#SBATCH -c 1
#SBATCH --time=0-12:00           # Time format: D-HH:MM
#SBATCH --mem=8gb

snakemake \
    --use-conda \
    --rerun-incomplete \
    --keep-going \
    -j 128 \
    --cluster-config cluster.yaml \
    --configfile config.yaml \
    --cluster "sbatch -A {cluster.account} -c {cluster.cpus} -t {cluster.time} -N {cluster.nodes} --mem {cluster.mem} --output {cluster.output} --error {cluster.error}" \
    --cluster-cancel scancel \
    --omit-from generate_cluster_ecdfs \
    2>&1 | tee spidr.log
