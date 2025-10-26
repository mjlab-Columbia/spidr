#!/bin/sh
#
# Script for running SPIDR pipeline as a set of SLURM batch jobs
#
#
#SBATCH --account=mjlab
#SBATCH --job-name=spidr
#SBATCH --output=spidr.log
#SBATCH -c 1
#SBATCH --time=0-12:00           # Time format: D-HH:MM
#SBATCH --mem=8gb

# Get the major version to determine whether --cluster-config or --profile method should be used
SNAKEMAKE_MAJOR=$(snakemake --version | tr "." "\t" | cut -f 1)

DRY_RUN=""
if [[ "$1" == "--dry_run" ]]; then
    DRY_RUN="-n"
fi

# If we're on Snakemake 8 or later, use profiles, otherwise use cluster config
# Note: 
# --cluster-config will eventually be deprecated completely so it will (eventually) be completely removed
# The option to use it is included for backwards compatibility. Occasionally, dependency conflicts between
# bioconda and conda-forge cause older versions of snakemake to be installed.
if [[ ${SNAKEMAKE_MAJOR} -lt 8 ]]; then
    snakemake \
        --use-conda \
        --rerun-incomplete \
        --keep-going \
        -j 128 \
        --cluster-config cluster.yaml \
        --printshellcmds \
        --configfile config.yaml \
        --cluster "sbatch -A {cluster.account} -c {cluster.cpus} -t {cluster.time} -N {cluster.nodes} --mem {cluster.mem} --output {cluster.output} --error {cluster.error}" \
        --cluster-cancel scancel \
        ${DRY_RUN} \
        2>&1 | tee spidr.log
elif [[ ${SNAKEMAKE_MAJOR} -ge 8 ]]; then
    snakemake \
        --profile profiles/linux_hpc \
        --configfile config.yaml \
        --omit-from generate_cluster_ecdfs \
        ${DRY_RUN} \
        2>&1 | tee spidr.log
fi

