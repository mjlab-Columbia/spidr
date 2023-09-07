#!/bin/sh

# Default value is "both" which runs both pipeline sections sequentially
pipeline_to_run="both"

# Check if a command-line argument is provided
if [ $# -gt 0  ]; then
        pipeline_to_run=$1
fi

# Run the selected pipeline section(s)
if [ "$pipeline_to_run" = "generate"  ] || [ "$pipeline_to_run" = "both"  ]; then
    snakemake \
    --snakefile generate-targets.smk \
        --use-conda \
        --rerun-incomplete \
        -j 128 \
        --cluster-config cluster.generate-targets.yaml \
        --configfile config.generate-targets.yaml \
        --cluster "sbatch -A {cluster.account} -c {cluster.cpus} -t {cluster.time} -N {cluster.nodes} --mem {cluster.mem} --output {cluster.output} --error {cluster.error}" \
        --cluster-cancel scancel \
        2>&1 | tee run-log-generate-targets.txt 
fi

if [ "$pipeline_to_run" = "postprocess"  ] || [ "$pipeline_to_run" = "both"  ]; then
   snakemake \
       --snakefile postprocess-targets.smk \
       --use-conda \
       --rerun-incomplete \
       -j 128 \
       --cluster-config cluster.postprocess-targets.yaml \
       --configfile config.postprocess-targets.yaml \
       --cluster "sbatch -A {cluster.account} -c {cluster.cpus} -t {cluster.time} -N {cluster.nodes} --mem {cluster.mem} --output {cluster.output} --error {cluster.error}" \
       --cluster-cancel scancel \
       2>&1 | tee run-log-postprocess-targets.txt
fi

# Inform the user if an invalid argument is provided
if [ "$pipeline_to_run" != "generate"  ] && [ "$pipeline_to_run" != "postprocess"  ] && [ "$pipeline_to_run" != "both"  ]; then
   echo "Invalid argument. Usage: $0 [generate|postprocess|both]"
fi

