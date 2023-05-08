# Environment variables
export ROOT=/burg/mjlab/projects/spidr_preprocessing_dg/test_paired_end_pipeline
export SNAKEMAKE_OUTPUT_CACHE=$ROOT/cache

snakemake \
--snakefile Snakefile \
--use-conda \
--rerun-incomplete \
-j 128 \
--cluster-config cluster.yaml \
--configfile config.yaml \
--cluster "sbatch -A {cluster.account} -c {cluster.cpus} -t {cluster.time} -N {cluster.nodes} --mem {cluster.mem} --output {cluster.output} --error {cluster.error}"
