from json import load
from os import path, makedirs
from sys import exit
import datetime

configfile: "config.yaml"

required_config_keys = [
    'email',
    'bID',
    'experiments',
    'num_tags',
    'conditions',
    'assembly',
    'output_dir',
    'temp_dir',
    'num_chunks',
    'rounds_format',
    'min_oligos',
    'proportion',
    'max_size',
    'bowtie2_index',
    'star_index',
    'read1_format',
    'read2_format',
    'read1_start_offset',
    'read2_start_offset',
]

optional_config_keys = [
    'cutadapt_oligos'
]

# Check that all the required yaml fields are present
try:
    assert all([i in config.keys() for i in required_config_keys])
    out_dir = config['output_dir']
except AssertionError:
    for key in required_config_keys:
        if key not in config.keys():
            print(f"Missing: {key} in config.yaml")
    
    exit()

# For optional fields, augment the value or fall back to defaults
if 'cutadapt_oligos' not in config.keys():
    config['cutadapt_oligos'] = "-g GGTGGTCTTT -g GCCTCTTGTT"
else:
    config['cutadapt_oligos'] = "-g file:" + config['cutadapt_oligos']

# Create cluster subdirectory within logs/ directory manually
makedirs(path.join(out_dir, "workup", "logs", "cluster"), exist_ok=True)

# Create directory for benchmark tsv files to be stored
makedirs(path.join("benchmarks"), exist_ok=True)

################################################################################
#Get experiment files
###############################################################################

# Prep experiments from fastq directory using fastq2json_updated.py, now load json file
FILES = load(open(config['experiments']))
ALL_EXPERIMENTS = sorted(FILES.keys())

ALL_FASTQ = []
for SAMPLE, file in FILES.items():
    ALL_FASTQ.extend([path.abspath(i) for i in file.get('R1')])
    ALL_FASTQ.extend([path.abspath(i) for i in file.get('R2')])

NUM_CHUNKS = [f"{i:03}" for i in range(0, config['num_chunks'])]

OUTPUTS = expand(
    [
        path.join(out_dir, "workup", "count_fully_barcoded_reads", "{experiment}_R1.part_{splitid}.pre_alignment_barcode_count.txt"),
        path.join(out_dir, "workup", "count_fully_barcoded_reads", "{experiment}_R2.part_{splitid}.pre_alignment_barcode_count.txt"),
        path.join(out_dir, "workup", "plot_cdna_length_histogram", "{experiment}.cdna_histogram.pdf"),
        path.join(out_dir, "workup", "splitbams_all_conditions", "{experiment}.done"),
        path.join(out_dir, "workup", "splitbams_by_condition", "{experiment}.{condition}.done"),
        path.join(out_dir, "workup", "cat_ligation_efficiency", "ligation_efficiency.txt"),
        path.join(out_dir, "workup", "generate_cluster_statistics", "cluster_statistics.txt"),
        path.join(out_dir, "workup", "condition-clusters", "RPM_read_distribution.pdf"),
        path.join(out_dir, "workup", "condition-clusters", "RPM_cluster_distribution.pdf"),
        path.join(out_dir, "workup", "condition-clusters", "BPM_read_distribution.pdf"),
        path.join(out_dir, "workup", "condition-clusters", "BPM_cluster_distribution.pdf"),
        path.join(out_dir, "workup", "split_incorrect_clusters", "RPM_read_distribution.pdf"),       
        path.join(out_dir, "workup", "split_incorrect_clusters", "RPM_cluster_distribution.pdf"),       
        path.join(out_dir, "workup", "split_incorrect_clusters", "BPM_read_distribution.pdf"),
        path.join(out_dir, "workup", "split_incorrect_clusters", "BPM_cluster_distribution.pdf"),
        path.join(out_dir, "workup", "count_barcoded_reads_post_alignment", "{experiment}.post_alignment_barcoded_count.txt"),
        path.join(out_dir, "workup", "count_barcoded_reads_in_clusters", "{experiment}.barcoded_reads_assigned_to_clusters.txt"),
        path.join(out_dir, "workup", "count_barcoded_reads_in_bams", "{experiment}.barcoded_reads_assigned_to_bams.txt"),
        path.join(out_dir, "workup", "generate_cluster_ecdfs", "Max_representation_ecdf.pdf"),
        path.join(out_dir, "workup", "generate_cluster_ecdfs", "Max_representation_counts.pdf"),
        path.join(out_dir, "workup", "qc", "{experiment}.bowtie2_qc.log"),
        path.join(out_dir, "workup", "qc", "{experiment}.part_{splitid}.barcode_table.tsv.gz"),
        path.join(out_dir, "workup", "qc", "{experiment}.thresh_and_split_condition.{condition}.log"),
        path.join(out_dir, "workup", "qc", "{experiment}.thresh_and_split_no_condition.ALL_CONDITIONS.log")
    ],
    experiment = ALL_EXPERIMENTS,
    condition = config['conditions'],
    splitid = NUM_CHUNKS
)

rule all:
    input: 
        OUTPUTS

#Send and email if an error occurs during execution
onerror:
    shell('mail -s "an error occurred" ' + config['email'] + ' < {log}')


wildcard_constraints:
    experiment = "[^\.]+"


rule split_fastq_read1:
    input:
        r1 = lambda wildcards: FILES[wildcards.experiment]['R1']
    output:
        expand(
            path.join(out_dir, "workup", "split_fastq_read1", "{{experiment}}_R1.part_{splitid}.fastq.gz"),
            splitid=NUM_CHUNKS
        )
    params:
        output_format = " --output ".join(
            expand(
                path.join(out_dir, "workup", "split_fastq_read1", "{{experiment}}_R1.part_{splitid}.fastq.gz"),
                splitid=NUM_CHUNKS
            )
        )
    log:
        path.join(out_dir, "workup", "logs", "split_fastq_read1", "{experiment}.log")
    conda:
        "envs/fastqsplitter.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 8,
        mem_mb = 64000,
        time = "02:00:00"
    benchmark:
        "benchmarks/{experiment}.split_fastq_read1.tsv"
    shell:
        '''
        (fastqsplitter \
            --input {input.r1} \
            --output {params.output_format}) &> {log}
        '''


rule split_fastq_read2:
    input:
        r2 = lambda wildcards: FILES[wildcards.experiment]['R2']
    output:
        expand(
            path.join(out_dir, "workup", "split_fastq_read2", "{{experiment}}_R2.part_{splitid}.fastq.gz"),
            splitid=NUM_CHUNKS
        )
    params:
      output_format = " --output ".join(
        expand(
            path.join(out_dir, "workup", "split_fastq_read2", "{{experiment}}_R2.part_{splitid}.fastq.gz"),
            splitid=NUM_CHUNKS
        )
      )
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.split_fastq_read2.log")
    conda:
        "envs/fastqsplitter.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 8,
        mem_mb = 64000,
        time = "02:00:00"
    benchmark:
        "benchmarks/{experiment}.split_fastq_read2.tsv"
    shell:
        '''
        (fastqsplitter \
            --input {input.r2} \
            --output {params.output_format}) &> {log}
        '''


rule trim_sequencing_adapters:
    input:
        [path.join(out_dir, "workup", "split_fastq_read1", "{experiment}_R1.part_{splitid}.fastq.gz"), 
         path.join(out_dir, "workup", "split_fastq_read2", "{experiment}_R2.part_{splitid}.fastq.gz")]
    output:
        path.join(out_dir, "workup", "trim_sequencing_adapters", "{experiment}_R1.part_{splitid}_val_1.fq.gz"),
        path.join(out_dir, "workup", "trim_sequencing_adapters", "{experiment}_R1.part_{splitid}.fastq.gz_trimming_report.txt"),
        path.join(out_dir, "workup", "trim_sequencing_adapters", "{experiment}_R2.part_{splitid}_val_2.fq.gz"),
        path.join(out_dir, "workup", "trim_sequencing_adapters", "{experiment}_R2.part_{splitid}.fastq.gz_trimming_report.txt")
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.trim_sequencing_adapters.log")
    conda:
        "envs/trim_galore.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 15,
        mem_mb = 50000,
        time = "12:00:00"
    benchmark:
        "benchmarks/{experiment}.{splitid}.trim_sequencing_adapters.tsv"
    shell:
        '''
        (trim_galore \
            --paired \
            --gzip \
            --cores 4 \
            --quality 20 \
            --fastqc \
            -o {out_dir}workup/trim_sequencing_adapters/ \
            {input}) &> {log}
        '''


rule identify_barcodes:
    input:
        r1 = path.join(out_dir, "workup", "trim_sequencing_adapters", "{experiment}_R1.part_{splitid}_val_1.fq.gz"),
        r2 = path.join(out_dir, "workup", "trim_sequencing_adapters", "{experiment}_R2.part_{splitid}_val_2.fq.gz")
    output:
        r1_barcoded = path.join(out_dir, "workup", "identify_barcodes", "{experiment}_R1.part_{splitid}.barcoded.fastq.gz"),
        r2_barcoded = path.join(out_dir, "workup", "identify_barcodes", "{experiment}_R2.part_{splitid}.barcoded.fastq.gz")
    params:
        bid_config = config['bID'],
        read1_format = config["read1_format"],
        read2_format = config["read2_format"],
        read1_start_offset = config["read1_start_offset"],
        read2_start_offset = config["read2_start_offset"],
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.identify_barcodes.log")
    conda:
        "envs/barcoding.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 64000,
        time = "12:00:00"
    benchmark:
        "benchmarks/{experiment}.{splitid}.identify_barcodes.tsv"
    shell:
        """
        (python scripts/python/identify_barcodes.py \
            --input_read1 {input.r1} \
            --input_read2 {input.r2} \
            --output_read1 {output.r1_barcoded} \
            --output_read2 {output.r2_barcoded} \
            --read1_format '{params.read1_format}' \
            --read2_format '{params.read2_format}' \
            --read1_start_offset {params.read1_start_offset} \
            --read2_start_offset {params.read2_start_offset} \
            --config {params.bid_config}) &> {log}
        """


rule count_barcoded_reads_pre_alignment:
    """
    Count the number of reads with a full combinatorial barcode after initial barcode identification

    Explanation of shell logic:
    gzip - Decompresses fastq.gz and streams it to stdout
    awk - extract header lines where barcodes are stored
    grep - extract barcode from header @NB551203:643:HNFCCAFX3:1:11101:21927:1057::[BEAD_ILF3][NYBot67_Stg][ROUND5_H4][ROUND4_E6][ROUND3_E5][ROUND2_C3][ROUND1_CNTRL_A1] --> [BEAD_ILF3][NYBot67_Stg][ROUND5_H4][ROUND4_E6][ROUND3_E5][ROUND2_C3][ROUND1_CNTRL_A1]
    sed - replace ][ with tabs and ends with empty string to make each line tab-delimited
    cut - extract every column except the first one
    wc - count number of lines
    """
    input:
        r1_barcoded = path.join(out_dir, "workup", "identify_barcodes", "{experiment}_R1.part_{splitid}.barcoded.fastq.gz"),
        r2_barcoded = path.join(out_dir, "workup", "identify_barcodes", "{experiment}_R2.part_{splitid}.barcoded.fastq.gz")
    output:
        r1_count = path.join(out_dir, "workup", "count_fully_barcoded_reads", "{experiment}_R1.part_{splitid}.pre_alignment_barcode_count.txt"),
        r2_count = path.join(out_dir, "workup", "count_fully_barcoded_reads", "{experiment}_R2.part_{splitid}.pre_alignment_barcode_count.txt")
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.count_barcoded_reads_pre_alignment.log")
    conda:
        "envs/coreutils.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 64000,
        time = "12:00:00"
    benchmark:
        "benchmarks/{experiment}.{splitid}.count_barcoded_reads_pre_alignment.tsv"
    shell:
        """
        (gzip -dc {input.r1_barcoded} \
            | awk 'NR % 4 == 1' \
            | grep -oP '::\K(\[[^\]]+\])+' \
            | sed -E 's/\]\[/\t/g; s/\[|\]//g' \
            | cut -f 1 --complement \
            | grep -v 'NOT_FOUND' \
            | wc -l > {output.r1_count}) &> {log}

        (gzip -dc {input.r2_barcoded} \
            | awk 'NR % 4 == 1' \
            | grep -oP '::\K(\[[^\]]+\])+' \
            | sed -E 's/\]\[/\t/g; s/\[|\]//g' \
            | cut -f 1 --complement \
            | grep -v 'NOT_FOUND' \
            | wc -l > {output.r2_count}) &>> {log}
        """


rule generate_barcode_table:
    input:
        path.join(out_dir, "workup", "identify_barcodes", "{experiment}_R1.part_{splitid}.barcoded.fastq.gz"),
    output:
        table = path.join(out_dir, "workup", "qc", "{experiment}.part_{splitid}.barcode_table.tsv.gz"),
        barcode_qc = path.join(out_dir, "workup", "qc", "{experiment}.part_{splitid}.barcode_qc.txt")
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.generate_barcode_table.log")
    conda:
        "envs/barcoding.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 16000,
        time = "01:00:00"
    benchmark:
        "benchmarks/{experiment}.{splitid}.generate_barcode_table.tsv"
    shell:
        """
        (gzip -dc {input} \
            | awk 'NR%4==1' \
            | grep -Eo '(\[[^]]*\])+$' \
            | sed 's/\]\[/\t/g; s/^\[//; s/\]$//'\
            | gzip > {output.table}) &> {log}

        (echo 'Number of fully barcoded reads: ' \
            $(gzip -dc {output.table} \
                | cut -f 1 --complement \
                | grep -v 'NOT_FOUND' \
                | wc -l) > {output.barcode_qc}) &>> {log}

        (echo 'Number of incompletely barcoded reads: ' \
            $(gzip -dc {output.table} \
                | cut -f 1 --complement \
                | grep 'NOT_FOUND' \
                | wc -l) > {output.barcode_qc}) &>> {log}
        """


rule get_ligation_efficiency:
    input:
        r1 = path.join(out_dir, "workup", "identify_barcodes", "{experiment}_R1.part_{splitid}.barcoded.fastq.gz")
    output:
        path.join(out_dir, "workup",  "get_ligation_efficiency", "{experiment}.part_{splitid}.ligation_efficiency.txt")
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.get_ligation_efficiency.log")
    conda:
        "envs/sprite.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 8000,
        time = "01:00:00"
    shell:
        """
        (python scripts/python/get_ligation_efficiency.py {input.r1} > {output}) &> {log}
        """


rule cat_ligation_efficiency:
    input:
        expand(
            path.join(out_dir, "workup",  "get_ligation_efficiency", "{experiment}.part_{splitid}.ligation_efficiency.txt"), 
            experiment=ALL_EXPERIMENTS, 
            splitid=NUM_CHUNKS
        )
    output:
        path.join(out_dir, "workup", "cat_ligation_efficiency", "ligation_efficiency.txt")
    log:
        path.join(out_dir, "workup", "logs", "cat_ligation_efficiency.log")
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 20000,
        time = "05:00:00"
    shell:
        """
        (tail -n +1 {input} > {output}) &> {log}
        """


rule split_reads_read1:
    '''
    split bpm and rpm will also remove incomplete barcodes
    '''
    input:
        path.join(out_dir, "workup", "identify_barcodes", "{experiment}_R1.part_{splitid}.barcoded.fastq.gz")
    output:
        rpm = path.join(out_dir, "workup", "split_reads_read1", "{experiment}_R1.part_{splitid}.barcoded_rpm.fastq.gz"),
        bpm = path.join(out_dir, "workup", "split_reads_read1", "{experiment}_R1.part_{splitid}.barcoded_bpm.fastq.gz"),
        short = path.join(out_dir, "workup", "split_reads_read1", "{experiment}_R1.part_{splitid}.barcoded_short.fastq.gz")
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.split_reads_read1.log")
    conda:
        "envs/python.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 48000,
        time = "04:00:00"
    benchmark:
        "benchmarks/{experiment}.{splitid}.split_reads_read1.tsv"
    shell:
        """
        (python scripts/python/split_rpm_bpm_fq.py \
            --input {input} \
            --rpm_output {output.rpm} \
            --bpm_output {output.bpm} \
            --short_output {output.short}) &> {log}
        """


rule split_reads_read2:
    '''
    split bpm and rpm will also remove incomplete barcodes
    '''
    input:
        path.join(out_dir, "workup", "identify_barcodes", "{experiment}_R2.part_{splitid}.barcoded.fastq.gz")
    output:
        rpm = path.join(out_dir, "workup", "split_reads_read2", "{experiment}_R2.part_{splitid}.barcoded_rpm.fastq.gz"),
        bpm = path.join(out_dir, "workup", "split_reads_read2", "{experiment}_R2.part_{splitid}.barcoded_bpm.fastq.gz"),
        short = path.join(out_dir, "workup", "split_reads_read2", "{experiment}_R2.part_{splitid}.barcoded_short.fastq.gz")
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.split_reads_read2.log")
    conda:
        "envs/python.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 48000,
        time = "04:00:00"
    benchmark:
        "benchmarks/{experiment}.{splitid}.split_reads_read2.tsv"
    shell:
        """
        (python scripts/python/split_rpm_bpm_fq.py \
            --input {input} \
            --rpm_output {output.rpm} \
            --bpm_output {output.bpm} \
            --short_output {output.short}) &> {log}
        """


rule trim_rpm_reads:
    input:
        read1 = path.join(out_dir, "workup", "split_reads_read1", "{experiment}_R1.part_{splitid}.barcoded_rpm.fastq.gz"),
        read2 = path.join(out_dir, "workup", "split_reads_read2", "{experiment}_R2.part_{splitid}.barcoded_rpm.fastq.gz")
    output:
        r1 = path.join(out_dir, "workup", "trim_rpm_reads", "{experiment}_R1.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz"),
        r2 = path.join(out_dir, "workup", "trim_rpm_reads", "{experiment}_R2.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz"),
        qc = path.join(out_dir, "workup", "trim_rpm_reads", "{experiment}_R1.part_{splitid}.barcoded_rpm.RDtrim.qc.txt")
    params:
        adapters_r1 = "-a ATCAGCACTTAGCGTCAG",
        adapters_r2 = "-G CTGACGCTAAGTGCTGAT",
        others = "--minimum-length 20"
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.trim_rpm_reads.log")
    conda:
        "envs/cutadapt.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 48000,
        time = "12:00:00"
    benchmark:
        "benchmarks/{experiment}.{splitid}.trim_rpm_reads.tsv"
    shell:
        '''
        (cutadapt \
            {params.adapters_r1} \
            {params.adapters_r2} \
            {params.others} \
            -o {output.r1} \
            -p {output.r2} \
            -j 0 \
            {input.read1} {input.read2} > {output.qc}) &> {log}
        
        (fastqc {output.r1}) &>> {log}
        fastqc {output.r2} &>> {log}
        '''


rule calculate_cdna_length:
    """
    Take the cDNA reads which have been stripped of primers and other stuff and count the length of each cDNA molecule
    """
    input:
        fq1 = path.join(out_dir, "workup", "trim_rpm_reads", "{experiment}_R1.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz"),
        fq2 = path.join(out_dir, "workup", "trim_rpm_reads", "{experiment}_R2.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz")
    output:
        fq1 = path.join(out_dir, "workup", "calculate_cdna_length", "{experiment}_R1.part_{splitid}.txt.gz"),
        fq2 = path.join(out_dir, "workup", "calculate_cdna_length", "{experiment}_R2.part_{splitid}.txt.gz")
    conda:
        "envs/python.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 16000,
        time = "00:30:00"
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.calculate_cdna_length.log")
    benchmark:
        "benchmarks/{experiment}.part_{splitid}.calculate_cdna_length.tsv"
    shell:
        """
        (gzip -dc {input.fq1} \
            | awk 'NR % 4 == 2' \
            | awk '{{print length}}' \
            | gzip > {output.fq1}) &> {log}

        (gzip -dc {input.fq2} \
            | awk 'NR % 4 == 2' \
            | awk '{{print length}}' \
            | gzip > {output.fq2}) &>> {log}
        """


rule aggregate_cdna_lengths_across_splits:
    """
    Combine cDNA lengths calculated per chunk into a single file
    """
    input:
        fq1 = expand(
            path.join(out_dir, "workup", "calculate_cdna_length", "{{experiment}}_R1.part_{splitid}.txt.gz"),
            splitid = NUM_CHUNKS
        ),
        fq2 = expand(
            path.join(out_dir, "workup", "calculate_cdna_length", "{{experiment}}_R2.part_{splitid}.txt.gz"),
            splitid = NUM_CHUNKS
        )
    output:
        fq1 = path.join(out_dir, "workup", "aggregate_cdna_lengths_across_splits", "{experiment}.cdna_lengths_R1.txt.gz"),
        fq2 = path.join(out_dir, "workup", "aggregate_cdna_lengths_across_splits", "{experiment}.cdna_lengths_R2.txt.gz")
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.aggregate_cdna_lengths_across_splits.log")
    conda:
        "envs/coreutils.yaml"
    threads:
        1
    benchmark:
        "benchmarks/{experiment}.aggregate_cdna_lengths_across_splits.tsv"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 16000,
        time = "00:30:00"
    shell:
        """
        (zcat -dck {input.fq1} | gzip > {output.fq1}) &> {log}
        (zcat -dck {input.fq2} | gzip > {output.fq2}) &> {log}
        """


rule plot_cdna_length_histogram:
    """
    Plot a histogram of cDNA lengths based on the cDNA lengths across all chunks
    """
    input:
        fq1 = path.join(out_dir, "workup", "aggregate_cdna_lengths_across_splits", "{experiment}.cdna_lengths_R1.txt.gz"),
        fq2 = path.join(out_dir, "workup", "aggregate_cdna_lengths_across_splits", "{experiment}.cdna_lengths_R2.txt.gz")
    output:
        path.join(out_dir, "workup", "plot_cdna_length_histogram", "{experiment}.cdna_histogram.pdf")
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.aggregate_cdna_lengths_across_splits.log")
    conda:
        "envs/python.yaml"
    threads:
        1
    benchmark:
        "benchmarks/{experiment}.aggregate_cdna_lengths_across_splits.tsv"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 16000,
        time = "00:30:00"
    shell:
        """
        (python scripts/python/plot_cdna_histogram.py \
            --read1 {input.fq1} \
            --read2 {input.fq2} \
            --output {output}) &> {log}
        """


rule trim_bead_oligo_reads:
    '''
    Trim 9mer oligo sequence from bead barcode
    '''
    input:
        path.join(out_dir, "workup", "split_reads_read1", "{experiment}_R1.part_{splitid}.barcoded_bpm.fastq.gz")
    output:
        fastq = path.join(out_dir, "workup", "trim_bead_oligo_reads", "{experiment}_R1.part_{splitid}.barcoded_bpm.RDtrim.fastq.gz"),
        qc = path.join(out_dir, "workup", "trim_bead_oligo_reads", "{experiment}_R1.part_{splitid}.barcoded_bpm.RDtrim.qc.txt")
    params:
        adapters_r1 = config['cutadapt_oligos']
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.trim_bead_oligo_reads.log")
    threads: 
        10
    conda:
        "envs/cutadapt.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 48000,
        time = "12:00:00"
    benchmark:
        "benchmarks/{experiment}.{splitid}.trim_bead_oligo_reads.tsv"
    shell:
        '''
        (cutadapt \
            {params.adapters_r1} \
            -o {output.fastq} \
            -j 0 \
            {input} > {output.qc}) &> {log}
        '''


rule align_bowtie2:
    input:
        fq1 = path.join(out_dir, "workup", "trim_rpm_reads", "{experiment}_R1.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz"),
        fq2 = path.join(out_dir, "workup", "trim_rpm_reads", "{experiment}_R2.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz")
    output:
        bam = temp(path.join(out_dir, "workup", "align_bowtie2", "{experiment}.part_{splitid}.bowtie2.bam")),
        mapped = path.join(out_dir, "workup", "align_bowtie2", "{experiment}.part_{splitid}.bowtie2.sorted.mapped.bam"),
        unmapped = path.join(out_dir, "workup", "align_bowtie2", "{experiment}.part_{splitid}.bowtie2.sorted.unmapped.bam"),
        index = path.join(out_dir, "workup", "align_bowtie2", "{experiment}.part_{splitid}.bowtie2.sorted.mapped.bam.bai")
    params:
        BOWTIE2_INDEX = config['bowtie2_index'][config['assembly']]
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.align_bowtie2.log")
    threads: 
        4
    conda:
        "envs/bowtie2.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 48000,
        time = "01:00:00"
    benchmark:
        "benchmarks/{experiment}.{splitid}.align_bowtie2.tsv"
    shell:
        '''
        (bowtie2 \
            -p 4 \
            -t \
            -x {params.BOWTIE2_INDEX} \
            -1 {input.fq1} \
            -2 {input.fq2} \
            | samtools view -bS -> {output.bam}) &> {log}
        
        (samtools view -b -f 4 {output.bam} | samtools sort -n -o {output.unmapped}) &>> {log}
        (samtools view -b -F 4 {output.bam} | samtools sort > {output.mapped}) &>> {log}
        (samtools index {output.mapped}) &>> {log}
        '''


rule collate_bowtie2_qc:
    input:
        expand(
            path.join(out_dir, "workup", "logs", "{{experiment}}.{splitid}.align_bowtie2.log"),
            splitid=NUM_CHUNKS
        )
    output:
        path.join(out_dir, "workup", "qc", "{experiment}.bowtie2_qc.log")
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.collate_bowtie2_qc.log")
    conda:
        "envs/bowtie2.yaml"
    threads:
        1
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 4000,
        time = "00:30:00"
    benchmark:
        "benchmarks/{experiment}.collate_bowtie2_qc.tsv"
    shell:
        """
        (for log in {input}; do
            echo $(basename $log .align_bowtie2.log) >> {output};
            cat $log >> {output};
            echo "" >> {output};
        done) &> {log}
        """


rule convert_bam_to_fastq:
    input: 
        path.join(out_dir, "workup", "align_bowtie2", "{experiment}.part_{splitid}.bowtie2.sorted.unmapped.bam")
    output: 
        r1 = path.join(out_dir, "workup", "convert_bam_to_fastq", "{experiment}.part_{splitid}.bowtie2.unmapped_R1.fq.gz"),
        r2 = path.join(out_dir, "workup", "convert_bam_to_fastq", "{experiment}.part_{splitid}.bowtie2.unmapped_R2.fq.gz")
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.convert_bam_to_fastq.log")
    threads: 
        1
    conda:
        "envs/samtools.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 32000,
        time = "01:00:00"
    benchmark:
        "benchmarks/{experiment}.{splitid}.convert_bam_to_fastq.tsv"
    shell:
        '''
        (samtools fastq \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 /dev/null \
            -s /dev/null \
            -@ {threads} \
            -n {input}) &> {log}
        '''


rule align_star:
    input:
        r1 = path.join(out_dir, "workup", "convert_bam_to_fastq", "{experiment}.part_{splitid}.bowtie2.unmapped_R1.fq.gz"),
        r2 = path.join(out_dir, "workup", "convert_bam_to_fastq", "{experiment}.part_{splitid}.bowtie2.unmapped_R2.fq.gz")
    output:
        sorted_by_coord = path.join(out_dir, "workup", "align_star", "{experiment}.part_{splitid}.Aligned.sortedByCoord.out.bam"),
        index = path.join(out_dir, "workup", "align_star", "{experiment}.part_{splitid}.Aligned.sortedByCoord.out.bam.bai")
    params:
        STAR_OPTIONS = " ".join(["--readFilesCommand zcat",
                                "--alignEndsType EndToEnd",
                                "--outFilterScoreMin 10",
                                "--outFilterMultimapNmax 1",
                                "--outFilterMismatchNmax 10",
                                "--alignIntronMax 100000",
                                "--alignMatesGapMax 1300",
                                "--alignIntronMin 80",
                                "--alignSJDBoverhangMin 5",
                                "--alignSJoverhangMin 8",
                                "--chimSegmentMin 20",
                                "--alignSJstitchMismatchNmax 5 -1 5 5",
                                "--outSAMunmapped Within",
                                "--outReadsUnmapped Fastx",
                                "--outSAMtype BAM SortedByCoordinate"]),
        prefix = path.join(out_dir, "workup", "align_star", "{experiment}.part_{splitid}."),
        STAR_INDEX = config['star_index'][config['assembly']]
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.align_star.log")
    threads:
        10
    conda:
        "envs/star.yaml"
    benchmark:
        "benchmarks/{experiment}.{splitid}.align_star.tsv"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 10,
        mem_mb = 100000,
        time = "12:00:00"
    shell:
        '''
        mkdir -p workup/align_star

        (STAR \
            --genomeDir {params.STAR_INDEX} \
            --readFilesIn {input.r1} {input.r2} \
            --runThreadN {threads} {params.STAR_OPTIONS} \
            --outFileNamePrefix {params.prefix}) &> {log}

        (samtools index {output.sorted_by_coord}) &>> {log}
        '''


rule add_chromosome_info_bowtie2:
    input:
        bt2 = path.join(out_dir, "workup", "align_bowtie2", "{experiment}.part_{splitid}.bowtie2.sorted.mapped.bam")
    output:
        bt2 = path.join(out_dir, "workup", "add_chromosome_info_bowtie2", "{experiment}.part_{splitid}.bowtie2.sorted.mapped.chr.bam")
    params:
        assembly = config['assembly']
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.add_chromosome_info_bowtie2.log"),
    conda:
        "envs/sprite.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb= 48000,
        time = "02:00:00"
    benchmark:
        "benchmarks/{experiment}.{splitid}.add_chromosome_info_bowtie2.tsv"
    shell:
        '''
        (python scripts/python/add_chr_bt2.py \
            -i {input.bt2} \
            -o {output.bt2} \
            --assembly {params.assembly}) &> {log}
        '''


rule add_chromosome_info_star:
    input:
        star = path.join(out_dir, "workup", "align_star", "{experiment}.part_{splitid}.Aligned.sortedByCoord.out.bam"),
        index = path.join(out_dir, "workup", "align_star", "{experiment}.part_{splitid}.Aligned.sortedByCoord.out.bam.bai"),
    output:
        star = path.join(out_dir, "workup", "add_chromosome_info_star", "{experiment}.part_{splitid}.Aligned.out.sortedByCoord.chr.bam"),
    params:
        assembly = config['assembly']
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.add_chromosome_info_star.log"),
    conda:
        "envs/sprite.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 48000,
        time = "01:00:00"
    benchmark:
        "benchmarks/{experiment}.{splitid}.add_chromosome_info_star.tsv"
    shell:
        '''
        (python scripts/python/ensembl2ucsc.py \
            -i {input.star} \
            -o {output.star} \
            --assembly {params.assembly}) &> {log}
        '''


rule merge_rna_bams:
    input:
        bt2 = expand(
            path.join(out_dir, "workup", "add_chromosome_info_bowtie2", "{{experiment}}.part_{splitid}.bowtie2.sorted.mapped.chr.bam"),
            splitid=NUM_CHUNKS
        ),
        star = expand(
            path.join(out_dir, "workup", "add_chromosome_info_star", "{{experiment}}.part_{splitid}.Aligned.out.sortedByCoord.chr.bam"), 
            splitid=NUM_CHUNKS
        )
    output:
        path.join(out_dir, "workup", "merge_rna_bams", "{experiment}.merged.RPM.bam")
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.merge_rna_bams.log")
    conda:
        "envs/samtools.yaml"
    threads:
        8
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 48000,
        time = "02:00:00"
    benchmark:
        "benchmarks/{experiment}.merge_rna_bams.tsv"
    shell:
        '''
        (samtools merge -@ {threads} {output} {input.bt2} {input.star}) &> {log}
        (samtools index {output}) &>> {log}
        '''


rule count_barcoded_reads_post_alignment:
    input:
        path.join(out_dir, "workup", "merge_rna_bams", "{experiment}.merged.RPM.bam")
    output:
        path.join(out_dir, "workup", "count_barcoded_reads_post_alignment", "{experiment}.post_alignment_barcoded_count.txt")
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.count_barcoded_reads_post_alignment.log")
    conda:
        "envs/samtools.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 4,
        mem_mb = 16000,
        time = "00:30:00"
    shell:
        """
        (samtools view --no-header {input} \
            | awk '{{print $1}}' \
            | sort -u -T {resources.tmpdir} --parallel {resources.cpus} --buffer-size {resources.mem_mb}M \
            | wc -l > {output}) &> {log}
        """


rule convert_fastq_to_bam:
    input:
        path.join(out_dir, "workup", "trim_bead_oligo_reads", "{experiment}_R1.part_{splitid}.barcoded_bpm.RDtrim.fastq.gz")
    output:
        sorted = path.join(out_dir, "workup", "convert_fastq_to_bam", "{experiment}.part_{splitid}.BPM.bam"),
        bam = temp(path.join(out_dir, "workup", "convert_fastq_to_bam", "{experiment}.part_{splitid}.BPM.unsorted.bam"))
    params:
        bid_config = config['bID']
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.make_bam.log")
    conda:
        "envs/pysam.yaml"
    threads:
        8
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 48000,
        time = "12:00:00"
    benchmark:
        "benchmarks/{experiment}.{splitid}.convert_fastq_to_bam.tsv"
    shell:
        '''
        (python scripts/python/fastq_to_bam.py \
            --input {input} \
            --output {output.bam} \
            --config {params.bid_config}) &> {log}

        (samtools sort -@ {threads} -o {output.sorted} {output.bam}) &>> {log}
        '''


rule merge_bead_bams:
    input:
        expand(
            path.join(out_dir, "workup", "convert_fastq_to_bam", "{{experiment}}.part_{splitid}.BPM.bam"), 
            splitid = NUM_CHUNKS
        )
    output:
        path.join(out_dir, "workup", "merge_bead_bams", "{experiment}.merged.BPM.bam")
    conda:
        "envs/sprite.yaml"
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.merge_bead_bams.log")
    threads:
        8
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 48000,
        time = "24:00:00"
    benchmark:
        "benchmarks/{experiment}.merge_bead_bams.tsv"
    shell:
        '''
        (samtools merge -@ {threads} {output} {input}) &> {log}
        '''


rule make_clusters:
    input:
        rpm = path.join(out_dir, "workup", "add_chromosome_info_star", "{experiment}.part_{splitid}.Aligned.out.sortedByCoord.chr.bam"),
        bpm = path.join(out_dir, "workup", "convert_fastq_to_bam", "{experiment}.part_{splitid}.BPM.bam"),
        bt2 = path.join(out_dir, "workup", "add_chromosome_info_bowtie2", "{experiment}.part_{splitid}.bowtie2.sorted.mapped.chr.bam")
    output:
        unsorted = temp(path.join(out_dir, "workup", "make_clusters", "{experiment}.part_{splitid}.unsorted.clusters")),
        sorted = path.join(out_dir, "workup", "make_clusters", "{experiment}.part_{splitid}.clusters")
    params:
        temp_dir = config['temp_dir'],
        num_tags = config['num_tags']
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.make_clusters.log")
    conda:
        "envs/sprite.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 32000,
        time = "01:00:00"
    benchmark:
        "benchmarks/{experiment}.{splitid}.make_clusters.tsv"
    shell:
        '''
        (python scripts/python/get_clusters.py \
            -i {input.bpm} {input.bt2} {input.rpm}\
            -o {output.unsorted} \
            -n {params.num_tags})  &> {log}

        (sort -k 1 -T {params.temp_dir} {output.unsorted} > {output.sorted}) &> {log}
        '''


rule merge_clusters:
    input:
        expand(
            path.join(out_dir, "workup", "make_clusters", "{{experiment}}.part_{splitid}.clusters"), 
            splitid=NUM_CHUNKS
        )
    output:
        mega = temp(path.join(out_dir, "workup", "merge_clusters", "{experiment}.duplicated.clusters")),
        final = path.join(out_dir, "workup", "merge_clusters", "{experiment}.clusters")
    params:
        temp_dir = config['temp_dir']
    conda:
       "envs/sprite.yaml"
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.merge_clusters.log")
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 64000,
        time = "01:00:00"
    benchmark:
        "benchmarks/{experiment}.merge_clusters.tsv"
    shell:
        '''
        (sort -k 1 -T {params.temp_dir} -m {input} > {output.mega}) &> {log}
        (python scripts/python/merge_clusters.py -i {output.mega} -o {output.final}) &>> {log}
        '''        


rule split_incorrect_clusters:
    input:
        clusters = path.join(out_dir, "workup", "merge_clusters", "{experiment}.clusters")
    output:
        complete_clusters = path.join(out_dir, "workup", "split_incorrect_clusters", "{experiment}.complete.clusters"),
        incomplete_clusters = path.join(out_dir, "workup", "split_incorrect_clusters", "{experiment}.incomplete.clusters")
    params:
        rounds_format = config["rounds_format"]
    conda:
        "envs/sprite.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 32000,
        time = "00:15:00"
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.split_incorrect_clusters.log")
    benchmark:
        "benchmarks/{experiment}.split_incorrect_clusters.tsv"
    shell:
        '''
        (python scripts/python/split_incorrect_clusters.py \
            --clusters {input.clusters} \
            --complete_output {output.complete_clusters} \
            --incomplete_output {output.incomplete_clusters} \
            --format {params.rounds_format}) &> {log}
        '''


rule get_bpm_rpm_counts:
    input:
        path.join(out_dir, "workup", "split_incorrect_clusters", "{experiment}.complete.clusters"),
    output:
        path.join(out_dir, "workup", "get_bpm_rpm_counts", "{experiment}.bpm_rpm_counts.tsv.gz")
    conda:
        "envs/python.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 16000,
        time = "00:30:00"
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.get_bpm_rpm_counts.log")
    benchmark:
        "benchmarks/{experiment}.get_bpm_rpm_counts.tsv"
    shell:
        """
        python scripts/python/get_bpm_rpm_counts.py --clusters {input} --output {output}
        """


rule count_barcoded_reads_in_clusters:
    input:
        path.join(out_dir, "workup", "get_bpm_rpm_counts", "{experiment}.bpm_rpm_counts.tsv.gz")
    output:
        path.join(out_dir, "workup", "count_barcoded_reads_in_clusters", "{experiment}.barcoded_reads_assigned_to_clusters.txt")
    conda:
        "envs/python.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 16000,
        time = "00:30:00"
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.get_bpm_rpm_counts.log")
    benchmark:
        "benchmarks/{experiment}.count_barcoded_reads_in_clusters.tsv"
    shell:
        """
        (python scripts/python/count_barcoded_reads_in_clusters.py \
            --bpm_rpm_counts {input} \
            --output {output}) &> {log}
        """


rule generate_cluster_statistics:
    input:
        expand(
            path.join(out_dir, "workup", "split_incorrect_clusters", "{experiment}.complete.clusters"), 
            experiment = ALL_EXPERIMENTS
        )
    output:
        path.join(out_dir, "workup", "generate_cluster_statistics", "cluster_statistics.txt")
    params:
        dir = path.join(out_dir, "workup", "generate_cluster_statistics")
    conda:
        "envs/sprite.yaml"
    log:
        path.join(out_dir, "workup", "logs", "generate_cluster_statistics.log")
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 50000,
        time = "12:00:00"
    shell:
        """
        (python scripts/python/generate_cluster_statistics.py \
            --directory {params.dir} \
            --pattern .clusters > {output}) &> {log}
        """


rule generate_cluster_ecdfs:
    input:
        expand(
            path.join(out_dir, "workup", "merge_clusters", "{experiment}.clusters"), 
            experiment=ALL_EXPERIMENTS
        )
    output:
        ecdf = path.join(out_dir, "workup", "generate_cluster_ecdfs", "Max_representation_ecdf.pdf"),
        counts = path.join(out_dir, "workup", "generate_cluster_ecdfs", "Max_representation_counts.pdf")
    params:
        input_dir = path.join(out_dir, "workup", "merge_clusters"),
        output_dir = path.join(out_dir, "workup", "generate_cluster_ecdfs")
    conda:
        "envs/plotting.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 50000,
        time = "12:00:00"
    log:
        path.join(out_dir, "workup", "logs", "generate_cluster_ecdfs.log")        
    shell:
        """
        (python scripts/python/max_representation_ecdfs_perlib.py \
            --input_directory {params.input_dir} \
            --output_directory {params.output_dir} \
            --pattern .clusters \
            --xlim 30) &> {log}
        """


rule split_on_first_tag:
    input:
        complete_clusters = path.join(out_dir, "workup", "split_incorrect_clusters", "{experiment}.complete.clusters")
    output:
        expand(
            path.join(out_dir, "workup", "condition-clusters", "{{experiment}}.{condition}.clusters"),
            condition = config['conditions']
        )
    conda:
        "envs/sprite.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 32000,
        time = "00:15:00"
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.splitonfirsttag.log")
    benchmark:
        "benchmarks/{experiment}.split_on_first_tag.tsv"
    shell:
        '''
        (python scripts/python/split_on_first_tag.py \
            --complete_clusters {input.complete_clusters} \
            --output_dir workup/condition-clusters) &> {log}
        '''


# Profile size distribution of clusters
rule get_size_distribution:
    input:
        expand(
            path.join(out_dir, "workup", "condition-clusters", "{experiment}.{condition}.clusters"), 
            experiment = ALL_EXPERIMENTS, 
            condition = config['conditions']
        ),
        expand(
            path.join(out_dir, "workup", "split_incorrect_clusters", "{experiment}.complete.clusters"), 
            experiment=ALL_EXPERIMENTS
        ),
    output:
        rpm = path.join(out_dir, "workup", "condition-clusters", "RPM_read_distribution.pdf"),
        rpm2 = path.join(out_dir, "workup", "condition-clusters", "RPM_cluster_distribution.pdf"),
        bpm = path.join(out_dir, "workup", "condition-clusters", "BPM_read_distribution.pdf"),
        bpm2 = path.join(out_dir, "workup", "condition-clusters", "BPM_cluster_distribution.pdf"),
        no_condition_rpm = path.join(out_dir, "workup", "split_incorrect_clusters", "RPM_read_distribution.pdf"),       
        no_condition_rpm2 = path.join(out_dir, "workup", "split_incorrect_clusters", "RPM_cluster_distribution.pdf"),       
        no_condition_bpm = path.join(out_dir, "workup", "split_incorrect_clusters", "BPM_read_distribution.pdf"),
        no_condition_bpm2 = path.join(out_dir, "workup", "split_incorrect_clusters", "BPM_cluster_distribution.pdf")   
    params:
        condition_dir = path.join(out_dir, "workup", "condition-clusters"),
        no_condition_dir = path.join("workup", "split_incorrect_clusters")
    log:
        path.join(out_dir, "workup", "logs", "get_size_distribution.log")
    conda:
        "envs/sprite.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 50000,
        time = "12:00:00"
    shell:
        '''
        (python scripts/python/get_bead_size_distribution.py --directory {params.no_condition_dir} --pattern .clusters --readtype BPM) &> {log}
        (python scripts/python/get_bead_size_distribution.py --directory {params.no_condition_dir} --pattern .clusters --readtype RPM) &>> {log}

        (python scripts/python/get_bead_size_distribution.py --directory {params.condition_dir} --pattern .clusters --readtype BPM) &>> {log}
        (python scripts/python/get_bead_size_distribution.py --directory {params.condition_dir} --pattern .clusters --readtype RPM) &>> {log}
        '''


# Generate bam files for individual targets based on assignments from clusterfile
rule thresh_and_split_condition:
    input:
        bam = path.join(out_dir, "workup", "merge_rna_bams", "{experiment}.merged.RPM.bam"),
        clusters = path.join(out_dir, "workup", "condition-clusters", "{experiment}.{condition}.clusters")
    output:
        bam = path.join(out_dir, "workup", "splitbams_by_condition", "{experiment}.{condition}.bam"),
        touch = touch(path.join(out_dir, "workup", "splitbams_by_condition", "{experiment}.{condition}.done")),
        log = path.join(out_dir, "workup", "qc", "{experiment}.thresh_and_split_condition.{condition}.log")
    params:
        directory = "workup/splitbams_by_condition",
        max_size = config["max_size"],
        proportion  = config["proportion"],
        min_oligos = config["min_oligos"],
        num_tags = config["num_tags"]
    conda:
        "envs/sprite.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 32000,
        time = "00:15:00"
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.{condition}.splitbams.log")
    benchmark:
        out_dir + "benchmarks/{experiment}.{condition}.thresh_and_split_control.tsv"
    shell:
        '''
        (python scripts/python/threshold_tag_and_split.py \
            -i {input.bam} \
            -c {input.clusters} \
            -o {output.bam} \
            -d {params.directory} \
            -l {output.log} \
            --min_oligos {params.min_oligos} \
            --proportion {params.proportion} \
            --max_size {params.max_size} \
            --num_tags {params.num_tags}) &> {log}
        '''


rule thresh_and_split_no_condition:
    input:
        bam = path.join(out_dir, "workup", "merge_rna_bams", "{experiment}.merged.RPM.bam"),
        clusters = path.join(out_dir, "workup", "split_incorrect_clusters", "{experiment}.complete.clusters")
    output:
        bam = path.join(out_dir, "workup", "splitbams_all_conditions", "{experiment}.ALL_CONDITIONS.bam"),
        touch = touch(path.join(out_dir, "workup", "splitbams_all_conditions", "{experiment}.done")),
        log = path.join(out_dir, "workup", "qc", "{experiment}.thresh_and_split_no_condition.ALL_CONDITIONS.log")
    params:
        directory = "workup/splitbams_all_conditions",
        max_size = config["max_size"],
        proportion  = config["proportion"],
        min_oligos = config["min_oligos"],
        num_tags = config["num_tags"]
    conda:
        "envs/sprite.yaml"
    resources:
        tmpdir = config["temp_dir"],
        cpus = 1,
        mem_mb = 32000,
        time = "00:15:00"
    log:
        path.join("workup", "logs", "{experiment}.merged.splitbams.log")
    benchmark:
        path.join("benchmarks", "{experiment}.merged.thresh_and_split_control.tsv")
    shell:
        '''
        (python scripts/python/threshold_tag_and_split.py \
            -i {input.bam} \
            -c {input.clusters} \
            -o {output.bam} \
            -d {params.directory} \
            -l {output.log} \
            --min_oligos {params.min_oligos} \
            --proportion {params.proportion} \
            --max_size {params.max_size} \
            --num_tags {params.num_tags}) &> {log}
        '''


rule count_barcoded_reads_in_bams:
    """
    Count all reads across all bams that are not from ambiguous.bam, none.bam, uncertain.bam, 
    which contain reads that did not meet cluster cutoffs or the pooled bam which is the aggregate
    of all 
    """
    input:
        path.join(out_dir, "workup", "splitbams_all_conditions", "{experiment}.done"),
    output:
        path.join(out_dir, "workup", "count_barcoded_reads_in_bams", "{experiment}.barcoded_reads_assigned_to_bams.txt")
    params:
        directory = path.join(out_dir, "workup", "splitbams_all_conditions"),
    log:
        path.join(out_dir, "workup", "logs", "{experiment}.count_barcoded_reads_in_bams.log")
    conda:
        "envs/samtools.yaml"
    shell:
        """
        ALL_BAMS=$(ls workup/splitbams_all_conditions/*.bam | grep -E '.*\.ALL_CONDITIONS_.+\.bam' | grep -v -E 'ambiguous|none|uncertain')

        (for bam in $ALL_BAMS; do samtools view -c $bam; done | awk '{{sum += $1}} END {{print sum}}' > {output}) &> {log}
        """


rule generate_splitbam_statistics:
    input:
        expand([path.join(out_dir, "workup", "splitbams_all_conditions", "{experiment}.done")], experiment=ALL_EXPERIMENTS),
        expand([path.join(out_dir, "workup", "splitbams_by_condition", "{experiment}.done")], experiment=ALL_EXPERIMENTS)
    output:
        all_conditions = path.join(out_dir, "workup", "splitbams_all_conditions", "splitbam_statistics.txt"),
        by_condition = path.join(out_dir, "workup", "splitbams_by_condition", "splitbam_statistics.txt")
    params:
        all_conditions = path.join(out_dir, "workup", "splitbams_all_conditions"),
        by_condition = path.join(out_dir, "workup", "splitbams_by_condition")
    log:
        path.join(out_dir, "workup", "logs", "generate_splitbam_statistics.log")
    conda:
        "envs/sprite.yaml"
    resources:
        tmpdir = config["temp_dir"]
    shell:
        """
        (for f in {params.all_conditions}/*bam; do echo $f; samtools view -c $f; done > {output.all_conditions}) &> {log}
        (for f in {params.by_condition}/*bam; do echo $f; samtools view -c $f; done > {output.by_condition}) &> {log}
        """
