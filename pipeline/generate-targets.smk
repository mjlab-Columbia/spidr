import json
import os 
import sys
import numpy as np
import datetime
from pathlib import Path

configfile: "config.generate-targets.yaml"

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
    'sequence_length'
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
            print(f"Missing: {key} in {configfile}")
    
    sys.exit()

# For optional fields, augment the value or fall back to defaults
if 'cutadapt_oligos' not in config.keys():
    config['cutadapt_oligos'] = "-g GGTGGTCTTT -g GCCTCTTGTT"
else:
    config['cutadapt_oligos'] = "-g file:" + config['cutadapt_oligos']

# Create cluster subdirectory within logs/ directory manually
os.makedirs(
    os.path.join(out_dir, "workup", "logs", "cluster"),
    exist_ok=True
)

# Create directory for benchmark tsv files to be stored
os.makedirs(
    os.path.join("benchmarks"),
    exist_ok=True
)

################################################################################
#Get experiment files
###############################################################################

# Prep experiments from fastq directory using fastq2json_updated.py, now load json file
FILES = json.load(open(config['experiments']))
ALL_EXPERIMENTS = sorted(FILES.keys())

ALL_FASTQ = []
for SAMPLE, file in FILES.items():
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R1')])
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R2')])

NUM_CHUNKS = [f"{i:03}" for i in np.arange(0, config['num_chunks'])]

OUTPUTS = expand(
    [
        os.path.join(out_dir, "workup", "splitbams-all-conditions", "{experiment}.done"),
        os.path.join(out_dir, "workup", "splitbams-by-condition", "{experiment}.{condition}.done"),
        os.path.join(out_dir, "workup", "ligation_efficiency.txt"),
        os.path.join(out_dir, "workup", "{experiment}.part_{splitid}.ligation_efficiency.txt"),
        os.path.join(out_dir, "workup", "clusters", "cluster_statistics.txt"),
        os.path.join(out_dir, "workup", "condition-clusters", "RPM_read_distribution.pdf"),
        os.path.join(out_dir, "workup", "condition-clusters", "RPM_cluster_distribution.pdf"),
        os.path.join(out_dir, "workup", "condition-clusters", "BPM_read_distribution.pdf"),
        os.path.join(out_dir, "workup", "condition-clusters", "BPM_cluster_distribution.pdf"),
        os.path.join(out_dir, "workup", "clusters", "RPM_read_distribution.pdf"),       
        os.path.join(out_dir, "workup", "clusters", "RPM_cluster_distribution.pdf"),       
        os.path.join(out_dir, "workup", "clusters", "BPM_read_distribution.pdf"),
        os.path.join(out_dir, "workup", "clusters", "BPM_cluster_distribution.pdf"),
        os.path.join(out_dir, "workup", "clusters", "Max_representation_ecdf.pdf"),
        os.path.join(out_dir, "workup", "clusters", "Max_representation_counts.pdf")
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
            os.path.join(out_dir, "workup", "splitfq", "{{experiment}}_R1.part_{splitid}.fastq"),
            splitid=NUM_CHUNKS
        )
    params:
        dir = out_dir + "workup/splitfq",
        prefix_r1 = "{experiment}_R1.part_0",
        prefix_r2 = "{experiment}_R2.part_0",
        num_chunks = config['num_chunks']
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.split_fastq_read1.log")
    conda:
        "envs/sprite.yaml"
    threads: 
        8
    benchmark:
        "benchmarks/{experiment}.split_fastq_read1.tsv"
    shell:
        '''
        mkdir -p {params.dir}
        
        (bash scripts/bash/split_fastq.sh \
            {input.r1} \
            {params.num_chunks} \
            {params.dir} \
            {params.prefix_r1}) &> {log}
        '''


rule split_fastq_read2:
    input:
        r2 = lambda wildcards: FILES[wildcards.experiment]['R2']
    output:
        expand(
            os.path.join(out_dir, "workup", "splitfq", "{{experiment}}_R2.part_{splitid}.fastq"),
            splitid=NUM_CHUNKS
        )
    params:
        dir = out_dir + "workup/splitfq",
        prefix_r1 = "{experiment}_R1.part_0",
        prefix_r2 = "{experiment}_R2.part_0",
        num_chunks = config['num_chunks']
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.split_fastq_read2.log")
    conda:
        "envs/sprite.yaml"
    threads: 
        8
    benchmark:
        "benchmarks/{experiment}.split_fastq_read2.tsv"
    shell:
        '''
        mkdir -p {params.dir}

        (bash scripts/bash/split_fastq.sh \
            {input.r2} \
            {params.num_chunks} \
            {params.dir} \
            {params.prefix_r2}) &> {log}
        '''


rule compress_fastq_read1:
    input:
        r1 = os.path.join(out_dir, "workup", "splitfq", "{experiment}_R1.part_{splitid}.fastq"),
    output:
        r1 = os.path.join(out_dir, "workup", "splitfq", "{experiment}_R1.part_{splitid}.fastq.gz"),
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.compress_fastq_read1.log")
    conda:
        "envs/sprite.yaml"
    threads:
        8
    benchmark:
        "benchmarks/{experiment}.{splitid}.compress_fastq_read1.tsv"
    shell:
        '''
        (pigz -p {threads} {input.r1}) &> {log}
        '''


rule compress_fastq_read2:
    input:
        r2 = os.path.join(out_dir, "workup", "splitfq", "{experiment}_R2.part_{splitid}.fastq")
    output:
        r2 = os.path.join(out_dir, "workup", "splitfq", "{experiment}_R2.part_{splitid}.fastq.gz")
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.compress_fastq_read2.log")
    conda:
        "envs/sprite.yaml"
    threads:
        8
    benchmark:
        "benchmarks/{experiment}.{splitid}.compress_fastq_read2.tsv"
    shell:
        '''
        (pigz -p {threads} {input.r2}) &> {log}
        '''


rule trim_sequencing_adapters:
    input:
        [os.path.join(out_dir, "workup", "splitfq", "{experiment}_R1.part_{splitid}.fastq.gz"), 
         os.path.join(out_dir, "workup", "splitfq", "{experiment}_R2.part_{splitid}.fastq.gz")]
    output:
        os.path.join(out_dir, "workup", "trimmed", "{experiment}_R1.part_{splitid}_val_1.fq.gz"),
        os.path.join(out_dir, "workup", "trimmed", "{experiment}_R1.part_{splitid}.fastq.gz_trimming_report.txt"),
        os.path.join(out_dir, "workup", "trimmed", "{experiment}_R2.part_{splitid}_val_2.fq.gz"),
        os.path.join(out_dir, "workup", "trimmed", "{experiment}_R2.part_{splitid}.fastq.gz_trimming_report.txt")
    threads:
        10
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.trim_sequencing_adapters.log")
    conda:
        "envs/sprite.yaml"
    benchmark:
        "benchmarks/{experiment}.{splitid}.trim_sequencing_adapters.tsv"
    shell:
        '''
        (trim_galore \
            --paired \
            --gzip \
            --cores {threads} \
            --quality 20 \
            --fastqc \
            -o {out_dir}workup/trimmed/ \
            {input}) &> {log}
        '''


rule identify_barcodes:
    input:
        r1 = os.path.join(out_dir, "workup", "trimmed", "{experiment}_R1.part_{splitid}_val_1.fq.gz"),
        r2 = os.path.join(out_dir, "workup", "trimmed", "{experiment}_R2.part_{splitid}_val_2.fq.gz")
    output:
        r1_barcoded = os.path.join(out_dir, "workup", "fastqs", "{experiment}_R1.part_{splitid}.barcoded.fastq.gz"),
        r2_barcoded = os.path.join(out_dir, "workup", "fastqs", "{experiment}_R2.part_{splitid}.barcoded.fastq.gz")
    params:
        bid_config = config['bID'],
        read1_format = config["read1_format"],
        read2_format = config["read2_format"],
        start_offset = config["read1_start_offset"],
        sequence_length = config["sequence_length"]
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.identify_barcodes.log")
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
            --start_offset {params.start_offset} \
            --sequence_length {params.sequence_length} \
            --config {params.bid_config}) &> {log}
        """



rule get_ligation_efficiency:
    input:
        r1 = os.path.join(out_dir, "workup", "fastqs", "{experiment}_R1.part_{splitid}.barcoded.fastq.gz")
    output:
        os.path.join(out_dir, "workup", "{experiment}.part_{splitid}.ligation_efficiency.txt")
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.get_ligation_efficiency.log")
    conda:
        "envs/sprite.yaml"
    shell:
        """
        (python scripts/python/get_ligation_efficiency.py {input.r1} > {output}) &> {log}
        """


rule cat_ligation_efficiency:
    input:
        expand(
            os.path.join(out_dir, "workup", "{experiment}.part_{splitid}.ligation_efficiency.txt"), 
            experiment=ALL_EXPERIMENTS, 
            splitid=NUM_CHUNKS
        )
    output:
        out_dir + "workup/ligation_efficiency.txt"
    log:
        os.path.join(out_dir, "workup", "logs", "cat_ligation_efficiency.log")
    shell:
        """
        (tail -n +1 {input} > {output}) &> {log}
        """


rule split_reads_read1:
    '''
    split bpm and rpm will also remove incomplete barcodes
    '''
    input:
        os.path.join(out_dir, "workup", "fastqs", "{experiment}_R1.part_{splitid}.barcoded.fastq.gz")
    output:
        os.path.join(out_dir, "workup", "fastqs", "{experiment}_R1.part_{splitid}.barcoded_rpm.fastq.gz"),
        os.path.join(out_dir, "workup", "fastqs", "{experiment}_R1.part_{splitid}.barcoded_bpm.fastq.gz"),
        os.path.join(out_dir, "workup", "fastqs", "{experiment}_R1.part_{splitid}.barcoded_short.fastq.gz")
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.split_reads_read1.log")
    conda:
       "envs/sprite.yaml"
    benchmark:
        "benchmarks/{experiment}.{splitid}.split_reads_read1.tsv"
    shell:
        """
        (python scripts/python/split_rpm_bpm_fq.py --r1 {input}) &> {log}
        """


rule split_reads_read2:
    '''
    split bpm and rpm will also remove incomplete barcodes
    '''
    input:
        os.path.join(out_dir, "workup", "fastqs", "{experiment}_R2.part_{splitid}.barcoded.fastq.gz")
    output:
        os.path.join(out_dir, "workup", "fastqs", "{experiment}_R2.part_{splitid}.barcoded_rpm.fastq.gz"),
        os.path.join(out_dir, "workup", "fastqs", "{experiment}_R2.part_{splitid}.barcoded_bpm.fastq.gz"),
        os.path.join(out_dir, "workup", "fastqs", "{experiment}_R2.part_{splitid}.barcoded_short.fastq.gz")
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.split_reads_read2.log")
    conda:
       "envs/sprite.yaml"
    benchmark:
        "benchmarks/{experiment}.{splitid}.split_reads_read2.tsv"
    shell:
        """
        (python scripts/python/split_rpm_bpm_fq.py --r1 {input}) &> {log}
        """


rule trim_rpm_reads:
    input:
        read1 = os.path.join(out_dir, "workup", "fastqs", "{experiment}_R1.part_{splitid}.barcoded_rpm.fastq.gz"),
        read2 = os.path.join(out_dir, "workup", "fastqs", "{experiment}_R2.part_{splitid}.barcoded_rpm.fastq.gz")
    output:
        r1 = os.path.join(out_dir, "workup", "trimmed", "{experiment}_R1.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz"),
        r2 = os.path.join(out_dir, "workup", "trimmed", "{experiment}_R2.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz"),
        qc = os.path.join(out_dir, "workup", "trimmed", "{experiment}_R1.part_{splitid}.barcoded_rpm.RDtrim.qc.txt")
    params:
        adapters_r1 = "-a ATCAGCACTTAGCGTCAG",
        adapters_r2 = "-G CTGACGCTAAGTGCTGAT",
        others = "--minimum-length 20"
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.RPM.cutadapt.log")
    threads: 
        10
    conda:
        "envs/sprite.yaml"
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
            -j {threads} \
            {input.read1} {input.read2} > {output.qc}) &> {log}
        
        (fastqc {output.r1}) &>> {log}
        fastqc {output.r2} &>> {log}
        '''


rule trim_bead_oligo_reads:
    '''
    Trim 9mer oligo sequence from bead barcode
    '''
    input:
        os.path.join(out_dir, "workup", "fastqs", "{experiment}_R1.part_{splitid}.barcoded_bpm.fastq.gz")
    output:
        fastq = os.path.join(out_dir, "workup", "trimmed", "{experiment}_R1.part_{splitid}.barcoded_bpm.RDtrim.fastq.gz"),
        qc = os.path.join(out_dir, "workup", "trimmed", "{experiment}_R1.part_{splitid}.barcoded_bpm.RDtrim.qc.txt")
    params:
        adapters_r1 = config['cutadapt_oligos']
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.BPM.cutadapt.log")
    threads: 
        10
    conda:
        "envs/sprite.yaml"
    benchmark:
        "benchmarks/{experiment}.{splitid}.trim_bead_oligo_reads.tsv"
    shell:
        '''
        (cutadapt \
            {params.adapters_r1} \
            -o {output.fastq} \
            -j {threads} \
            {input} > {output.qc}) &> {log}
        '''


rule align_bowtie2:
    input:
        fq1 = os.path.join(out_dir, "workup", "trimmed", "{experiment}_R1.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz"),
        fq2 = os.path.join(out_dir, "workup", "trimmed", "{experiment}_R2.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz")
    output:
        bam = temp(os.path.join(out_dir, "workup", "alignments", "{experiment}.part_{splitid}.bowtie2.bam")),
        mapped = os.path.join(out_dir, "workup", "alignments", "{experiment}.part_{splitid}.bowtie2.sorted.mapped.bam"),
        unmapped = os.path.join(out_dir, "workup", "alignments", "{experiment}.part_{splitid}.bowtie2.sorted.unmapped.bam"),
        index = os.path.join(out_dir, "workup", "alignments", "{experiment}.part_{splitid}.bowtie2.sorted.mapped.bam.bai")
    params:
        BOWTIE2_INDEX = config['bowtie2_index'][config['assembly']]
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.bt2.log")
    threads: 
        4
    conda:
        "envs/sprite.yaml"
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


rule convert_bam_to_fastq:
    input: 
        os.path.join(out_dir, "workup", "alignments", "{experiment}.part_{splitid}.bowtie2.sorted.unmapped.bam")
    output: 
        r1 = os.path.join(out_dir, "workup", "fastqs", "{experiment}.part_{splitid}.bowtie2.unmapped_R1.fq.gz"),
        r2 = os.path.join(out_dir, "workup", "fastqs", "{experiment}.part_{splitid}.bowtie2.unmapped_R2.fq.gz")
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.bam2fq.log")
    threads: 
        1
    conda:
        "envs/sprite.yaml"
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
        r1 = os.path.join(out_dir, "workup", "fastqs", "{experiment}.part_{splitid}.bowtie2.unmapped_R1.fq.gz"),
        r2 = os.path.join(out_dir, "workup", "fastqs", "{experiment}.part_{splitid}.bowtie2.unmapped_R2.fq.gz")
    output:
        sam = temp(os.path.join(out_dir, "workup", "alignments", "{experiment}.part_{splitid}.Aligned.out.sam")),
        sorted = os.path.join(out_dir, "workup", "alignments", "{experiment}.part_{splitid}.Aligned.out.sorted.bam"),
        filtered = temp(os.path.join(out_dir, "workup", "alignments", "{experiment}.part_{splitid}.Aligned.out.bam"))
    params:
        STAR_OPTIONS = "--readFilesCommand zcat --alignEndsType EndToEnd --outFilterScoreMin 10 --outFilterMultimapNmax 1 --outFilterMismatchNmax 10 --alignIntronMax 100000 --alignMatesGapMax 1300 --alignIntronMin 80 --alignSJDBoverhangMin 5 --alignSJoverhangMin 8 --chimSegmentMin 20 --alignSJstitchMismatchNmax 5 -1 5 5 --outSAMunmapped Within --outReadsUnmapped Fastx", prefix = out_dir + "workup/alignments/{experiment}.part_{splitid}.",
        STAR_INDEX = config['star_index'][config['assembly']]
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.star.log")
    threads:
        10
    conda:
        "envs/sprite.yaml"
    benchmark:
        "benchmarks/{experiment}.{splitid}.align_star.tsv"
    shell:
        '''
        (STAR \
            --genomeDir {params.STAR_INDEX} \
            --readFilesIn {input.r1} {input.r2} \
            --runThreadN {threads} {params.STAR_OPTIONS} \
            --outFileNamePrefix {params.prefix}) &> {log}

        (samtools view -@ {threads} -bS {output.sam} > {output.filtered}) &>> {log}
        (samtools sort -@ {threads} -o {output.sorted} {output.filtered}) &>> {log}
        (samtools index {output.sorted}) &>> {log}
        '''


rule add_chromosome_info_bowtie2:
    input:
        bt2 = os.path.join(out_dir, "workup", "alignments", "{experiment}.part_{splitid}.bowtie2.sorted.mapped.bam")
    output:
        bt2 = os.path.join(out_dir, "workup", "alignments", "{experiment}.part_{splitid}.bowtie2.sorted.mapped.chr.bam")
    params:
        assembly = config['assembly']
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.add_chr.log"),
    conda:
        "envs/sprite.yaml"
    benchmark:
        "benchmarks/{experiment}.{splitid}.add_chr.tsv"
    shell:
        '''
        (python scripts/python/add_chr_bt2.py \
            -i {input.bt2} \
            -o {output.bt2} \
            --assembly {params.assembly}) &> {log}
        '''


rule add_chromosome_info_star:
    input:
        star = os.path.join(out_dir, "workup", "alignments", "{experiment}.part_{splitid}.Aligned.out.sorted.bam"),
    output:
        star = os.path.join(out_dir, "workup", "alignments", "{experiment}.part_{splitid}.Aligned.out.sorted.chr.bam"),
    params:
        assembly = config['assembly']
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.add_chr.log"),
    conda:
        "envs/sprite.yaml"
    benchmark:
        "benchmarks/{experiment}.{splitid}.add_chr.tsv"
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
            os.path.join(out_dir + "workup/alignments/{{experiment}}.part_{splitid}.bowtie2.sorted.mapped.chr.bam"),
            splitid=NUM_CHUNKS
        ),
        star = expand(
            os.path.join(out_dir + "workup/alignments/{{experiment}}.part_{splitid}.Aligned.out.sorted.chr.bam"), 
            splitid=NUM_CHUNKS
        )
    output:
        os.path.join(out_dir, "workup", "alignments", "{experiment}.merged.RPM.bam")
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.merge_rna_bams.log")
    conda:
        "envs/sprite.yaml"
    threads:
        8
    benchmark:
        "benchmarks/{experiment}.merge_rna_bams.tsv"
    shell:
        '''
        (samtools merge -@ {threads} {output} {input.bt2} {input.star}) &> {log}
        (samtools index {output}) &>> {log}
        '''


rule convert_fastq_to_bam:
    input:
        os.path.join(out_dir, "workup", "trimmed", "{experiment}_R1.part_{splitid}.barcoded_bpm.RDtrim.fastq.gz")
    output:
        sorted = os.path.join(out_dir, "workup", "alignments", "{experiment}.part_{splitid}.BPM.bam"),
        bam = temp(os.path.join(out_dir, "workup", "alignments", "{experiment}.part_{splitid}.BPM.unsorted.bam"))
    params:
        bid_config = config['bID']
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.make_bam.log")
    conda:
        "envs/sprite.yaml"
    threads:
        8
    benchmark:
        "benchmarks/{experiment}.{splitid}.convert_fastq_to_bam.tsv"
    shell:
        '''
        python scripts/python/fastq_to_bam.py --input {input} --output {output.bam} --config {params.bid_config} &> {log}
        (samtools sort -@ {threads} -o {output.sorted} {output.bam}) &>> {log}
        '''


rule merge_bead_bams:
    input:
        expand(
            os.path.join(out_dir, "workup", "alignments", "{{experiment}}.part_{splitid}.BPM.bam"), 
            splitid = NUM_CHUNKS
        )
    output:
        os.path.join(out_dir, "workup", "alignments", "{experiment}.merged.BPM.bam")
    conda:
        "envs/sprite.yaml"
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.merge_bead_bams.log")
    threads:
        8
    benchmark:
        "benchmarks/{experiment}.merge_bead_bams.tsv"
    shell:
        '''
        (samtools merge -@ {threads} {output} {input}) &> {log}
        '''


rule make_clusters:
    input:
        rpm = os.path.join(out_dir, "workup", "alignments", "{experiment}.part_{splitid}.Aligned.out.sorted.chr.bam"),
        bpm = os.path.join(out_dir, "workup", "alignments", "{experiment}.part_{splitid}.BPM.bam"),
        bt2 = os.path.join(out_dir, "workup", "alignments", "{experiment}.part_{splitid}.bowtie2.sorted.mapped.chr.bam")
    output:
        unsorted = temp(os.path.join(out_dir, "workup", "clusters", "{experiment}.part_{splitid}.unsorted.clusters")),
        sorted = os.path.join(out_dir, "workup", "clusters", "{experiment}.part_{splitid}.clusters")
    params:
        temp_dir = config['temp_dir'],
        num_tags = config['num_tags']
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.{splitid}.make_clusters.log")
    conda:
        "envs/sprite.yaml"
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
            os.path.join(out_dir, "workup", "clusters", "{{experiment}}.part_{splitid}.clusters"), 
            splitid=NUM_CHUNKS
        )
    output:
        mega = temp(os.path.join(out_dir, "workup", "clusters", "{experiment}.duplicated.clusters")),
        final = os.path.join(out_dir, "workup", "clusters", "{experiment}.clusters")
    params:
        temp_dir = config['temp_dir']
    conda:
       "envs/sprite.yaml"
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.merge_clusters.log")
    benchmark:
        "benchmarks/{experiment}.merge_clusters.tsv"
    shell:
        '''
        (sort -k 1 -T {params.temp_dir} -m {input} > {output.mega}) &> {log}
        (python scripts/python/merge_clusters.py -i {output.mega} -o {output.final}) &>> {log}
        '''        


rule generate_cluster_statistics:
    input:
        expand(
            os.path.join(out_dir, "workup", "clusters", "{experiment}.complete.clusters"), 
            experiment = ALL_EXPERIMENTS
        )
    output:
        os.path.join(out_dir, "workup", "clusters", "cluster_statistics.txt")
    params:
        dir = out_dir + "workup/clusters"
    conda:
        "envs/sprite.yaml"
    log:
        os.path.join(out_dir, "workup", "logs", "generate_cluster_statistics.log")
    shell:
        """
        (python scripts/python/generate_cluster_statistics.py \
            --directory {params.dir} \
            --pattern .clusters > {output}) &> {log}
        """


rule generate_cluster_ecdfs:
    input:
        expand(
            os.path.join(out_dir, "workup", "clusters", "{experiment}.clusters"), 
            experiment=ALL_EXPERIMENTS
        )
    output:
        ecdf = os.path.join(out_dir, "workup", "clusters", "Max_representation_ecdf.pdf"),
        counts = os.path.join(out_dir, "workup", "clusters", "Max_representation_counts.pdf")
    params:
        dir = out_dir + "workup/clusters"
    conda:
        "envs/plotting.yaml"
    log:
        os.path.join(out_dir, "workup", "logs", "generate_cluster_ecdfs.log")        
    shell:
        """
        (python scripts/python/max_representation_ecdfs_perlib.py \
            --directory {params.dir} \
            --pattern .clusters \
            --xlim 30) &> {log}
        """


# Profile size distribution of clusters
rule get_size_distribution:
    input:
        expand(
            os.path.join(out_dir, "workup", "condition-clusters", "{experiment}.{condition}.clusters"), 
            experiment = ALL_EXPERIMENTS, 
            condition = config['conditions']
        ),
        expand(
            os.path.join(out_dir, "workup", "clusters", "{experiment}.complete.clusters"), 
            experiment=ALL_EXPERIMENTS
        ),
    output:
        rpm = os.path.join(out_dir, "workup", "condition-clusters", "RPM_read_distribution.pdf"),
        rpm2 = os.path.join(out_dir, "workup", "condition-clusters", "RPM_cluster_distribution.pdf"),
        bpm = os.path.join(out_dir, "workup", "condition-clusters", "BPM_read_distribution.pdf"),
        bpm2 = os.path.join(out_dir, "workup", "condition-clusters", "BPM_cluster_distribution.pdf"),
        no_condition_rpm = os.path.join(out_dir, "workup", "clusters", "RPM_read_distribution.pdf"),       
        no_condition_rpm2 = os.path.join(out_dir, "workup", "clusters", "RPM_cluster_distribution.pdf"),       
        no_condition_bpm = os.path.join(out_dir, "workup", "clusters", "BPM_read_distribution.pdf"),
        no_condition_bpm2 = os.path.join(out_dir, "workup", "clusters", "BPM_cluster_distribution.pdf")   
    params:
        condition_dir = "workup/condition-clusters",
        no_condition_dir = "workup/clusters"
    log:
        os.path.join(out_dir, "workup", "logs", "get_size_distribution.log")
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        (python scripts/python/get_bead_size_distribution.py --directory {params.no_condition_dir} --pattern .clusters --readtype BPM) &> {log}
        (python scripts/python/get_bead_size_distribution.py --directory {params.no_condition_dir} --pattern .clusters --readtype RPM) &>> {log}

        (python scripts/python/get_bead_size_distribution.py --directory {params.condition_dir} --pattern .clusters --readtype BPM) &>> {log}
        (python scripts/python/get_bead_size_distribution.py --directory {params.condition_dir} --pattern .clusters --readtype RPM) &>> {log}
        '''


rule split_incorrect_clusters:
    input:
        clusters = os.path.join(out_dir, "workup", "clusters", "{experiment}.clusters")
    output:
        complete_clusters = os.path.join(out_dir, "workup", "clusters", "{experiment}.complete.clusters"),
        incomplete_clusters = os.path.join(out_dir, "workup", "clusters", "{experiment}.incomplete.clusters")
    params:
        rounds_format = config["rounds_format"]
    conda:
        "envs/sprite.yaml"
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.incorrectclusters.log")
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


rule split_on_first_tag:
    input:
        complete_clusters = os.path.join(out_dir, "workup", "clusters", "{experiment}.complete.clusters")
    output:
        expand(
            os.path.join(out_dir, "workup", "condition-clusters", "{{experiment}}.{condition}.clusters"),
            condition = config['conditions']
        )
    conda:
        "envs/sprite.yaml"
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.splitonfirsttag.log")
    benchmark:
        "benchmarks/{experiment}.split_on_first_tag.tsv"
    shell:
        '''
        (python scripts/python/split_on_first_tag.py \
            --complete_clusters {input.complete_clusters} \
            --output_dir workup/condition-clusters) &> {log}
        '''


# Generate bam files for individual targets based on assignments from clusterfile
rule thresh_and_split_condition:
    input:
        bam = os.path.join(out_dir, "workup", "alignments", "{experiment}.merged.RPM.bam"),
        clusters = os.path.join(out_dir, "workup", "condition-clusters", "{experiment}.{condition}.clusters")
    output:
        bam = os.path.join(out_dir, "workup", "splitbams-by-condition", "{experiment}.{condition}.bam"),
        touch = touch(os.path.join(out_dir, "workup", "splitbams-by-condition", "{experiment}.{condition}.done"))
    params:
        directory = "workup/splitbams-by-condition",
        max_size = config["max_size"],
        proportion  = config["proportion"],
        min_oligos = config["min_oligos"],
        num_tags = config["num_tags"]
    conda:
        "envs/sprite.yaml"
    log:
        os.path.join(out_dir, "workup", "logs", "{experiment}.{condition}.splitbams.log")
    benchmark:
        out_dir + "benchmarks/{experiment}.{condition}.thresh_and_split_control.tsv"
    shell:
        '''
        mkdir -p splitbams_tmpdir
        export TMPDIR=splitbams_tmpdir/

        (python scripts/python/threshold_tag_and_split.py \
            -i {input.bam} \
            -c {input.clusters} \
            -o {output.bam} \
            -d {params.directory} \
            --min_oligos {params.min_oligos} \
            --proportion {params.proportion} \
            --max_size {params.max_size} \
            --num_tags {params.num_tags}) &> {log}
        '''


rule thresh_and_split_no_condition:
    input:
        bam = os.path.join(out_dir, "workup", "alignments", "{experiment}.merged.RPM.bam"),
        clusters = os.path.join(out_dir, "workup", "clusters", "{experiment}.complete.clusters")
    output:
        bam = os.path.join(out_dir, "workup", "splitbams-all-conditions", "{experiment}.ALL_CONDITIONS.bam"),
        touch = touch(os.path.join(out_dir, "workup", "splitbams-all-conditions", "{experiment}.done"))
    params:
        directory = "workup/splitbams-all-conditions",
        max_size = config["max_size"],
        proportion  = config["proportion"],
        min_oligos = config["min_oligos"],
        num_tags = config["num_tags"]
    conda:
        "envs/sprite.yaml"
    log:
        out_dir + "workup/logs/{experiment}.merged.splitbams.log"
    benchmark:
        out_dir + "benchmarks/{experiment}.merged.thresh_and_split_control.tsv"
    shell:
        '''
        mkdir -p splitbams_tmpdir
        export TMPDIR=splitbams_tmpdir/

        (python scripts/python/threshold_tag_and_split.py \
            -i {input.bam} \
            -c {input.clusters} \
            -o {output.bam} \
            -d {params.directory} \
            --min_oligos {params.min_oligos} \
            --proportion {params.proportion} \
            --max_size {params.max_size} \
            --num_tags {params.num_tags}) &> {log}
        '''


rule generate_splitbam_statistics:
    input:
        expand([os.path.join(out_dir, "workup", "splitbams-all-conditions", "{experiment}.done")], experiment=ALL_EXPERIMENTS),
        expand([os.path.join(out_dir, "workup", "splitbams-by-condition", "{experiment}.done")], experiment=ALL_EXPERIMENTS)
    output:
        all_conditions = os.path.join(out_dir, "workup", "splitbams-all-conditions", "splitbam_statistics.txt"),
        by_condition = os.path.join(out_dir, "workup", "splitbams-by-condition", "splitbam_statistics.txt")
    params:
        all_conditions = os.path.join(out_dir, "workup", "splitbams-all-conditions"),
        by_condition = os.path.join(out_dir, "workup", "splitbams-by-condition")
    log:
        os.path.join(out_dir, "workup", "logs", "generate_splitbam_statistics.log")
    conda:
        "envs/sprite.yaml"
    shell:
        """
        (for f in {params.all_conditions}/*bam; do echo $f; samtools view -c $f; done > {output.all_conditions}) &> {log}
        (for f in {params.by_condition}/*bam; do echo $f; samtools view -c $f; done > {output.by_condition}) &> {log}
        """
