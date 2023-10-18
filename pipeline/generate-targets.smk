'''
Aim: A Snakemake workflow to process SPIDR paired end data (in-progress)
'''

import os 
import sys
import numpy as np
import datetime
from pathlib import Path

##############################################################################
#Initialize settings
##############################################################################

#Copy config file into logs
v = datetime.datetime.now()
run_date = v.strftime('%Y.%m.%d.')

try:
    config_path = config["config_path"]
except:
    config_path = 'config.generate-targets.yaml'

configfile: config_path

try:
    email = config['email']
except:
    email = None
    print("Will not send email on error")

##############################################################################
#Location of scripts
##############################################################################

barcode_id_jar = "scripts/java/BarcodeIdentification_v1.2.0.jar"
lig_eff = "scripts/python/get_ligation_efficiency.py"
split_bpm_rpm = "scripts/python/split_rpm_bpm_fq.py"
add_chr = "scripts/python/ensembl2ucsc.py"
get_clusters = "scripts/python/get_clusters.py"
merge_clusters = "scripts/python/merge_clusters.py"
label_clusters = "scripts/python/label_clusters.py"
fq_to_bam = "scripts/python/fastq_to_bam.py"
add_RG_to_bam = "scripts/python/add_tag_to_bam.py"
split_fastq = "scripts/split_fastq.sh"
add_chr_bt2 = "scripts/python/add_chr_bt2.py"

cluster_counts = "scripts/python/generate_cluster_statistics.py"
cluster_sizes = "scripts/python/get_bead_size_distribution.py"
cluster_ecdfs = "scripts/python/max_representation_ecdfs_perlib.py"
split_incorrect_clusters = "scripts/split_incorrect_clusters.py"
tag_and_split = "scripts/threshold_tag_and_split.py"
random_downsample_clip = "scripts/java/RandomDownsample_CLIP.jar"
split_on_first_tag = "scripts/python/split_on_first_tag.py"

##############################################################################
#General settings
##############################################################################

try:
    bid_config = config['bID']
    print('Using BarcodeID config', bid_config)
except:
    bid_config = 'config.txt'
    print('Config "bID" not specified, looking for config at:', bid_config)

try:
    num_tags = config['num_tags']
    print('Using', num_tags, 'tags')
except:
    num_tags = "6"
    print('Config "num_tags" not specified, using:', num_tags)

try:
    conditions = config['conditions']
except:
    print("No conditions specified in config file. No defaults.")
    sys.exit()

try:
    assembly = config['assembly']
    assert assembly in ['mm10', 'hg38'], 'Only "mm10" or "hg38" currently supported'
    print('Using', assembly)
except:
    print('Config "assembly" not specified, defaulting to "mm10"')
    assembly = 'mm10'

try:
    experiments = config['experiments']
    print('Using experiments file:', experiments)
except:
    experiments = './experiments.json'
    print('Defaulting to working directory for experiments json file')

try:
    out_dir = config['output_dir']
    print('All data will be written to:', out_dir)
except:
    out_dir = ''
    print('Defaulting to working directory as output directory')

try:
    temp_dir = config['temp_dir']
    print("Using temporary directory:", temp_dir)
except:
    temp_dir = '/central/scratch/'
    print('Defaulting to central scratch as temporary directory')

try:
    num_chunks = config['num_chunks']
except:
    num_chunks = 2

try:
    rounds_format = config['rounds_format']
except:
    rounds_format = "format_6_rounds.txt"

##############################################################################
# Load Post Clustering Setting
##############################################################################

try:
    generate_splitbams = config["generate_splitbams"]
except:
    generate_splitbams = False

try:
    min_oligos = config["min_oligos"]
except:
    min_oligos = 2

try:
    proportion = config["proportion"]
except:
    proportion = 0.8

try:
    max_size = config["max_size"]
except:
    max_size = 10000

if generate_splitbams:
    print("Will generate bam files for individual targets using:", file=sys.stderr)
    print("\t min_oligos: ", min_oligos, file=sys.stderr)
    print("\t proportion: ", proportion, file=sys.stderr)
    print("\t max_size: ", max_size, file=sys.stderr)
else:
    print("Will not generate bam files for individual targets.", file=sys.stderr)

##############################################################################
##Trimming Sequences
##############################################################################

try:
    adapters = "-g file:" + config['cutadapt_dpm']
    print('Using cutadapt sequence file', adapters)
except:
    adapters = "-g GGTGGTCTTT -g GCCTCTTGTT \
        -g CCAGGTATTT -g TAAGAGAGTT -g TTCTCCTCTT -g ACCCTCGATT"
    print("No file provided for cutadapt. Using standard cutadapt sequences")

try:
    oligos = "-g file:" + config['cutadapt_oligos']
    print('Using bead oligo file', oligos)
except:
    oligos = "-g GGTGGTCTTT -g GCCTCTTGTT"
    print("Using junk oligos. FIX ME")

################################################################################
#Aligner Indexes
################################################################################

#RNA aligner

try:
    bowtie2_index = config['bowtie2_index'][config['assembly']]
except:
    print('Bowtie2 index not specified in config.yaml')
    sys.exit() #no default, exit

try:
    star_index = config['star_index'][config['assembly']]
except:
    print('Star index not specified in config.yaml')
    sys.exit() #no default, exit

################################################################################
#make output directories (aren't created automatically on cluster)
################################################################################

# Path(out_dir + "workup/logs/cluster").mkdir(parents=True, exist_ok=True)
os.makedirs(out_dir + "workup/logs/cluster", exist_ok=True)
out_created = os.path.exists(out_dir + "workup/logs/cluster")
print('Output logs path created:', out_created)

# Create directory for benchmark tsv files to be stored
os.makedirs("benchmarks", exist_ok=True)

################################################################################
#Get experiment files
###############################################################################

#Prep experiments from fastq directory using fastq2json_updated.py, now load json file
FILES = json.load(open(experiments))
ALL_EXPERIMENTS = sorted(FILES.keys())

ALL_FASTQ = []
for SAMPLE, file in FILES.items():
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R1')])
    ALL_FASTQ.extend([os.path.abspath(i) for i in file.get('R2')])

CONFIG = [out_dir + "workup/logs/config_" + run_date + "yaml"]

NUM_CHUNKS = [f"{i:03}" for i in np.arange(0, num_chunks)]

################################################################################
#Trimming
################################################################################

SPLIT_FQ = expand(out_dir + "workup/splitfq/{experiment}_{read}.part_{splitid}.fastq.gz", experiment=ALL_EXPERIMENTS, read = ["R1", "R2"], splitid=NUM_CHUNKS)

TRIM = expand([out_dir + "workup/trimmed/{experiment}_R1.part_{splitid}_val_1.fq.gz", out_dir + "workup/trimmed/{experiment}_R2.part_{splitid}_val_2.fq.gz"], experiment=ALL_EXPERIMENTS, splitid=NUM_CHUNKS)
TRIM_LOG = expand(out_dir + "workup/trimmed/{experiment}_{read}.part_{splitid}.fastq.gz_trimming_report.txt", experiment=ALL_EXPERIMENTS, read = ["R1", "R2"], splitid = NUM_CHUNKS)
TRIM_RD = expand([out_dir + "workup/trimmed/{experiment}_R1.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz",
                  out_dir + "workup/trimmed/{experiment}_R1.part_{splitid}.barcoded_bpm.RDtrim.fastq.gz",
                  out_dir + "workup/trimmed/{experiment}_R2.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz"],
                  experiment=ALL_EXPERIMENTS, splitid=NUM_CHUNKS)

################################################################################
#Logging
################################################################################

LE_LOG_ALL = [out_dir + "workup/ligation_efficiency.txt"]

################################################################################
#Barcoding
################################################################################

BARCODEID = expand(out_dir + "workup/fastqs/{experiment}_{read}.part_{splitid}.barcoded.fastq.gz", experiment=ALL_EXPERIMENTS, read = ["R1", "R2"], splitid=NUM_CHUNKS)

SPLIT_RPM_BPM = expand([out_dir + "workup/fastqs/{experiment}_R1.part_{splitid}.barcoded_bpm.fastq.gz",
                    out_dir + "workup/fastqs/{experiment}_R1.part_{splitid}.barcoded_rpm.fastq.gz"], experiment=ALL_EXPERIMENTS, splitid=NUM_CHUNKS)

SPLIT_RPM_BPM2 = expand([out_dir + "workup/fastqs/{experiment}_R2.part_{splitid}.barcoded_bpm.fastq.gz",
                    out_dir + "workup/fastqs/{experiment}_R2.part_{splitid}.barcoded_rpm.fastq.gz"], experiment=ALL_EXPERIMENTS, splitid=NUM_CHUNKS)

################################################################################
#RNA workup
################################################################################

BT2_RNA_ALIGN = expand([out_dir + "workup/alignments/{experiment}.part_{splitid}.bowtie2.sorted.mapped.bam",
                        out_dir + "workup/alignments/{experiment}.part_{splitid}.bowtie2.sorted.unmapped.bam",
                        out_dir + "workup/alignments/{experiment}.part_{splitid}.bowtie2.sorted.mapped.bam.bai"],
                        experiment=ALL_EXPERIMENTS, splitid=NUM_CHUNKS)

BAM_TO_FQ = expand([out_dir + "workup/fastqs/{experiment}.part_{splitid}.bowtie2.unmapped_R1.fq.gz",
                    out_dir + "workup/fastqs/{experiment}.part_{splitid}.bowtie2.unmapped_R2.fq.gz"],
                    experiment=ALL_EXPERIMENTS, splitid=NUM_CHUNKS)

STAR_ALIGN = expand(out_dir + "workup/alignments/{experiment}.part_{splitid}.Aligned.out.sorted.bam", experiment=ALL_EXPERIMENTS, splitid=NUM_CHUNKS)

UNALIGNED = expand([out_dir + "workup/unmapped/{experiment}_R1.part_{splitid}.unaligned.fastq.gz",
                    out_dir + "workup/unmapped/{experiment}_R2.part_{splitid}.unaligned.fastq.gz"],
                    experiment=ALL_EXPERIMENTS, splitid=NUM_CHUNKS)

MERGE_RNA = expand(out_dir + "workup/alignments/{experiment}.merged.RPM.bam", experiment=ALL_EXPERIMENTS)

CHR_RPM = expand([out_dir + "workup/alignments/{experiment}.part_{splitid}.Aligned.out.sorted.chr.bam",
                  out_dir + "workup/alignments/{experiment}.part_{splitid}.bowtie2.sorted.mapped.chr.bam"], 
                  experiment=ALL_EXPERIMENTS, splitid=NUM_CHUNKS)

################################################################################
#Bead workup
################################################################################

FQ_TO_BAM = expand(out_dir + "workup/alignments/{experiment}.part_{splitid}.BPM.bam", experiment=ALL_EXPERIMENTS, splitid=NUM_CHUNKS)

MERGE_BEAD = expand(out_dir + "workup/alignments/{experiment}.merged.BPM.bam", experiment=ALL_EXPERIMENTS)

################################################################################
#Clustering
################################################################################

CLUSTERS = expand(out_dir + "workup/clusters/{experiment}.part_{splitid}.clusters", experiment=ALL_EXPERIMENTS, splitid=NUM_CHUNKS)
CLUSTERS_MERGED = expand(out_dir + "workup/clusters/{experiment}.clusters", experiment=ALL_EXPERIMENTS)
CLUSTERS_MERGED_COMPLETE = expand(out_dir + "workup/clusters/{experiment}.complete.clusters", experiment=ALL_EXPERIMENTS)
CONDITION_CLUSTERS = expand(
    out_dir + "workup/condition-clusters/{experiment}.{condition}.clusters", 
    experiment=ALL_EXPERIMENTS, 
    condition=conditions
)

##############################################################################
# Post Clustering
##############################################################################

THRESH_AND_SPLIT_NO_CONDITION = expand(
        [out_dir + "workup/splitbams-all-conditions/{experiment}.done"], 
        experiment=ALL_EXPERIMENTS
)

THRESH_AND_SPLIT_CONDITION = expand(
        [out_dir + "workup/splitbams-by-condition/{experiment}.{condition}.done"], 
        experiment=ALL_EXPERIMENTS, 
        condition=conditions
)

#COUNTS = [out_dir + "workup/clusters/cluster_statistics.txt"]
#
#SIZES = [out_dir + "workup/condition-clusters/DPM_read_distribution.pdf",
#         out_dir + "workup/condition-clusters/DPM_cluster_distribution.pdf",
#         out_dir + "workup/condition-clusters/BPM_cluster_distribution.pdf",
#         out_dir + "workup/condition-clusters/BPM_read_distribution.pdf"]
#
#ECDFS = [out_dir + "workup/clusters/Max_representation_ecdf.pdf",
#         out_dir + "workup/clusters/Max_representation_counts.pdf"]
#
#MULTI_QC = [out_dir + "workup/qc/multiqc_report.html"]
 
################################################################################
################################################################################
#RULE ALL
################################################################################
################################################################################

rule all:
    input: CONFIG + SPLIT_FQ + ALL_FASTQ + TRIM + TRIM_LOG + TRIM_RD + BARCODEID + LE_LOG_ALL + MERGE_BEAD + FQ_TO_BAM + CLUSTERS + CLUSTERS_MERGED_COMPLETE + CONDITION_CLUSTERS + SPLIT_RPM_BPM + SPLIT_RPM_BPM2 + BT2_RNA_ALIGN + STAR_ALIGN + CHR_RPM + MERGE_RNA + THRESH_AND_SPLIT_CONDITION + THRESH_AND_SPLIT_NO_CONDITION 

#Send and email if an error occurs during execution
onerror:
    shell('mail -s "an error occurred" ' + email + ' < {log}')

wildcard_constraints:
    experiment = "[^\.]+"


##################################################################################
#Trimming and barcode identification
##################################################################################

rule splitfq:
    input:
        r1 = lambda wildcards: FILES[wildcards.experiment]['R1'],
        r2 = lambda wildcards: FILES[wildcards.experiment]['R2']
    output:
        temp(expand([(out_dir + "workup/splitfq/{{experiment}}_R1.part_{splitid}.fastq"), (out_dir + "workup/splitfq/{{experiment}}_R2.part_{splitid}.fastq")],  splitid=NUM_CHUNKS))
    params:
        dir = out_dir + "workup/splitfq",
        prefix_r1 = "{experiment}_R1.part_0",
        prefix_r2 = "{experiment}_R2.part_0"
    log:
        out_dir + "workup/logs/{experiment}.splitfq.log"
    conda:
        "envs/sprite.yaml"
    threads: 
        8
    benchmark:
        "benchmarks/{experiment}.splitfq.tsv"
    shell:
        '''
        mkdir -p {params.dir}
        bash {split_fastq} {input.r1} {num_chunks} {params.dir} {params.prefix_r1}
        bash {split_fastq} {input.r2} {num_chunks} {params.dir} {params.prefix_r2}
        '''

rule compress_fastq:
    input:
        r1 = out_dir + "workup/splitfq/{experiment}_R1.part_{splitid}.fastq",
        r2 = out_dir + "workup/splitfq/{experiment}_R2.part_{splitid}.fastq"
    output:
        r1 = out_dir + "workup/splitfq/{experiment}_R1.part_{splitid}.fastq.gz",
        r2 = out_dir + "workup/splitfq/{experiment}_R2.part_{splitid}.fastq.gz"
    conda:
        "envs/sprite.yaml"
    threads:
        8
    benchmark:
        "benchmarks/{experiment}.{splitid}.compress_fastq.tsv"
    shell:
        '''
        pigz -p {threads} {input.r1}
        pigz -p {threads} {input.r2}
        '''        

#Trim adaptors
#multiple cores requires pigz to be installed on the system
rule adaptor_trimming_pe:
    input:
        [out_dir + "workup/splitfq/{experiment}_R1.part_{splitid}.fastq.gz", 
         out_dir + "workup/splitfq/{experiment}_R2.part_{splitid}.fastq.gz"]
    output:
         out_dir + "workup/trimmed/{experiment}_R1.part_{splitid}_val_1.fq.gz",
         out_dir + "workup/trimmed/{experiment}_R1.part_{splitid}.fastq.gz_trimming_report.txt",
         out_dir + "workup/trimmed/{experiment}_R2.part_{splitid}_val_2.fq.gz",
         out_dir + "workup/trimmed/{experiment}_R2.part_{splitid}.fastq.gz_trimming_report.txt"
    threads:
        10
    log:
        out_dir + "workup/logs/{experiment}.{splitid}.trim_galore.log"
    conda:
        "envs/sprite.yaml"
    benchmark:
        "benchmarks/{experiment}.{splitid}.adaptor_trimming_pe.tsv"
    shell:
        '''
        trim_galore \
        --paired \
        --gzip \
        --cores {threads} \
        --quality 20 \
        --fastqc \
        -o {out_dir}workup/trimmed/ \
        {input} &> {log}
        '''

#Identify barcodes using BarcodeIdentification_v1.2.0.jar
rule barcode_id:
    input:
        r1 = out_dir + "workup/trimmed/{experiment}_R1.part_{splitid}_val_1.fq.gz",
        r2 = out_dir + "workup/trimmed/{experiment}_R2.part_{splitid}_val_2.fq.gz"
    output:
        r1_barcoded = out_dir + "workup/fastqs/{experiment}_R1.part_{splitid}.barcoded.fastq.gz",
        r2_barcoded = out_dir + "workup/fastqs/{experiment}_R2.part_{splitid}.barcoded.fastq.gz"
    log:
        out_dir + "workup/logs/{experiment}.{splitid}.bID.log"
    benchmark:
        "benchmarks/{experiment}.{splitid}.barcode_id.tsv"
    shell:
        "java -jar {barcode_id_jar} \
        --input1 {input.r1} --input2 {input.r2} \
        --output1 {output.r1_barcoded} --output2 {output.r2_barcoded} \
        --config {bid_config} &> {log}"

#Get ligation efficiency
rule get_ligation_efficiency:
    input:
        r1 = out_dir + "workup/fastqs/{experiment}_R1.part_{splitid}.barcoded.fastq.gz" 
    output:
        temp(out_dir + "workup/{experiment}.part_{splitid}.ligation_efficiency.txt")
    conda:
        "envs/sprite.yaml"
    shell:
        "python {lig_eff} {input.r1} > {output}"


rule cat_ligation_efficiency:
    input:
        expand(out_dir + "workup/{experiment}.part_{splitid}.ligation_efficiency.txt", experiment=ALL_EXPERIMENTS, splitid=NUM_CHUNKS)
    output:
        out_dir + "workup/ligation_efficiency.txt"
    shell:
        "tail -n +1 {input} > {output}"

rule split_bpm_rpm:
    '''
    split bpm and rpm will also remove incomplete barcodes
    '''
    input:
        out_dir + "workup/fastqs/{experiment}_R1.part_{splitid}.barcoded.fastq.gz"
    output:
        out_dir + "workup/fastqs/{experiment}_R1.part_{splitid}.barcoded_rpm.fastq.gz",
        out_dir + "workup/fastqs/{experiment}_R1.part_{splitid}.barcoded_bpm.fastq.gz",
        out_dir + "workup/fastqs/{experiment}_R1.part_{splitid}.barcoded_short.fastq.gz"
    log:
        out_dir + "workup/logs/{experiment}.{splitid}.BPM_RPM.log"
    conda:
       "envs/sprite.yaml"
    benchmark:
        "benchmarks/{experiment}.{splitid}.split_bpm_rpm.tsv"
    shell:
        "python {split_bpm_rpm} --r1 {input} &> {log}"

rule split_bpm_rpm2:
    '''
    split bpm and rpm will also remove incomplete barcodes
    '''
    input:
        out_dir + "workup/fastqs/{experiment}_R2.part_{splitid}.barcoded.fastq.gz"
    output:
        out_dir + "workup/fastqs/{experiment}_R2.part_{splitid}.barcoded_rpm.fastq.gz",
        out_dir + "workup/fastqs/{experiment}_R2.part_{splitid}.barcoded_bpm.fastq.gz",
        out_dir + "workup/fastqs/{experiment}_R2.part_{splitid}.barcoded_short.fastq.gz"
    log:
        out_dir + "workup/logs/{experiment}.{splitid}.BPM_RPM.log"
    conda:
       "envs/sprite.yaml"
    benchmark:
        "benchmarks/{experiment}.{splitid}.split_bpm_rpm2.tsv"
    shell:
        "python {split_bpm_rpm} --r1 {input} &> {log}"

################################################################################
#Cutadapt
################################################################################

rule cutadapt_rpm:
    input:
        read1 = out_dir + "workup/fastqs/{experiment}_R1.part_{splitid}.barcoded_rpm.fastq.gz",
        read2 = out_dir + "workup/fastqs/{experiment}_R2.part_{splitid}.barcoded_rpm.fastq.gz"
    output:
        r1=out_dir + "workup/trimmed/{experiment}_R1.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz",
        r2=out_dir + "workup/trimmed/{experiment}_R2.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz",
        qc=out_dir + "workup/trimmed/{experiment}_R1.part_{splitid}.barcoded_rpm.RDtrim.qc.txt"
    params:
        adapters_r1 = "-a ATCAGCACTTAGCGTCAG",
        adapters_r2 = "-G CTGACGCTAAGTGCTGAT",
        others = "--minimum-length 20"
    log:
        out_dir + "workup/logs/{experiment}.{splitid}.RPM.cutadapt.log"
    threads: 
        10
    conda:
        "envs/sprite.yaml"
    benchmark:
        "benchmarks/{experiment}.{splitid}.cutadapt_rpm.tsv"
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
        
         fastqc {output.r1}
         fastqc {output.r2}
        '''

rule cutadapt_oligo:
    '''
    Trim 9mer oligo sequence from bead barcode
    '''
    input:
        out_dir + "workup/fastqs/{experiment}_R1.part_{splitid}.barcoded_bpm.fastq.gz"
    output:
        fastq=out_dir + "workup/trimmed/{experiment}_R1.part_{splitid}.barcoded_bpm.RDtrim.fastq.gz",
        qc=out_dir + "workup/trimmed/{experiment}_R1.part_{splitid}.barcoded_bpm.RDtrim.qc.txt"
    params:
        adapters_r1 = oligos
    log:
        out_dir + "workup/logs/{experiment}.{splitid}.BPM.cutadapt.log"
    threads: 
        10
    conda:
        "envs/sprite.yaml"
    benchmark:
        "benchmarks/{experiment}.{splitid}.cutadapt_oligo.tsv"
    shell:
        '''
        (cutadapt \
         {params.adapters_r1} \
         -o {output.fastq} \
         -j {threads} \
         {input} > {output.qc}) &> {log}
        '''

################################################################################
#RNA alignment
################################################################################

rule bowtie2_align:
    input:
        fq1 = out_dir + "workup/trimmed/{experiment}_R1.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz",
        fq2 = out_dir + "workup/trimmed/{experiment}_R2.part_{splitid}.barcoded_rpm.RDtrim.fastq.gz"
    output:
        bam = temp(out_dir + "workup/alignments/{experiment}.part_{splitid}.bowtie2.bam"),
        mapped = out_dir + "workup/alignments/{experiment}.part_{splitid}.bowtie2.sorted.mapped.bam",
        unmapped = out_dir + "workup/alignments/{experiment}.part_{splitid}.bowtie2.sorted.unmapped.bam",
        index = out_dir + "workup/alignments/{experiment}.part_{splitid}.bowtie2.sorted.mapped.bam.bai"
    log:
        out_dir + "workup/logs/{experiment}.{splitid}.bt2.log"
    threads: 
        10
    conda:
        "envs/sprite.yaml"
    benchmark:
        "benchmarks/{experiment}.{splitid}.bowtie2_align.tsv"
    shell:
        '''
        (bowtie2 \
        -p 10 \
        -t \
        -x {bowtie2_index} \
        -1 {input.fq1} -2 {input.fq2} | \
        samtools view -bS -> {output.bam}) &> {log}
        samtools view -b -f 4 {output.bam} | samtools sort -n -o {output.unmapped}
        samtools view -b -F 4 {output.bam} | samtools sort > {output.mapped}
        samtools index {output.mapped}
        '''

rule bam_to_fq:
    input: out_dir + "workup/alignments/{experiment}.part_{splitid}.bowtie2.sorted.unmapped.bam"
    output: 
        r1 = out_dir + "workup/fastqs/{experiment}.part_{splitid}.bowtie2.unmapped_R1.fq.gz",
        r2 = out_dir + "workup/fastqs/{experiment}.part_{splitid}.bowtie2.unmapped_R2.fq.gz"
    log:
        out_dir + "workup/logs/{experiment}.{splitid}.bam2fq.log"
    threads: 
        1
    conda:
        "envs/sprite.yaml"
    benchmark:
        "benchmarks/{experiment}.{splitid}.bam_to_fq.tsv"
    shell:
        '''
        samtools fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -@ {threads} -n {input} 
        '''

rule star_align:
    input:
        r1 = out_dir + "workup/fastqs/{experiment}.part_{splitid}.bowtie2.unmapped_R1.fq.gz",
        r2 = out_dir + "workup/fastqs/{experiment}.part_{splitid}.bowtie2.unmapped_R2.fq.gz"
    output:
        sam = temp(out_dir + "workup/alignments/{experiment}.part_{splitid}.Aligned.out.sam"),
        sorted = out_dir + "workup/alignments/{experiment}.part_{splitid}.Aligned.out.sorted.bam",
        filtered = temp(out_dir + "workup/alignments/{experiment}.part_{splitid}.Aligned.out.bam")
    params:
        STAR_OPTIONS = "--readFilesCommand zcat --alignEndsType EndToEnd --outFilterScoreMin 10 --outFilterMultimapNmax 1 --outFilterMismatchNmax 10 --alignIntronMax 100000 --alignMatesGapMax 1300 --alignIntronMin 80 --alignSJDBoverhangMin 5 --alignSJoverhangMin 8 --chimSegmentMin 20 --alignSJstitchMismatchNmax 5 -1 5 5 --outSAMunmapped Within --outReadsUnmapped Fastx", prefix = out_dir + "workup/alignments/{experiment}.part_{splitid}."
    log:
        out_dir + "workup/logs/{experiment}.{splitid}.star.log"
    threads:
        10
    conda:
        "envs/sprite.yaml"
    benchmark:
        "benchmarks/{experiment}.{splitid}.star_align.tsv"
    shell:
        '''
        (STAR \
        --genomeDir {star_index} \
        --readFilesIn {input.r1} {input.r2} \
        --runThreadN {threads} {params.STAR_OPTIONS} \
        --outFileNamePrefix {params.prefix}) &> {log}

        samtools view -@ {threads} -bS {output.sam} > {output.filtered}
        samtools sort -@ {threads} -o {output.sorted} {output.filtered}
        samtools index {output.sorted}
        '''

rule compress_unaligned:
    input:
        r1 = out_dir + "workup/alignments/{experiment}.part_{splitid}.Unmapped.out.mate1",
        r2 = out_dir + "workup/alignments/{experiment}.part_{splitid}.Unmapped.out.mate2"
    output:
        fq1 = out_dir + "workup/unmapped/{experiment}_R1.part_{splitid}.unaligned.fastq.gz",
        fq2 = out_dir + "workup/unmapped/{experiment}_R2.part_{splitid}.unaligned.fastq.gz"
    conda:
        "envs/sprite.yaml"
    threads:
        8
    benchmark:
        "benchmarks/{experiment}.{splitid}.compress_unaligned.tsv"
    shell:
        '''
        pigz -p {threads} {input.r1} {input.r2}
       
        mv {input.r1}.gz {output.fq1}
        mv {input.r2}.gz {output.fq2}
        '''

rule add_chr:
    input:
        star = out_dir + "workup/alignments/{experiment}.part_{splitid}.Aligned.out.sorted.bam",
        bt2 = out_dir + "workup/alignments/{experiment}.part_{splitid}.bowtie2.sorted.mapped.bam"
    output:
        star = out_dir + "workup/alignments/{experiment}.part_{splitid}.Aligned.out.sorted.chr.bam",
        bt2 = out_dir + "workup/alignments/{experiment}.part_{splitid}.bowtie2.sorted.mapped.chr.bam"
    log:
        out_dir + "workup/logs/{experiment}.{splitid}.add_chr.log",
    conda:
        "envs/sprite.yaml"
    benchmark:
        "benchmarks/{experiment}.{splitid}.add_chr.tsv"
    shell:
        '''
        python {add_chr} -i {input.star} -o {output.star} --assembly {assembly} &> {log}
        python {add_chr_bt2} -i {input.bt2} -o {output.bt2} --assembly {assembly} 
        '''

rule merge_rna:
    input:
        bt2 = expand(out_dir + "workup/alignments/{{experiment}}.part_{splitid}.bowtie2.sorted.mapped.chr.bam", splitid=NUM_CHUNKS),
        star = expand(out_dir + "workup/alignments/{{experiment}}.part_{splitid}.Aligned.out.sorted.chr.bam", splitid=NUM_CHUNKS)
    output:
        out_dir + "workup/alignments/{experiment}.merged.RPM.bam"
    conda:
        "envs/sprite.yaml"
    threads:
        8
    log:
        out_dir + "workup/logs/{experiment}.merge_bams.log"
    benchmark:
        "benchmarks/{experiment}.merge_rna.tsv"
    shell:
        '''
        (samtools merge -@ {threads} {output} {input.bt2} {input.star}) &> {log}
        samtools index {output}
        '''

##############################################################################
#Workup Bead Oligo
##############################################################################

rule fastq_to_bam:
    input:
        out_dir + "workup/trimmed/{experiment}_R1.part_{splitid}.barcoded_bpm.RDtrim.fastq.gz"
    output:
        sorted = out_dir + "workup/alignments/{experiment}.part_{splitid}.BPM.bam",
        bam = temp(out_dir + "workup/alignments/{experiment}.part_{splitid}.BPM.unsorted.bam")
    log:
        out_dir + "workup/logs/{experiment}.{splitid}.make_bam.log"
    conda:
        "envs/sprite.yaml"
    threads:
        8
    benchmark:
        "benchmarks/{experiment}.{splitid}.fastq_to_bam.tsv"
    shell:
        '''
        python {fq_to_bam} --input {input} --output {output.bam} --config {bid_config} &> {log}
        samtools sort -@ {threads} -o {output.sorted} {output.bam}

        '''

rule merge_beads:
    input:
        expand(out_dir + "workup/alignments/{{experiment}}.part_{splitid}.BPM.bam", splitid = NUM_CHUNKS)
    output:
        out_dir + "workup/alignments/{experiment}.merged.BPM.bam"
    conda:
        "envs/sprite.yaml"
    log:
        out_dir + "workup/logs/{experiment}.merge_beads.log"
    threads:
        8
    benchmark:
        "benchmarks/{experiment}.merge_beads.tsv"
    shell:
        '''
        (samtools merge -@ {threads} {output} {input}) >& {log}
        '''

###########################################################################
#Make clusters
###########################################################################

rule make_clusters:
    input:
        rpm=out_dir + "workup/alignments/{experiment}.part_{splitid}.Aligned.out.sorted.chr.bam",
        bpm=out_dir + "workup/alignments/{experiment}.part_{splitid}.BPM.bam",
        bt2=out_dir + "workup/alignments/{experiment}.part_{splitid}.bowtie2.sorted.mapped.chr.bam"
    output:
        unsorted = temp(out_dir + "workup/clusters/{experiment}.part_{splitid}.unsorted.clusters"),
        sorted = out_dir + "workup/clusters/{experiment}.part_{splitid}.clusters"
    log:
        out_dir + "workup/logs/{experiment}.{splitid}.make_clusters.log"
    conda:
        "envs/sprite.yaml"
    benchmark:
        "benchmarks/{experiment}.{splitid}.make_clusters.tsv"
    shell:
        '''
        (python {get_clusters} \
        -i {input.bpm} {input.bt2} {input.rpm}\
        -o {output.unsorted} \
        -n {num_tags})  &> {log}

        sort -k 1 -T {temp_dir} {output.unsorted} > {output.sorted}
        '''

rule merge_clusters:
    input:
        expand(out_dir + "workup/clusters/{{experiment}}.part_{splitid}.clusters", splitid=NUM_CHUNKS)
    output:
        mega = temp(out_dir + "workup/clusters/{experiment}.duplicated.clusters"),
        final = out_dir + "workup/clusters/{experiment}.clusters"
    conda:
       "envs/sprite.yaml"
    log:
        out_dir + "workup/logs/{experiment}.merge_clusters.log"
    benchmark:
        "benchmarks/{experiment}.merge_clusters.tsv"
    shell:
        '''
         sort -k 1 -T {temp_dir} -m {input} > {output.mega}
        (python {merge_clusters} -i {output.mega} -o {output.final}) &> {log}
        '''        

##############################################################################
# Profile clusters
##############################################################################

# Generate simple statistics for clusters
rule generate_cluster_statistics:
    input:
        expand([out_dir + "workup/clusters/{experiment}.complete.clusters"], experiment=ALL_EXPERIMENTS)
    output:
        out_dir + "workup/clusters/cluster_statistics.txt"
    params:
        dir = out_dir + "workup/clusters"
    conda:
        "envs/sprite.yaml"
    shell:
        '''
        python {cluster_counts} --directory {params.dir} --pattern .clusters > {output}
        '''

# Generate ecdfs of oligo distribution
#rule generate_cluster_ecdfs:
#    input:
#        expand([out_dir + "workup/clusters/{experiment}.clusters"], experiment=ALL_EXPERIMENTS)
#    output:
#        ecdf = out_dir + "workup/clusters/Max_representation_ecdf.pdf",
#        counts = out_dir + "workup/clusters/Max_representation_counts.pdf"
#    params:
#        dir = out_dir + "workup/clusters"
#    conda:
#        "envs/plotting.yaml"
#    shell:
#        '''
#        python {cluster_ecdfs} --directory {params.dir} --pattern .clusters --xlim 30
#        '''
#
## Profile size distribution of clusters
#rule get_size_distribution:
#    input:
#        expand([out_dir + "workup/condition-clusters/{experiment}.{condition}.clusters"], experiment=ALL_EXPERIMENTS, condition=conditions),
#        expand([out_dir + "workup/clusters/{experiment}.complete.clusters"], experiment=ALL_EXPERIMENTS),
#    output:
#        # FIXME: should these be changed to RPM equivalents?
#        dpm = out_dir + "workup/condition-clusters/DPM_read_distribution.pdf",
#        dpm2 = out_dir + "workup/condition-clusters/DPM_cluster_distribution.pdf",
#        bpm = out_dir + "workup/condition-clusters/BPM_read_distribution.pdf",
#        bpm2 = out_dir + "workup/condition-clusters/BPM_cluster_distribution.pdf",
#        no_condition_dpm = out_dir + "workup/clusters/DPM_read_distribution.pdf",
#        no_condition_dpm2 = out_dir + "workup/clusters/DPM_cluster_distribution.pdf",
#        no_condition_bpm = out_dir + "workup/clusters/BPM_read_distribution.pdf",
#        no_condition_bpm2 = out_dir + "workup/clusters/BPM_cluster_distribution.pdf"
#    params:
#        condition_dir = "workup/condition-clusters",
#        no_condition_dir = "workup/clusters"
#    conda:
#        "envs/sprite.yaml"
#    shell:
#        '''
#        python {cluster_sizes} --directory {params.no_condition_dir} --pattern .clusters --readtype BPM
#        python {cluster_sizes} --directory {params.no_condition_dir} --pattern .clusters --readtype DPM
#
#        python {cluster_sizes} --directory {params.condition_dir} --pattern .clusters --readtype BPM
#        python {cluster_sizes} --directory {params.condition_dir} --pattern .clusters --readtype DPM
#        '''
#
################################################################################
# Logging and MultiQC
################################################################################
# Copy config.yaml into logs folder with run date
rule log_config:
    input:
        config_path
    output:
        out_dir + "workup/logs/config_" + run_date + "yaml"
    shell:
        "cp {input} {output}"

#rule multiqc:
#    input:
#        expand([out_dir + "workup/clusters/{experiment}.complete.clusters"], experiment=ALL_EXPERIMENTS, condition=conditions),
#        expand([out_dir + "workup/condition-clusters/{experiment}.{condition}.clusters"], experiment=ALL_EXPERIMENTS, condition=conditions)
#    output:
#        out_dir + "workup/qc/multiqc_report.html"
#    log:
#        out_dir + "workup/logs/multiqc.log"
#    conda:
#        "envs/sprite.yaml"
#    shell:
#        "(multiqc --force {out_dir}workup -o {out_dir}workup/qc) &> {log}"

##############################################################################
# Remove incorrect clusters
##############################################################################
rule split_incorrect_clusters:
    input:
        clusters = out_dir + "workup/clusters/{experiment}.clusters"
    output:
        complete_clusters = out_dir + "workup/clusters/{experiment}.complete.clusters",
        incomplete_clusters = out_dir + "workup/clusters/{experiment}.incomplete.clusters"
    conda:
        "envs/sprite.yaml"
    log:
        out_dir + "workup/logs/{experiment}.incorrectclusters.log"
    benchmark:
        "benchmarks/{experiment}.split_incorrect_clusters.tsv"
    shell:
        '''
        (python {split_incorrect_clusters} \
        --clusters {input.clusters} \
        --complete_output {output.complete_clusters} \
        --incomplete_output {output.incomplete_clusters} \
    	--format {rounds_format}) &> {log}
        '''

##############################################################################
# Split based on first tag
##############################################################################
rule split_on_first_tag:
    input:
        complete_clusters = out_dir + "workup/clusters/{experiment}.complete.clusters"
    output:
        expand(
            out_dir + "workup/condition-clusters/{{experiment}}.{condition}.clusters",
            condition=conditions
        )
    conda:
        "envs/sprite.yaml"
    log:
        out_dir + "workup/logs/{experiment}.splitonfirsttag.log"
    benchmark:
        "benchmarks/{experiment}.split_on_first_tag.tsv"
    shell:
        '''
        (python {split_on_first_tag} \
            --complete_clusters {input.complete_clusters} \
            --output_dir workup/condition-clusters) &> {log}
        '''

##############################################################################
# Splitbams
##############################################################################

# Generate bam files for individual targets based on assignments from clusterfile
rule thresh_and_split_condition:
    input:
        bam = out_dir + "workup/alignments/{experiment}.merged.RPM.bam",
        clusters = out_dir + "workup/condition-clusters/{experiment}.{condition}.clusters"
    output:
        bam = out_dir + "workup/splitbams-by-condition/{experiment}.{condition}.bam",
        touch = touch(out_dir + "workup/splitbams-by-condition/{experiment}.{condition}.done")
    conda:
        "envs/sprite.yaml"
    log:
        out_dir + "workup/logs/{experiment}.{condition}.splitbams.log"
    benchmark:
        out_dir + "benchmarks/{experiment}.{condition}.thresh_and_split_control.tsv"
    params:
        directory = "workup/splitbams-by-condition"
    shell:
        '''
        (python {tag_and_split} \
            -i {input.bam} \
            -c {input.clusters} \
            -o {output.bam} \
            -d {params.directory} \
            --min_oligos {min_oligos} \
            --proportion {proportion} \
            --max_size {max_size} \
            --num_tags {num_tags}) &> {log}
        '''

# Generate bam files for individual targets based on assignments from clusterfile
rule thresh_and_split_no_condition:
    input:
        bam = out_dir + "workup/alignments/{experiment}.merged.RPM.bam",
        clusters = out_dir + "workup/clusters/{experiment}.complete.clusters"
    output:
        bam = out_dir + "workup/splitbams-all-conditions/{experiment}.ALL_CONDITIONS.bam",
        touch = touch(out_dir + "workup/splitbams-all-conditions/{experiment}.done")
    conda:
        "envs/sprite.yaml"
    log:
        out_dir + "workup/logs/{experiment}.merged.splitbams.log"
    benchmark:
        out_dir + "benchmarks/{experiment}.merged.thresh_and_split_control.tsv"
    params:
        directory = "workup/splitbams-all-conditions"
    shell:
        '''
        (python {tag_and_split} \
            -i {input.bam} \
            -c {input.clusters} \
            -o {output.bam} \
            -d {params.directory} \
            --min_oligos {min_oligos} \
            --proportion {proportion} \
            --max_size {max_size} \
            --num_tags {num_tags}) &> {log}
        '''

# Generate summary statistics of individiual bam files
rule generate_splitbam_statistics:
    input:
        expand([out_dir + "workup/splitbams-all-conditions/{experiment}.done"], experiment=ALL_EXPERIMENTS),
        expand([out_dir + "workup/splitbams-by-condition/{experiment}.done"], experiment=ALL_EXPERIMENTS)
    output:
        all_conditions = out_dir + "workup/splitbams-all-conditions/splitbam_statistics.txt",
        by_condition = out_dir + "workup/splitbams-by-condition/splitbam_statistics.txt"
    params:
        all_conditions = out_dir + "workup/splitbams-all-conditions",
        by_condition = out_dir + "workup/splitbams-by-condition"
    conda:
        "envs/sprite.yaml"
    shell:
        """
        for f in {params.all_conditions}/*bam; do echo $f; samtools view -c $f; done > {output.all_conditions}
        for f in {params.by_condition}/*bam; do echo $f; samtools view -c $f; done > {output.by_condition}
        """

