# SPIDR Pipeline

Repository to host all the code for generating bam files from the fastq files of a SPIDR experiment. If you would like to contribute or make changes, please create a pull request and see the [Development](#development) section.

## Pipeline Usage
### Prerequisites
1. Ensure you have `mamba` installed.
   ```bash
   mamba --version
   ```
   If not installed, follow the installation instructions in the pipeline documentation and consult the `mamba` installation instructions [here](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) for more details.

2. Ensure you have `snakemake` installed in your base environment:
   ```bash
   snakemake --version
   ```
   If not installed:
   ```bash
   mamba install -c bioconda snakemake==7.32.4
   ```

### Running the Pipeline
The following steps are for example data available on Columbia's Ginsburg HPC. Please see `documentation/technical/mjlab_spidr_pipeline_v1.pdf` for more details on how to run the pipeline on your own data.

1. Create a directory for your run:
   ```bash
   mkdir -p /burg/mjlab/projects/spidr-runs/<UNI>/<DIRNAME>
   cd /burg/mjlab/projects/spidr-runs/<UNI>/<DIRNAME>
   ```

2. Clone the SPIDR repository:
   ```bash
   git clone https://github.com/mjlab-Columbia/spidr.git
   cd spidr
   ```

3. Create the `experiments.json` file:
   ```bash
   python scripts/python/fastq2json.py --fastq_dir <path_to_read_files>
   ```

4. Set up configuration files:
   - Copy example barcode configuration file:
     ```bash
     cp /burg/mjlab/projects/spidr-barcoding-files/config_6_rounds_mTOR.txt ./config_6_rounds_mTOR.txt
     ```
   - Copy barcode format file:
     ```bash
     cp /burg/mjlab/projects/spidr-barcoding-files/format_6_rounds_mTOR.txt ./format_6_rounds_mTOR.txt
     ```
   - Copy and edit config.yaml:
     ```bash
     cp /burg/mjlab/projects/spidr-barcoding-files/config.yaml ./config.yaml
     ```

5. Request resources for an interactive job:
   ```bash
   srun --pty -t 0-04:00 -A mjlab --mem=32G -N 1 -c 4 /bin/bash
   ```

6. Run the pipeline:
   - Dry run to check configuration:
     ```bash
     bash run.sh --dry_run
     ```
   - Run in foreground:
     ```bash
     bash run.sh
     ```
   - Or run in background:
     ```bash
     sbatch run.sh
     ```

### Debugging
You can generate a visual representation of the pipeline using:
```bash
snakemake --rulegraph | dot -Tpdf > rulegraph.pdf
```

For more detailed information about the pipeline, please refer to the technical documentation in `documentation/technical/mjlab_spidr_pipeline_v1.pdf`.

### QC Summary Table
Each processed library produces a summary table at `workup/qc/summary/{experiment}.qc_summary.tsv`. Metrics are listed as rows and libraries as columns. Read counts are formatted as natural numbers; values of 1,000 or greater use comma delimiters (for example, `1,234,567`).

To merge summaries from multiple libraries:

```bash
bash scripts/bash/merge_summary_tables.sh \
  --input workup/qc/summary/*.tsv \
  --output workup/qc/summary/merged_{experiment}.tsv
```

| Metric | Description |
| --- | --- |
| `raw_fastq_reads` | Total R1 reads in the raw FASTQ files before any processing. |
| `fastp_reads_before_filtering` | R1 reads reported by fastp before deduplication (matches `raw_fastq_reads`; fastp's top-level `summary.total_reads` counts R1 and R2 together). |
| `fastp_reads_after_filtering` | R1 reads retained by fastp after deduplication. With `--dedup`, this reflects duplicate removal rather than quality filtering. |
| `fastp_duplication_rate` | Duplicate rate estimated by fastp (0–1). |
| `fastp_q20_rate_after_filtering` | Fraction of bases with Phred quality ≥ 20 after fastp filtering. |
| `fastp_q30_rate_after_filtering` | Fraction of bases with Phred quality ≥ 30 after fastp filtering. |
| `fastp_gc_content_after_filtering` | GC content after fastp filtering (0–1). |
| `filtered_fastq_reads` | Total R1 reads remaining after adapter trimming with Trim Galore. |
| `pre_alignment_fully_barcoded_reads_r1` | Sum of fully barcoded R1 reads across all chunks before alignment. |
| `total_bpm_reads` | Total bead-processing-module (BPM) reads across all chunks. |
| `total_rpm_reads` | Total RNA-processing-module (RPM) reads across all chunks. |
| `bowtie2_total_reads` | Total read pairs aligned to the ncRNA Bowtie2 index across all chunks. |
| `bowtie2_overall_alignment_rate_pct` | Read-weighted overall Bowtie2 concordant alignment rate (%). |
| `star_input_reads` | Total read pairs sent to STAR genome alignment across all chunks. |
| `star_uniquely_mapped_reads` | Total reads uniquely mapped by STAR across all chunks. |
| `star_mapped_multiple_loci` | STAR reads mapped to multiple loci across all chunks. |
| `star_mapped_too_many_loci` | STAR reads mapped to too many loci across all chunks. |
| `star_unmapped_too_many_mismatches` | STAR reads unmapped because of too many mismatches across all chunks. |
| `star_unmapped_too_short` | STAR reads unmapped because they were too short across all chunks. |
| `star_unmapped_other` | STAR reads unmapped for other reasons across all chunks. |
| `star_chimeric_reads` | STAR chimeric reads across all chunks. |
| `star_uniquely_mapped_pct` | Percentage of STAR input reads that were uniquely mapped. |
| `post_alignment_unique_rpm_reads` | Unique RPM read names in the merged RNA BAM after Bowtie2 and STAR alignment. |
| `bpm_duplication_rate` | BPM UMI duplication rate (duplicate UMIs / total BPM reads). |
| `barcoded_reads_assigned_to_clusters` | RPM reads assigned to clusters (sum of RPM counts in cluster file). |
| `barcoded_reads_assigned_to_bams_{condition}` | Reads assigned to split BAM files for a given experimental condition, excluding ambiguous/none/uncertain assignments. |
| `thresh_split_total_reads_{condition}` | Total reads processed when splitting clusters into condition-specific outputs for a given condition. |

## Development
### VSCode Setup
1. Install the Remote-SSH extension in VSCode
2. Configure your SSH connection to the remote host (Ginsburg)
3. Open the project folder on the remote host
4. Install recommended extensions when prompted

### Changing Folders in Remote-SSH
After connecting to the remote host and opening a folder, you can change folders by:
1. Opening the Command Palette (F1 or Cmd/Ctrl + Shift + P)
2. Typing "Remote-SSH: Open Folder"
3. Selecting the new folder path on the remote host

