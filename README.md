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

