import os
import sys
import click
from pathlib import Path
from typing import Iterable
from tqdm import tqdm
from pdb import set_trace


@click.command()
@click.option('--input_dir', '-i', type=click.Path(exists=True))
@click.option('--output_dir', '-o', type=click.Path(exists=True))
@click.option('--downsampling_script', '-d', type=click.Path(exists=True))
@click.option('--use_sbatch', type=bool, default=True)
@click.option('--dry_run', type=bool, default=False)
def main(input_dir: Path, output_dir: Path, downsampling_script: Path, use_sbatch: bool, dry_run: bool) -> None:
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)
    downsampling_script = os.path.abspath(downsampling_script)

    # Create iterable of bam files with full paths
    bams: Iterable[Path] = os.listdir(input_dir)
    bams: Iterable[Path] = [os.path.join(input_dir, bam) for bam in bams]

    # When using sbatch, there's no need for progress bar since jobs are dispatched concurrently
    if use_sbatch:
        progress = enumerate(bams)
    else:
        progress = tqdm(enumerate(bams), total=len(bams))

    for index, bam in progress:
        # If using srun to dispatch jobs, specify parameters here
        if use_sbatch:
            preamble: str = "sbatch --account=mjlab --time=00:30:00 --mem=16GB --nodes=1 --cpus-per-task=1"
        else:
            preamble: str = ""

        input_bam: Path = bam
        control_bams: Iterable[Path] = [bam for i, bam in enumerate(bams) if i != index]
        control_bams: str = " ".join(control_bams)
        command: str = f"java -jar --enable-preview {downsampling_script} {input_bam} {control_bams} {output_dir}"

        # The --wrap option encloses `command` in a simple /bin/sh script
        full_command: str = f"{preamble} --wrap='{command}'".strip()
        print(full_command, "\n\n") if dry_run else os.system(full_command)


if __name__ == "__main__":
    main()
