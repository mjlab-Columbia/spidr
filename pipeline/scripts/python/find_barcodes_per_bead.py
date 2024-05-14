import click
import os
from pdb import set_trace as st
import re
from collections import defaultdict
from typing import List
import gzip
from tqdm import tqdm

@click.command()
@click.option('--output', '-o', help='Output file')
@click.option('--header', '-h', help='Header of the barcode sequence in the fastq file', default=False)
@click.argument('input', type=click.Path(exists=True), nargs=-1)
def main(input: List[os.PathLike], header: bool, output: os.PathLike):
    """
    Find all bead IDs for each barcode and write to a file.
    Each line contains period separated barcode and a comma separated list of bead IDs

    Args:
        input (List[os.PathLike]): List of compressed or uncompressed barcoded fastq files from BPM reads
        output (os.PathLike): Output file path

    Returns:
        None

    Author: Darvesh Gorhe <dsg2157@columbia.edu>
    """

    # Iterate through fastq files and extract barcodes
    for file in tqdm(input, total=len(input)):
        if file.endswith('.gz'):
            with gzip.open(file, 'r') as f:
                lines = [line.decode("utf-8").strip() for line in f.readlines()]
                headers = lines[0::4]
                barcodes = [header.split('::')[-1] for header in headers]
        else:
            with open(file, 'r') as f:
                lines = f.readlines()
                headers = lines[0::4]
                barcodes = [header.split('::')[-1] for header in headers]

    # Clean up and sort barcodes
    full_barcodes = [b.replace('[', '').strip()[:-1] for b in barcodes]
    full_barcodes = [b.split(']') for b in full_barcodes]
    full_barcodes = sorted(full_barcodes)

    # Create a dictionary with barcodes as keys and a set of bead IDs as values
    barcode_dict = defaultdict(lambda: defaultdict(int))

    # For each barcode sequence, add the coreesponding bead ID and increment its count
    for b in full_barcodes:
        barcode_string = ".".join(b[1:])
        bead_id = b[0]
        barcode_dict[barcode_string][bead_id] += 1
    
    # Write the dictionary to a file
    with open(output, 'w') as f:
        if header:
            f.write(f"barcode\tbead_id_names\tbead_id_counts\n")

        for k, v in barcode_dict.items():
            sorted_bead_ids = dict(sorted(v.items()))
            bead_ids = ",".join(sorted_bead_ids.keys())
            bead_counts = ",".join([str(i) for i in sorted_bead_ids.values()])
            f.write(f"{k}\t{bead_ids}\t{bead_counts}\n")


if __name__ == '__main__':
    main()
