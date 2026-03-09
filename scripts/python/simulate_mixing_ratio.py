from collections import Counter, defaultdict
import click
from scipy.stats import bernoulli
from typing import Tuple
import numpy as np


def get_total_counts(counts_dict, print_values=True) -> Tuple[int, int]:
    total_bpm = 0
    total_rpm = 0
    for barcode in counts_dict.keys():
        if "B" in counts_dict[barcode]:
            total_bpm += counts_dict[barcode]["B"]
        if "R" in counts_dict[barcode]:
            total_rpm += counts_dict[barcode]["R"]

    print(f"Total BPM: {total_bpm}")
    print(f"Total RPM: {total_rpm}")
    print(f"BPM/RPM: {total_bpm / total_rpm}")

    return total_bpm, total_rpm


@click.command()
@click.option('--input', '-i', 'input_file', required=True, help='Input file path (CSV/TSV/GZ).')
@click.option('--output', '-o', 'output_file', required=True, help='Output file path.')
@click.option('--cdna_fraction', '-c', 'cdna_fraction', type=float,
              required=True, help='Proportion of original cDNA reads to keep')
@click.option('--bead_fraction', '-b', 'bead_fraction', type=float,
              required=True, help='Proportion of original bead reads to keep')
def main(input_file, output_file, cdna_fraction, bead_fraction):
    default_rng = np.random.default_rng(6)

    # FIXME: Support arbitrary mixing ratios, we want any combination of bead fractions and cdna fractions from 0 to 1
    assert bead_fraction == 1

    counts_dict = defaultdict(Counter)
    cluster_data = defaultdict(list)

    print(f"Reading {input_file}")
    with open(input_file) as f:
        all_lines = [line for line in f]
        for line in all_lines:
            parts = line.strip().split("\t")
            barcode = ".".join(parts[0].split(".")[:-1])
            tags = parts[1:]
            tag_types = [part[0] for part in parts[1:]]
            counts_dict[barcode].update(tag_types)
            cluster_data[barcode] = tags

    total_bpm, total_rpm = get_total_counts(counts_dict)

    # FIXME: Support arbitrary mixing ratios, right now this assumes that we're only removing cDNA reads
    omit_cdna_read = bernoulli.rvs(cdna_fraction, loc=0, size=total_rpm, random_state=default_rng)

    subsampled_counts_dict = defaultdict(Counter)
    subsampled_cluster_data = defaultdict(list)

    cdna_read_counter = 0
    barcodes = list(cluster_data.keys())

    print("Subsampling")
    for barcode in barcodes:
        reads = cluster_data[barcode]
        bpm_reads = [read for read in reads if read.startswith("B")]
        rpm_reads = [read for read in reads if read.startswith("R")]
        subsampled_rpm_reads = []

        for read in rpm_reads:
            if omit_cdna_read[cdna_read_counter] == 1:
                subsampled_rpm_reads.append(read)

            cdna_read_counter += 1

        updated_reads = bpm_reads + subsampled_rpm_reads
        updated_read_types = [read[0] for read in updated_reads]
        subsampled_cluster_data[barcode] = updated_reads
        subsampled_counts_dict[barcode].update(updated_read_types)

    get_total_counts(subsampled_counts_dict)

    print(f"Writing {output_file}")
    with open(output_file, "w") as f:
        for barcode, reads in subsampled_counts_dict.items():
            reads_to_write = "\t".join(reads)
            f.write(barcode + "\t" + reads_to_write + "\n")

    print("Done!")


if __name__ == "__main__":
    main()
