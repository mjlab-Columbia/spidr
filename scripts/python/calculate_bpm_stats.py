import csv
import gzip
from collections import Counter, defaultdict

import click
from tqdm import tqdm


def open_text(path, mode="r"):
    if path.endswith(".gz"):
        if "w" in mode:
            return gzip.open(path, "wt", encoding="utf-8", newline="")
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def format_output_line(barcode, umi_counts):
    num_umi = len(umi_counts)
    total_bpm = sum(umi_counts.values())
    num_duplicates = total_bpm - num_umi
    raw_counts = "\t".join(f"{umi}_{count}" for umi, count in umi_counts.items())
    return "\t".join(
        [barcode, str(num_umi), str(total_bpm), str(num_duplicates), raw_counts]
    )


@click.command()
@click.option(
    "--input",
    "-i",
    type=click.Path(exists=True),
    help="Tab-separated file with 2 cols. Col 1: barcode, Col 2: UMI sequence",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(),
    help="Output file containing duplication rate per cluster",
)
def main(input, output):
    """
    Script for gathering duplication rate from a 2 column file containing full barcodes and their UMI sequences

    The output will have a variable number of columns but the first 4 columns will be (in order):
    * full barcode
    * number of unique umis associated with the full barcode
    * total number of umis associated with the full barcode
    * number of duplicates (total - unique)

    This is followed by each UMI and its count in the firm <UMI sequence>_<UMI count>

    Last Updated: October 24th, 2025
    """

    print("Loading barcodes and UMIs")
    barcode_umis = defaultdict(Counter)

    with open_text(input) as file_in:
        reader = csv.reader(file_in, delimiter="\t")
        for row in tqdm(reader, desc="Counting UMIs per barcode", unit=" reads"):
            if len(row) < 2:
                continue
            barcode_umis[row[0]][row[1]] += 1

    with open_text(output, "w") as file_out:
        for barcode in tqdm(
            sorted(barcode_umis),
            total=len(barcode_umis),
            desc="Writing to output",
        ):
            file_out.write(format_output_line(barcode, barcode_umis[barcode]) + "\n")


if __name__ == "__main__":
    main()
