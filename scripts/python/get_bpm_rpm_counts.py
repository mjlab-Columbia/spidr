import click
from tqdm import tqdm
from collections import Counter
import pandas as pd


@click.command()
@click.option("--clusters", "-c", type=click.Path(exists=True), help="Cluster file")
@click.option("--output", "-o", type=click.Path(),
              help="Tab-delimited table with columns barcode, BPM counts, and RPM counts")
def main(clusters, output):
    counts_table = []
    with open(clusters, "r") as file_in:
        progress_bar = tqdm(file_in)

        for line in progress_bar:
            barcode, *reads = line.strip().split("\t")
            reads_type = [read[:3] for read in reads]
            reads_dict = dict(Counter(reads_type))
            entry = {
                "barcode": barcode,
                "BPM": reads_dict["BPM"] if "BPM" in reads_dict else 0,
                "RPM": reads_dict["RPM"] if "RPM" in reads_dict else 0
            }
            counts_table.append(entry)

    df = pd.DataFrame(counts_table)
    df.to_csv(output, sep="\t", index=False)

    return


if __name__ == "__main__":
    main()
