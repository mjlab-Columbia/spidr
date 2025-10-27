import pandas as pd
import click
from tqdm import tqdm
from collections import defaultdict
from pdb import set_trace
from typing import List


@click.command()
@click.option("--input", "-i", type=click.Path(exists=True), help="Tab-separated file with 2 cols. Col 1: barcode, Col 2: UMI sequence")
@click.option("--output", "-o", type=click.Path(), help="Output file containing duplication rate per cluster")
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
    if input.endswith("gz"):
        df = pd.read_csv(input, sep="\t", header=None, names=["barcode", "umi"], compression="gzip").set_index("barcode")
    else:
        df = pd.read_csv(input, sep="\t", header=None, names=["barcode", "umi"]).set_index("barcode")

    print("Sorting barcodes to speed up UMI search")
    df.sort_index(inplace=True)
    unique_barcodes = df.index.unique().tolist()
    input_progress_bar = tqdm(unique_barcodes, total=len(unique_barcodes), desc="Counting UMIs per barcode")
    barcode_dict = dict()
    output_entries = []

    for barcode in input_progress_bar:
        umis_array = df.loc[barcode].values

        # Ensure numpy array properly flattens depending on number of items returned
        if len(umis_array.shape) > 1:
            umis: List[str] = umis_array.squeeze().tolist()
        else:
            umis: List[str] = umis_array.tolist()

        umi_dict = defaultdict(int)

        for umi in umis:
            umi_dict[umi] += 1

        # We want to keep track of the overall duplication rate, but also report the details per cluster
        barcode_dict[barcode] = dict(umi_dict)
        num_umi = len(umi_dict.keys())
        total_bpm = sum(umi_dict.values())
        raw_counts = "\t".join([f"{umi}_{count}" for umi, count in umi_dict.items()])

        entry = {
            "barcode": barcode,
            "num_umi": num_umi,
            "total_bpm": total_bpm,
            "num_duplicates": total_bpm - num_umi,
            "raw_counts": raw_counts
        }
        output_entries.append(entry)


    with open(output, "w") as file_out:
        output_progress_bar = tqdm(output_entries, total=len(output_entries), desc="Writing to output")

        for entry in output_progress_bar:
            entry_string_cast = [str(item) if type(item) is not str else item for item in entry.values()]
            line_to_write = "\t".join(entry_string_cast) + "\n"
            file_out.write(line_to_write)


    return


if __name__ == "__main__":
    main()
