import os
import pandas as pd
from tqdm import tqdm

def get_files():
    # List of annotatated TSV files produced by annotator
    spidr_outputs = [os.path.join("spidr-annotated", file) for file in os.listdir("spidr-annotated")]
    encode_outputs = [os.path.join("encode-annotated", file) for file in os.listdir("encode-annotated")]
    return spidr_outputs + encode_outputs


def get_percentages_by_region(files):
    # Column names based on annotator GitHub
    columns = ["chromosome", "start", "stop", "name", "intensity", "strand", "gene_id", "gene_name", "genic_region_type", "all_overlapping_annotation"]

    # List to store percent_by_type dataframes
    pbt = []

    for file in tqdm(files, total=len(files)):
        # Read in each file as a dataframe
        df = pd.read_csv(file, sep="\t")
        df.columns = columns

        # Count by each region type and turn that into a percentage
        counts_by_type = df[["gene_id", "genic_region_type"]].groupby(by="genic_region_type").count()
        total_count = counts_by_type.sum().values.item()
        percent_by_type = (counts_by_type / total_count) * 100
        
        # Change the column to the file name without file extension for easier
        basename = os.path.basename(file).replace('.txt', '')
        percent_by_type.columns = [f"{basename}"]
        percent_by_type.fillna(value=0)
        pbt.append(percent_by_type)

    return pd.concat(pbt, axis=1, join='outer')

def clean_and_separate(output):
    spidr_cols = sorted([col for col in output.columns.tolist() if "_spidr_" in col])
    encode_cols = sorted([col for col in output.columns.tolist() if "_encode_" in col])

    spidr = output[spidr_cols]
    encode = output[encode_cols]

    new_spidr_cols = []
    for col in spidr.columns:
        if "Bethyl" in col or "CST" in col:
            new_spidr_cols.append("_".join(col.split("_")[0:2]))
        else:
            new_spidr_cols.append(col.split("_")[0])

    spidr.columns = new_spidr_cols

    new_encode_cols = []
    for col in encode.columns:
        if "Bethyl" in col or "CST" in col:
            new_encode_cols.append("_".join(col.split("_")[0:2]))
        else:
            new_encode_cols.append(col.split("_")[0])

    encode.columns = new_encode_cols

    common_cols = [col for col in new_encode_cols if col in new_spidr_cols]
    spidr = spidr[common_cols]
    encode = encode[common_cols]

    return spidr, encode

if __name__ == "__main__":
    files = get_files()
    output = get_percentages_by_region(files)
    spidr, encode = clean_and_separate(output)

    spidr.to_csv("processed-files/spidr_percent_by_region.csv")
    encode.to_csv("processed-files/encode_percent_by_region.csv")
    