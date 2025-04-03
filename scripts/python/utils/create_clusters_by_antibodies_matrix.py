import gzip
from pandas import DataFrame, concat, Series
from click import command, option, Path
from typing import List
from os import PathLike
from tqdm import tqdm


def parse_lines(lines: List[str]) -> DataFrame:
    progress_bar = tqdm(lines, total=len(lines))
    rows = []  # List[Series]

    for line in progress_bar:
        all_items = line.split("\t")
        barcode, beads_str, counts_str = all_items
        beads = beads_str.split(",")
        counts = [count.strip() for count in counts_str.split(",")]
        row = {bead: count for bead, count in zip(beads, counts)}
        row["cluster"] = barcode
        rows.append(Series(row))

    print("Concatenating rows")
    concat_rows = concat(rows, join="outer", axis=1)

    print("Transposing, indexing, and filling N/A values")
    df = concat_rows.T.set_index("cluster").fillna(0)
    return df


@command()
@option("--input", "-i", type=Path(exists=True), help="Path to output of find_beads_per_barcode.py")
@option("--output", "-o", type=Path(), help="Where to write matrix file as a tsv")
def main(input: PathLike, output: PathLike) -> None:
    if str(input).endswith(".gz"):
        with gzip.open(input) as f:
            lines_bytes = f.readlines()
            lines = [line.decode("utf-8").strip() for line in lines_bytes]
            df = parse_lines(lines)
            df.to_csv(output, sep="\t")
            print(f"Output saved to: {output}")
    else:
        with open(input, "r") as f:
            lines = f.readlines()
            df = parse_lines(lines)
            df.to_csv(output, sep="\t")
            print(f"Output saved to: {output}")


if __name__ == "__main__":
    main()
