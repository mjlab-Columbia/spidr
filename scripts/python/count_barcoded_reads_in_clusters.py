import click
import pandas as pd
from pdb import set_trace


@click.command()
@click.option("--bpm_rpm_counts", type=click.Path(exists=True), required=True, help="Table with number of BPM and RPM reads per cluster")
@click.option("--output", "-o", type=click.Path(), required=True)
def main(bpm_rpm_counts, output):
    """
    Count the total number of cDNA reads that were assigned to a cluster
    Note that the BPM/RPM count table doesn't contain the identity of the RPM reads (i.e. which antibody they came from).
    As a result, you cannot apply the min_oligos and proportion cutoffs on this table.
    """

    print(f"Reading {bpm_rpm_counts}")
    if bpm_rpm_counts.endswith(".tsv.gz"):
        df = pd.read_csv(bpm_rpm_counts, sep="\t", compression="gzip", index_col=0)
    elif bpm_rpm_counts.endswith(".tsv"):
        df = pd.read_csv(bpm_rpm_counts, sep="\t", index_col=0)
    else:
        raise NotImplementedError("Input file type not supported")

    total_cdna_reads_assigned_to_clusters = df["RPM"].sum().item()

    with open(output, "w") as output_file:
        output_file.write(str(total_cdna_reads_assigned_to_clusters))

    print("Done!")


if __name__ == "__main__":
    main()
