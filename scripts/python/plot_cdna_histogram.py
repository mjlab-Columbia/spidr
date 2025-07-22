import click
import gzip
import matplotlib.pyplot as plt


@click.command()
@click.option("--read1", "-1", type=click.Path(exists=True), help="Path to text file where each line is length of cDNA in pre-aligned fastq per line in read 1")
@click.option("--read2", "-2", type=click.Path(exists=True), help="Path to text file where each line is length of cDNA in pre-aligned fastq per line in read 1")
@click.option("--output", "-o", type=click.Path(), help="Path to output file. Extension will determine file pdf regardless of extension")
def main(read1, read2, output) -> None:

    if read1.endswith(".gz"):
        with gzip.open(read1, "rb") as read1_file:
            read1_lengths = [int(line) for line in read1_file.readlines()]
    else:
        with open(read1, "r") as read1_file:
            read1_lengths = [int(line) for line in read1_file.readlines()]

    if read2.endswith(".gz"):
        with gzip.open(read2, "rb") as read2_file:
            read2_lengths = [int(line) for line in read2_file.readlines()]
    else:
        with open(read2, "r") as read2_file:
            read2_lengths = [int(line) for line in read2_file.readlines()]

    fig, ax = plt.subplots(1, 2, sharex=True)
    ax[0].hist(read1_lengths)
    ax[1].hist(read2_lengths)

    fig.tight_layout()
    fig.savefig(output, format="pdf", dpi=300)

    return


if __name__ == "__main__":
    main()
