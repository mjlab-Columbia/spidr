import click
import matplotlib.pyplot as plt
import numpy as np


@click.command()
@click.option("--read1", "-1", type=click.Path(exists=True),
              help="Path to text file where each line is length of cDNA in pre-aligned fastq per line in read 1")
@click.option("--read2", "-2", type=click.Path(exists=True),
              help="Path to text file where each line is length of cDNA in pre-aligned fastq per line in read 1")
@click.option("--output", "-o", type=click.Path(),
              help="Path to output file. Extension will determine file pdf regardless of extension")
def main(read1, read2, output) -> None:

    read1_lengths = np.loadtxt(read1, dtype=np.uint8)
    read2_lengths = np.loadtxt(read2, dtype=np.uint8)

    fig, ax = plt.subplots(1, 2, figsize=(12, 5), sharex=False)
    ax[0].hist(read1_lengths, bins=50)
    ax[1].hist(read2_lengths, bins=50)

    ax[0].set_title("Read 1 cDNA Lengths")
    ax[0].set_ylabel("Frequency")
    ax[0].set_xlabel("cDNA Length")

    ax[1].set_title("Read 2 cDNA Lengths")
    ax[1].set_ylabel("Frequency")
    ax[1].set_xlabel("cDNA Length")

    fig.tight_layout()
    fig.savefig(output, format="pdf", dpi=300)

    return


if __name__ == "__main__":
    main()
