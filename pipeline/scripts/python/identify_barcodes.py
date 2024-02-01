import click
import pandas as pd
from pdb import set_trace as st
import os
from typing import Tuple, Dict, List
from tqdm import tqdm
import gzip
import numpy as np

# Number of bases between barcodes (equivalent to the length of odd overhang or even overhang)
# TODO: Move this into CLI options or config
SPACER = 7
LAXITY = 6
BARCODE_LENGTH = 17
NOT_FOUND = "NOT_FOUND"


def find_bead_id(chunk: Tuple[str],
                 start_offset: int,
                 sequence_length: int,
                 bead_hashmap: Dict[str, str]) -> str:
    """
    Find the bead id sequence in read 1 of paired-end read

    Args:
        chunk: Tuple[str] = tuple of 4 reads from a fastq (header, read, +, quality)
        start_offset: int = number of bases from the 5' end of read 1 to start looking for bead id
        sequence_length: int = length of bead id to search for (this might get deprecated)
        bead_hashmap: Dict[str, str] = dictionary of form { BEAD_SEQUENCE: BEAD_ID_NAME } to efficiently search for bead ids

    Returns:
        str = name of bead id found or NOT_FOUND
    """
    header, read, plus_line, quality = chunk
    potential_match = read[start_offset:(start_offset + sequence_length)]

    if potential_match in bead_hashmap:
        return bead_hashmap[potential_match]
    else:
        return NOT_FOUND


def hamming_distance(seq1: str, seq2: str) -> int:
    """
    Calculate Hamming distance between 2 string of the same length.
    Equivalent to the number of mismatches in 2 strings of the same length.

    Args:
        seq1: str = first sequence to compare
        seq2: str = second sequence to compare

    Returns:
        int: Hamming distance between strings
    """
    assert len(seq1) == len(seq2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))


def find_terminal_barcode(read: str,
                          barcode_df: pd.DataFrame,
                          start: int) -> Tuple[str, int]:
    """
    Identify the terminal barcode on read 2

    Args:
        read: str = read 2 from paired end reads
        barcode_df: pd.DataFrame = dataframe containing config entries for terminal barcodes
        start: int = 0-based index to start looking from in read

    Returns:
        str = terminal barcode name (e.g. "NYBot8") or NOT_FOUND
        int = new start position relative to 5' end of read 2
    """
    # For every possible terminal barcode length, check for terminal barcodes
    possible_term_sizes = barcode_df["sequence"].str.len().unique()
    possible_term_sizes = sorted(possible_term_sizes.tolist())

    for term_size in possible_term_sizes:
        cols = ["name", "sequence", "tolerance"]
        size_mask = barcode_df["sequence"].str.len() == term_size
        term_df = barcode_df[size_mask]
        seq_tol = term_df[cols].to_records(index=False)

        # Search for terminal barcodes within a Hamming distance of 2 from sequence
        for name, seq, tol in seq_tol:
            dist = hamming_distance(seq, read[:term_size])
            if dist <= tol:
                return name, start + term_size

    # If we've gone through all possible sizes and haven't returned anything
    return NOT_FOUND, start + term_size


def find_nonterminal_barcode(read: str,
                             barcode_df: pd.DataFrame,
                             start: int) -> Tuple[str, int]:
    """
    Find non-terminal barcodes in read 2

    Args:
        read: str = read 2 from paired end reads
        barcode_df: pd.DataFrame = dataframe containing config entries for specific barcode type (e.g. ODD)
        start: int = 0-based index to start looking from in read

    Returns:
        str = barcode name (e.g. "ROUND2_B2")
        int = new start position relative to 5' end of read 2
    """
    # Barcodes can be offset within the range [0, LAXITY]
    for offset in range(0, LAXITY + 1):
        possible_bc_sizes = barcode_df["sequence"].str.len().unique()
        possible_bc_sizes = sorted(possible_bc_sizes.tolist())

        # If there's multiple sizes of barcodes, search by each barcode size
        for barcode_size in possible_bc_sizes:
            window_start = start + offset
            window_end = window_start + barcode_size
            substring = read[window_start:window_end]

            # If we have less than a barcode's worth of bases from the substring, we won't find another barcode
            if len(substring) < barcode_size:
                return NOT_FOUND, start + BARCODE_LENGTH

            cols = ["name", "sequence", "tolerance"]
            mask = barcode_df["sequence"].str.len() == barcode_size
            bc_df = barcode_df[mask]
            seq_tol = bc_df[cols].to_records(index=False)

            # Search for barcodes within a Hamming distance of 2
            for name, seq, tol in seq_tol:
                dist = hamming_distance(seq, substring)
                if dist <= tol:
                    return name, start + BARCODE_LENGTH

            is_max_size = (barcode_size == max(possible_bc_sizes))
            is_max_laxity = (offset == LAXITY)
            if (is_max_size) and (is_max_laxity):
                return NOT_FOUND, start + BARCODE_LENGTH


def find_barcodes(read: str,
                  config_df: pd.DataFrame,
                  read2_format: List[str]) -> List[str]:
    """
    Find all barcode sequences in read 2

    Args:
        read: str = read 2 of paired end reads from 5' to 3' 
        config_df: pd.DataFrame = configuration file with sequence type, name, sequence bases, and tolerance
        read2_format: List[str] = barcode categories listed from 5' to 3' on read 2 (e.g. ['Y', 'SPACER', 'ODD'])

    Returns:
        List[str] = list of barcodes from name column of config_df
    """
    start = 0
    barcodes = []
    layout = read2_format.copy() + ["END"]
    barcode_type = layout.pop(0)

    while len(layout) > 0:
        type_mask = config_df["type"] == barcode_type
        barcode_df = config_df[type_mask]

        # Y is the label for a terminal barcode
        if barcode_type == "Y":
            bc, new_start = find_terminal_barcode(read, barcode_df, start)
            barcodes.append(bc)
            start = new_start
        elif barcode_type == "SPACER":
            start += SPACER
        elif barcode_type == "END":
            break
        else:
            bc, new_start = find_nonterminal_barcode(read, barcode_df, start)
            barcodes.append(bc)
            start = new_start

        # If there isn't a BAROCDE_LENGTH's worth of bases, it's impossible to find another barcode
        if (len(read) - start) < BARCODE_LENGTH:
            break

        # pop the first element from layout (e.g. ['A', 'B'] - -> 'A')
        barcode_type = layout.pop(0)

    return barcodes


def pad_barcodes(barcodes: List[str], expected_length: int) -> List[str]:
    """
    Pad barcodes to expected length

    Args:
        barcodes: List[str] = list of barcodes
        expected_length: int = number of barcodes expected in `barcodes`

    Returns:
        List[str] = list of barcodes padded up to `expected_length`
    """
    if len(barcodes) < expected_length:
        return barcodes + [NOT_FOUND] * (expected_length - len(barcodes))
    elif len(barcodes) == expected_length:
        return barcodes
    else:
        raise Exception("barcodes is greater than expected_length")


@ click.command()
@ click.option('--input_read1', type=click.Path(exists=True), help='Unbarcoded read 1 fastq file')
@ click.option('--input_read2', type=click.Path(exists=True), help='Unbarcoded read 2 fastq file')
@ click.option('--output_read1', help='Path to output barcoded read 1 fastq file')
@ click.option('--output_read2', help='Path to output barcoded read 2 fastq file')
@ click.option('--read1_format', help="Read 1 barcode format string (e.g. 'DPM')")
@ click.option('--read2_format', help="Read 2 barcode format (e.g. 'Y|SPACER|ODD|SPACER|EVEN')")
@ click.option('--start_offset', type=int, help="Bases to skip from 5' end before bead oligo search")
@ click.option('--sequence_length', type=int, help="Length of sequences to identify")
@ click.option('--config', type=click.Path(exists=True), help='Config file contains bead sequences')
@ click.option('--show_progress_bar', type=bool, help='Whether or not to show a tqdm progress bar', default=False)
def main(input_read1: os.PathLike, input_read2: os.PathLike, output_read1: os.PathLike, output_read2: os.PathLike, read1_format: str, read2_format: str, start_offset: str, sequence_length: int, config: os.PathLike) -> None:
    """
    Entry point for the program
    """
    # Turn format strings into lists of strings (e.g. "A|B" --> ["A", "B"])
    read1_format = [s.strip() for s in read1_format.split('|')]
    read2_format = [s.strip() for s in read2_format.split('|')]
    num_read1_barcodes = len([b for b in read1_format if b != "SPACER"])
    num_read2_barcodes = len([b for b in read2_format if b != "SPACER"])
    total_barcode_length = num_read1_barcodes + num_read2_barcodes

    # Load configuration dataframe
    config_df = pd.read_csv(config,
                            sep="\t",
                            skiprows=3,
                            usecols=[0, 1, 2, 3],
                            names=["type", "name", "sequence", "tolerance"])

    # Create a dictionary of the format { BEAD_SEQUENCE: BEAD_NAME } for all DPM entries
    bead_df = config_df[config_df["type"] == "DPM"][["name", "sequence"]]
    bead_hashmap = bead_df.set_index("sequence").to_dict()["name"]

    # Read from inputs and write to outputs in same order to preserve fastq header order
    with gzip.open(input_read1, "rb") as read1_in, \
            gzip.open(output_read1, "wb") as read1_out, \
            gzip.open(input_read2, "rb") as read2_in, \
            gzip.open(output_read2, "wb") as read2_out:

        # fastq files have 4 lines per read
        chunk_size = 4

        # Extract all lines from fastq
        read1_lines = [line.decode("utf-8").strip() for line in read1_in]
        read2_lines = [line.decode("utf-8").strip() for line in read2_in]
        num_chunks = len(read1_lines) / chunk_size

        # Iterators that return tuples of 4 lines at a time (header, read, +, quality)
        read1_iterator = zip(*[iter(read1_lines)] * chunk_size)
        read2_iterator = zip(*[iter(read2_lines)] * chunk_size)
        progress_bar = tqdm(zip(read1_iterator, read2_iterator),
                            total=num_chunks)

        if show_progress_bar:
            progress_bar = tqdm(zip(read1_iterator, read2_iterator),
                                total=num_chunks)
            counter = None
        else:
            progress_bar = zip(read1_iterator, read2_iterator)
            counter = 0

        # Iterate through read 1 and 2 in chunks of 4 (i.e. one read + metadata at a time)
        for chunk_read1, chunk_read2 in progress_bar:

            # Display progress if we're using a counter
            if counter is not None:
                counter += 1
                if counter % 10000 == 0:
                    progress_bar.set_description(f"Processed {counter} reads")

            header_r1, read_r1, plus_line_r1, quality_r1 = chunk_read1
            header_r2, read_r2, plus_line_r2, quality_r2 = chunk_read2

            # Bead ID is only found in read 1
            bead_id = find_bead_id(chunk_read1,
                                   start_offset,
                                   sequence_length,
                                   bead_hashmap)

            # Find all non-bead barcodes
            barcodes = find_barcodes(read_r2,
                                     config_df,
                                     read2_format)

            # Concatenate bead id and barcodes and pad with ["NOT_FOUND"] if necessary
            final_barcodes = [bead_id] + barcodes
            padded_barcodes = pad_barcodes(final_barcodes,
                                           total_barcode_length)
            barcode_string = "".join([f"[{bc}]" for bc in padded_barcodes])
            modified_header_r1 = header_r1 + "::" + barcode_string
            modified_header_r2 = header_r2 + "::" + barcode_string

            read_r1_out = "\n".join(
                [modified_header_r1, read_r1, plus_line_r1, quality_r1])
            read_r2_out = "\n".join(
                [modified_header_r2, read_r2, plus_line_r2, quality_r2])

            # Write output as bytes
            chunk_out_r1 = read_r1_out + "\n"
            chunk_out_r2 = read_r2_out + "\n"
            read1_out.write(chunk_out_r1.encode("utf-8"))
            read2_out.write(chunk_out_r2.encode("utf-8"))


if __name__ == '__main__':
    main()
