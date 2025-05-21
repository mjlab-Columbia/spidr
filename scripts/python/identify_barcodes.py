from click import command, option, Path
import pandas as pd
import os
from typing import Tuple, Dict, List, Set
from tqdm import tqdm
import gzip
from collections import deque

# Number of bases between barcodes (equivalent to the length of odd overhang or even overhang)
# TODO: Move this into CLI options or config
SPACER = 7
LAXITY = 6
BARCODE_LENGTH = 17
NOT_FOUND = "NOT_FOUND"


def find_bead_id(chunk: Tuple[str],
                 read1_start_offset: int,
                 bead_id_lengths: List[int],
                 bead_hashmap: Dict[str, str]) -> str:
    """
    Find the bead id sequence in read 1 of paired-end read

    Args:
        chunk: Tuple[str] = tuple of 4 reads from a fastq (header, read, +, quality)
        read1_start_offset: int = number of bases from the 5' end of read 1 to start looking for bead id
        bead_id_lengths: List[int] = possible lengths of bead ids
        bead_hashmap: Dict[str, str] = dictionary of form { BEAD_SEQUENCE: BEAD_ID_NAME } to search for bead ids

    Returns:
        str = name of bead id found or NOT_FOUND
    """
    header, read, plus_line, quality = chunk

    # Iterate through all possible bead id lengths and check if the sequence is in the hashmap
    for length in bead_id_lengths:
        potential_match = read[read1_start_offset:(read1_start_offset + length)]
        if potential_match in bead_hashmap:
            return bead_hashmap[potential_match]

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
    # return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))
    return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))


def get_all_seqs_within_k(seq: str, k: int) -> Set[str]:
    """
    Get all sequences within a Hamming distance of k from a given sequence

    Args:
        seq: str = sequence to search around
        k: int = Hamming distance to search for

    Returns:
        Set[str] = set of all sequences within a Hamming distance of k from `seq`
    """
    alphabet = "ACGT"
    sequence_length = len(seq)
    possible_mutations = set()

    # Initialize the queue with the original sequence and distance 0
    queue = deque([(seq, 0)])

    while queue:
        current_seq, current_distance = queue.popleft()

        # If we reach the desired Hamming distance, add it to the result
        if current_distance == k:
            possible_mutations.add(current_seq)
            continue

        # Generate all possible sequences with Hamming distance + 1
        for position in range(sequence_length):
            for nucleotide in alphabet:
                if nucleotide != current_seq[position]:
                    mutated_seq = current_seq[:position] + \
                        nucleotide + current_seq[position + 1:]
                    queue.append((mutated_seq, current_distance + 1))

    return possible_mutations


def create_hamming_hashmap(barcode_df: pd.DataFrame) -> Dict[str, Set[str]]:
    """
    Create a hashmap of all barcodes with a Hamming distance of k from known barcodes

    Args:
        barcode_df: pd.DataFrame = dataframe containing barcode sequences
        k: int = Hamming distance to search for
    """
    cols = ["name", "sequence", "tolerance"]
    basic_hashmap = barcode_df[cols].to_records(index=False)
    hamming_hashmap = {}
    for name, seq, tol in basic_hashmap:
        hamming_hashmap[name] = get_all_seqs_within_k(seq, tol)
        hamming_hashmap[name].add(seq)

    return hamming_hashmap


def find_terminal_barcode(read: str,
                          barcode_hashmap: Dict[str, Set[str]],
                          start: int,
                          possible_lengths: Set[int]) -> Tuple[str, int]:
    """
    Identify the terminal barcode on read 2

    Args:
        read: str = read 2 from paired end reads
        barcode_hashmap: Dict[str, Set[str]] = hashmap with config entries for term barcodes and their Hamming neighbors
        start: int = 0-based index to start looking from in read

    Returns:
        str = terminal barcode name (e.g. "NYBot8") or NOT_FOUND
        int = new start position relative to 5' end of read 2
    """
    # For every possible terminal barcode length, check for terminal barcodes
    for term_size in possible_lengths:
        for terminal_bc in barcode_hashmap:
            if read[start:(start + term_size)] in barcode_hashmap[terminal_bc]:
                return terminal_bc, start + term_size

    # If we've gone through all possible sizes and haven't returned anything
    return NOT_FOUND, start + term_size


def find_nonterminal_barcode(read: str,
                             barcode_hashmap: Dict[str, Set[str]],
                             start: int,
                             possible_lengths: List[int] = [15]) -> Tuple[str, int]:
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
    # To account for barcodes of different lengths, we need to search for each possible barcode size
    possible_bc_sizes = sorted(possible_lengths)
    max_bc_size = max(possible_bc_sizes)

    # Barcodes can be offset within the range [0, LAXITY]
    for offset in range(0, LAXITY + 1):
        # Window start doesn't change with barcode length, but the window end does
        window_start = start + offset

        # If there's multiple sizes of barcodes, search by each barcode size
        for barcode_size in possible_bc_sizes:
            window_end = window_start + barcode_size
            substring = read[window_start:window_end]

            # If we have less than a barcode's worth of bases from the substring, we won't find another barcode
            if len(substring) < barcode_size:
                return NOT_FOUND, start + BARCODE_LENGTH

            # If we find a barcode in the hashmap, return the barcode and the new start position
            for bc in barcode_hashmap:
                if substring in barcode_hashmap[bc]:
                    return bc, start + BARCODE_LENGTH

            is_max_size = (barcode_size == max_bc_size)
            is_max_laxity = (offset == LAXITY)
            if (is_max_size) and (is_max_laxity):
                return NOT_FOUND, start + BARCODE_LENGTH


def find_barcodes(read: str,
                  even_hashmap: Dict[str, Set[str]],
                  odd_hashmap: Dict[str, Set[str]],
                  term_hashmap: Dict[str, Set[str]],
                  read2_format: List[str],
                  read2_start_offset: int,
                  possible_term_barcode_lengths: Set[int],
                  possible_odd_barcode_lengths: Set[int],
                  possible_even_barcode_lengths: Set[int]) -> List[str]:
    """
    Find all barcode sequences in read 2

    Args:
        read: str = read 2 of paired end reads from 5' to 3'
        odd_hashmap: Dict[str, Set[str]] = hashmap of all odd barcodes and their Hamming neighbors
        even_hashmap: Dict[str, Set[str]] = hashmap of all even barcodes and their Hamming neighbors
        term_hashmap: Dict[str, Set[str]] = hashmap of all terminal barcodes and their Hamming neighbors
        read2_format: List[str] = barcode categories listed from 5' to 3' on read 2 (e.g. ['Y', 'SPACER', 'ODD'])

    Returns:
        List[str] = list of barcodes from name column of config_df
    """
    start = read2_start_offset
    barcodes = []
    layout = read2_format.copy() + ["END"]
    barcode_type = layout.pop(0)

    while len(layout) > 0:
        if barcode_type == "Y":
            bc, new_start = find_terminal_barcode(read, term_hashmap, start, possible_term_barcode_lengths)
            barcodes.append(bc)
            start = new_start
        elif barcode_type == "SPACER":
            start += SPACER
        elif barcode_type == "END":
            break
        elif barcode_type == "ODD":
            bc, new_start = find_nonterminal_barcode(read, odd_hashmap, start, possible_odd_barcode_lengths)
            barcodes.append(bc)
            start = new_start
        elif barcode_type == "EVEN":
            bc, new_start = find_nonterminal_barcode(read, even_hashmap, start, possible_even_barcode_lengths)
            barcodes.append(bc)
            start = new_start
        else:
            raise Exception("Invalid barcode type")

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


@command()
@option('--input_read1', type=Path(exists=True), help='Unbarcoded read 1 fastq file')
@option('--input_read2', type=Path(exists=True), help='Unbarcoded read 2 fastq file')
@option('--output_read1', help='Path to output barcoded read 1 fastq file')
@option('--output_read2', help='Path to output barcoded read 2 fastq file')
@option('--read1_format', help="Read 1 barcode format string (e.g. 'DPM')")
@option('--read2_format', help="Read 2 barcode format (e.g. 'Y|SPACER|ODD|SPACER|EVEN')")
@option('--read1_start_offset', type=int,
        help="Skip this many bases from 5' end on read1 before bead oligo search", default=0, show_default=True)
@option('--read2_start_offset', type=int,
        help="Skip this many bases from 5' on read2 before terminal barcode search", default=0, show_default=True)
@option('--config', type=Path(exists=True), help='Config file contains bead sequences')
@option('--show_progress_bar', type=bool, help='Whether or not to show a tqdm progress bar', default=False)
def main(input_read1: os.PathLike,
         input_read2: os.PathLike,
         output_read1: os.PathLike,
         output_read2: os.PathLike,
         read1_format: str,
         read2_format: str,
         read1_start_offset: str,
         read2_start_offset: str,
         config: os.PathLike,
         show_progress_bar: bool) -> None:
    """
    Barcode identification script for SPIDR pipeline

    Last Updated: 2024-05-29
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
    bead_id_lengths = bead_df["sequence"].str.len().unique()

    # Keep a dataframe for each barcode type
    odd_df = config_df[config_df["type"] == "ODD"].copy()
    even_df = config_df[config_df["type"] == "EVEN"].copy()
    term_df = config_df[config_df["type"] == "Y"].copy()
    odd_hashmap = create_hamming_hashmap(barcode_df=odd_df)
    even_hashmap = create_hamming_hashmap(barcode_df=even_df)
    term_hashmap = create_hamming_hashmap(barcode_df=term_df)
    odd_barcode_lengths = set(odd_df.sequence.map(len))
    even_barcode_lengths = set(odd_df.sequence.map(len))
    term_barcode_lengths = set(term_df.sequence.map(len))

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
                    print(f"Processed {counter} reads")

            header_r1, read_r1, plus_line_r1, quality_r1 = chunk_read1
            header_r2, read_r2, plus_line_r2, quality_r2 = chunk_read2

            # Remove the space and everything after it in the header, since it's not needed
            # Keeping the space causes issues with bowtie2 and STAR
            header_r1 = header_r1.split(" ")[0]
            header_r2 = header_r2.split(" ")[0]

            # Bead ID is only found in read 1
            bead_id = find_bead_id(chunk_read1,
                                   read1_start_offset,
                                   bead_id_lengths,
                                   bead_hashmap)

            # Find all non-bead barcodes
            barcodes = find_barcodes(read=read_r2,
                                     odd_hashmap=odd_hashmap,
                                     even_hashmap=even_hashmap,
                                     term_hashmap=term_hashmap,
                                     read2_format=read2_format,
                                     read2_start_offset=read2_start_offset,
                                     possible_term_barcode_lengths=term_barcode_lengths,
                                     possible_odd_barcode_lengths=odd_barcode_lengths,
                                     possible_even_barcode_lengths=even_barcode_lengths)

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
