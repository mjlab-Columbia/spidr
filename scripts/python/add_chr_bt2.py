#!/usr/bin/python
import argparse
import pysam
import assembly

def main():
    args = parse_arguments()
    do(args)

def parse_arguments():
    parser = argparse.ArgumentParser(description =
            "This program adds chr to alignments.")
    parser.add_argument('-i', '--input', action = 'store', metavar = 'FILE',
                        help = 'Input BAM file')
    parser.add_argument('-o', '--output', action = 'store', metavar = 'FILE',
                        help = 'Output BAM file containing reads that pass')
    parser.add_argument('--assembly', metavar = "ASSEMBLY",
                        action = "store",
                        choices = ["mixed", "mm9", "mm10", "hg19", "hg38", "none"],
                        default = "none",
                        help = "The genome assembly. (default mm9)")
    return parser.parse_args()

def do(args):

    input_file = pysam.AlignmentFile(args.input, "rb")
    output_file = pysam.AlignmentFile(args.output, "wb", template = input_file)


    filtered_count = 0
    out_count = 0
    with pysam.AlignmentFile(args.input, "rb") as input_file:
        if args.assembly == 'none':
            output_file = pysam.AlignmentFile(args.output, "wb", template = input_file)
            for read in input_file.fetch(until_eof = True):
                if "NOT_FOUND" in read.query_name:
                    read.query_name = read.query_name.replace('NOT_FOUND', 'RPM')
                    output_file.write(read)
                    out_count += 1
                else:
                    filtered_count += 1
        else:
            output_file = pysam.AlignmentFile(args.output, "wb", template = input_file)
            for read in input_file.fetch():
                if "NOT_FOUND" in read.query_name:
                    read.query_name = read.query_name.replace('NOT_FOUND', 'RPM')
                    output_file.write(read)
                    out_count += 1
                else:
                    filtered_count += 1

        output_file.close()

    print('Filtered reads:', filtered_count)
    print('Written out reads:', out_count)

if __name__ == "__main__":
    main()
