import argparse
import cluster as c

def main():
    args = parse_arguments()
    c.merge_clusters(args.input, args.output)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description = 'Generates a clusters file from a BAM file.')
    parser.add_argument('-i', '--input',
                        metavar = "FILE",
                        action = "store",
                        required=True,
                        help = "The input clusters file.")
    parser.add_argument('-o', '--output',
                        metavar = "FILE",
                        action = "store",
                        required=True,
                        help = "The output clusters file.")
    return parser.parse_args()

if __name__ == "__main__":
    main()
