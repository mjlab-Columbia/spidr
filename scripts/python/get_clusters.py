import argparse

import cluster as c


def main():
    args = parse_arguments()
    clusters = c.get_clusters(args.input, args.num_tags)
    c.write_clusters_to_file(clusters, args.output, labeled=False)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Generates a clusters file from a BAM file.')
    parser.add_argument('-i', '--input',
                        metavar="FILE",
                        action="store",
                        nargs='+',
                        required=True,
                        help="The input BAM file(s).")
    parser.add_argument('-o', '--output',
                        metavar="FILE",
                        action="store",
                        required=True,
                        help="The output clusters file.")
    parser.add_argument('-n', '--num_tags',
                        metavar='INT',
                        type=int,
                        action='store',
                        required=True,
                        help="The number of tags contained in the barcode " +
                        "of each BAM record.")
    return parser.parse_args()


if __name__ == "__main__":
    main()
