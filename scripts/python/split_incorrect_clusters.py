import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm
from collections import defaultdict


def main():
    args = parse_args()
    formatdict = load_format(args.format)
    read_path = args.clusters
    complete_out_path = args.complete_output
    incomplete_out_path = args.incomplete_output

    complete = 0
    incomplete = 0
    print(complete_out_path)
    print(incomplete_out_path)

    line_count = 0
    with open(read_path, 'r') as clusters:
        for line in clusters:
            line_count += 1

    # TODO: Vectorize this code with pandas
    with open(read_path, 'r') as clusters, \
            open(complete_out_path, 'wt') as complete_out, \
            open(incomplete_out_path, 'wt') as incomplete_out:
        for line in tqdm(clusters, total=line_count):
            cluster_barcode = line.strip('\n').split('\t', 1)[0]
            barcodes = cluster_barcode.split('.')[:-1]
            tags = np.array(barcodes)
            indexed = [i for i, t in enumerate(tags) if not i in formatdict[t]]

            if len(indexed) != 0:
                incomplete += 1
                incomplete_out.write(line)
            else:
                complete += 1
                complete_out.write(line)

    print('Total clusters: ', complete + incomplete)
    print('Clusters with incorrect barcodes: ', incomplete)
    print('Clusters with correct barcodes:', complete)


def parse_args():

    parser = argparse.ArgumentParser(description='Split clusters into possible and impossible/NOT_FOUND barcodes')
    parser.add_argument('--clusters',
                        dest='clusters',
                        type=str,
                        required=True,
                        help='Input clusters file')
    parser.add_argument('--complete_output',
                        dest='complete_output',
                        type=str,
                        required=True,
                        help='Path of output cluster file with complete clusters')
    parser.add_argument('--incomplete_output',
                        dest='incomplete_output',
                        type=str,
                        required=True,
                        help='Path of output cluster file with complete clusters')
    parser.add_argument('-f', '--format',
                        dest='format',
                        type=str,
                        required=True,
                        help='Allowed barcodes at each position')
    return parser.parse_args()


def load_format(formatfile):
    df = pd.read_csv(formatfile, header=None, sep='\t')
    df = pd.DataFrame(df.set_index(1)[0])

    result_dict = defaultdict(lambda: [-1])

    for index, value in df.iterrows():
        result_dict[index].append(int(value[0]))

    return result_dict


if __name__ == '__main__':
    main()
