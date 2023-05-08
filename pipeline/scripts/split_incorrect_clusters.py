import pandas as pd
import numpy as np
import re
import argparse
from collections import defaultdict, Counter
import tqdm
import gzip
import os


def main():
    args = parse_args()
    formatdict = load_format(args.format)
    read_path = args.clusters

    complete_out_path = read_path.replace('.clusters', '.complete.clusters')
    incomplete_out_path =  read_path.replace('.clusters', '.incomplete.clusters')
    counter = 0
    complete = 0
    incomplete = 0
    print(complete_out_path)
    print(incomplete_out_path)
    
    with open(read_path, 'r') as clusters, \
    open(complete_out_path, 'wt') as complete_out, \
    open(incomplete_out_path, 'wt') as incomplete_out:
        for line in clusters:
                counter +=1
                cluster_barcode = line.strip('\n').split('\t', 1)[0]
                barcodes = cluster_barcode.split('.')[:-1]
                if counter % 100000 == 0:
                    print(counter)
                tags = np.array(barcodes)
                indexed = [i for i, t in enumerate(tags) if i!=formatdict[t]]
                if len(indexed) != 0:
                    incomplete +=1
                    incomplete_out.write(line)
                else:
                    complete += 1
                    complete_out.write(line)

    print('Total clusters: ', counter)
    print('Clusters with incorrect barcodes: ', incomplete)
    print('Clusters with correct barcodes:', complete)


def parse_args():

    parser = argparse.ArgumentParser(description='Split clusters into possible and impossible/NOT_FOUND barcodes')
    parser.add_argument('--clusters',
                        dest='clusters',
                        type=str,
                        required=True,
                        help='Input clusters file')
    parser.add_argument('-f', '--format',
                        dest='format',
                        type=str,
                        required=True,
                        help='Allowed barcodes at each position')
    return parser.parse_args()


def load_format(formatfile):
    df = pd.read_csv(formatfile, header=None, sep='\t')
    return df.set_index(1)[0].to_dict(into=defaultdict(lambda: -1))


if __name__ == '__main__':
    main()
