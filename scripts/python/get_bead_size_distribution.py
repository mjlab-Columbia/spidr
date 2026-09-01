import os
import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm
from collections import defaultdict
from matplotlib import colormaps


def main():
    args = parse_arguments()
    cluster_counts = []
    read_counts = []
    for f in tqdm(args.files, desc="Processing cluster files", unit="file"):
        df1, df2 = count_statistics(f, args.readtype)
        cluster_counts.append(df1)
        read_counts.append(df2)
    cluster_df = pd.concat(cluster_counts, axis=1).transpose()
    read_df = pd.concat(read_counts, axis=1).transpose()
    cluster_fig = plot_profile(cluster_df)
    cluster_fig.savefig(
        os.path.join(args.output_dir, args.prefix + args.readtype + '_cluster_distribution.pdf'),
        bbox_inches='tight',
    )
    read_fig = plot_profile(read_df)
    read_fig.savefig(
        os.path.join(args.output_dir, args.prefix + args.readtype + '_read_distribution.pdf'),
        bbox_inches='tight',
    )


def count_statistics(filename, readtype):
    bins = np.array([0, 1, 5, 10, 20, 50, 100, 200])
    cluster_counts = defaultdict(int)
    read_counts = defaultdict(int)
    for bin in np.arange(1, len(bins) + 1):
        cluster_counts[bin] = 0
        read_counts[bin] = 0
    with open(filename, 'r') as f:
        line_count = sum(1 for _ in f)
    with open(filename, 'r') as f:
        for line in tqdm(f, total=line_count, desc=os.path.basename(filename), unit="cluster"):
            if readtype in line:
                reads = line.split('\t', 1)[1]
                counts = reads.count(readtype)
                bin = np.digitize(counts, bins, right=True)
                cluster_counts[bin] += 1
                read_counts[bin] += counts
    df_cluster_counts = pd.DataFrame.from_dict(cluster_counts, orient='index')
    df_cluster_counts.columns = [filename.rsplit('/', 1)[-1].rsplit('.clusters')[0]]
    df_read_counts = pd.DataFrame.from_dict(read_counts, orient='index')
    df_read_counts.columns = [filename.rsplit('/', 1)[-1].rsplit('.clusters')[0]]
    return df_cluster_counts, df_read_counts


def plot_profile(df):
    columns = ['1', '2-5', '6-10', '11-20', '21-50', '51-100', '101-200', '201+']
    df.columns = columns
    df = df.div(df.sum(axis=1), axis=0)
    plot = df.plot(kind='bar', stacked=True, ylabel='Proportion', cmap=colormaps.get_cmap('Dark2'))
    plot.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plot.set_ylabel = 'Proportion'
    return plot.get_figure()


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Generate the cluster size distribution plot.')

    parser.add_argument('--files',
                        nargs='+',
                        metavar="FILE",
                        required=True,
                        help="Cluster files for a single experiment")
    parser.add_argument('--output_dir',
                        metavar="DIR",
                        required=True,
                        help="Directory to write distribution plots")
    parser.add_argument('--prefix',
                        metavar="PREFIX",
                        default="",
                        help="Filename prefix for output plots (e.g. experiment name)")
    parser.add_argument('--readtype',
                        required=True,
                        help="The read type to count")
    return parser.parse_args()


if __name__ == "__main__":
    main()
