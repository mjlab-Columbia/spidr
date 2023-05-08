import pandas as pd
from collections import defaultdict, Counter
import seaborn as sns
import matplotlib.pyplot as plt
import cycler
import matplotlib as mpl
import numpy as np
import tqdm
import argparse
color = plt.cm.Dark2(np.linspace(0, 1, 8))
mpl.rcParams['axes.prop_cycle'] = cycler.cycler('color', color)

def main():
    args = parse_arguments()
    #ecdf_plot = max_representation_ecdf(args.clusters)
    #ecdf_plot.savefig(args.outprefix + '.max_representation_ecdf.png', bbox_inches='tight')
    ecdf_counts = max_representation_ecdf_counts(args.clusters, args.xlim)
    ecdf_counts.savefig(args.outprefix + '.max_representation_counts.png', bbox_inches='tight')


def max_representation_ecdf(clusterfile):
    results = []
    with open(clusterfile, 'r') as clusters:
        for line in tqdm.tqdm(clusters):
            barcode, *reads = line.rstrip('\n').split('\t')
            bead_reads = [read for read in reads if read.startswith('BPM')]
            if len(bead_reads) > 0:
                bead_labels = Counter([read.split(':')[0].split('_',1)[1] for read in bead_reads])
                candidate = bead_labels.most_common()[0]
                results.append(candidate[1]/len(bead_reads))
    ax = plt.figure().subplots()
    ax = sns.ecdfplot(results, linewidth = 3, ax=ax)
    ax.set(xlabel = 'Max Bead Representation Proportion', ylabel='Proportion of Clusters')
    return ax.get_figure()

def max_representation_ecdf_counts(clusterfile, xlimit):
    results = []
    with open(clusterfile, 'r') as clusters:
        for line in tqdm.tqdm(clusters):
            barcode, *reads = line.rstrip('\n').split('\t')
            bead_reads = [read for read in reads if read.startswith('BPM')]
            if len(bead_reads) > 0:
                bead_labels = Counter([read.split(':')[0].split('_',1)[1] for read in bead_reads])
                candidate = bead_labels.most_common()[0]
                results.append(candidate[1])
    ax = plt.figure().subplots()
    ax = sns.ecdfplot(results, linewidth = 3, ax=ax)
    ax.set(xlabel = 'Num Oligos for Most Common Type', ylabel='Proportion of Clusters')
    ax.set(xlim=(0, int(xlimit)))
    return ax.get_figure()


def parse_arguments():
    parser = argparse.ArgumentParser(
        description = 'Generate the maximum ecdf distribution plot.')
    parser.add_argument('--clusters',
                        metavar = "FILE",
                        action = "store",
                        help = "The clusters file")
    parser.add_argument('--outprefix',
                        action="store",
                        help = "The prefix to save plots with")
    parser.add_argument('--xlim',
                        action="store",
                        help = "The maximum x value on the counts plot")
    return parser.parse_args()

if __name__ == "__main__":
    main()



