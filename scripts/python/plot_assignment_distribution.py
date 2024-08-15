import argparse
from collections import defaultdict, Counter
import re
import tqdm
import pandas as pd

def parse_args():

    parser = argparse.ArgumentParser(description='Generate pie charts for read and cluster assignment distributions')
    parser.add_argument('-c', '--clusters',
                        dest='clusters',
                        type=str,
                        required=True,
                        help='Clusters from which to assign protein tag')
    parser.add_argument('--min_oligos',
                        action = "store",
                        type = int,
                        help = 'The minimum number of oligos to call a cluster')
    parser.add_argument('--proportion',
                        action = "store",
                        type = float,
                        help = 'The maximum representation proportion')
    parser.add_argument('--max_size',
                        action = "store",
                        type = int,
                        help = "The maximum cluster size to keep")
    parser.add_argument('--outprefix',
                        action="store",
                        help = 'Prefix used for saving plots')
    opts = parser.parse_args()

    return opts


def main():
    args = parse_args()
    labels,dpm_counts, total_count = assign_labels(args.clusters,args.min_oligos, args.proportion, args.max_size)
    df = build_count_df(labels, dpm_counts)
    print('Filtered = Clusters are too big')
    print('Uncertain = Too few oligos')
    print('Ambiguous = Too low proportion')
    print('None = No co-clustering oligo')
    print('Total number of clusters: ', total_count)
    print(df)
    df.to_csv(args.outprefix + '.counts', index=True, header=True, sep='\t')
    metadata = 'Min Oligos:' + str(args.min_oligos) + ', Proportion:' + str(args.proportion) + ', Max size: ' + str(args.max_size)
    plot_pies(df, metadata, args.outprefix)

def plot_pies(df, metadata, outprefix):
    df = df.drop(['ambiguous', 'uncertain','none'])
    ax_cluster = df.plot(kind='pie', title = metadata, y='Cluster Counts', labels = df.index, legend=False)
    ax_cluster.yaxis.set_visible(False)
    ax_cluster.get_figure().savefig(outprefix + '.clusters_pie.png', bbox_inches='tight')
    ax_reads = df.plot(kind='pie', title = metadata, y='Read Counts', labels = df.index, legend=False)
    ax_reads.yaxis.set_visible(False)
    ax_reads.get_figure().savefig(outprefix + '.reads_pie.png', bbox_inches='tight')
    ax_ratio = df.plot(kind='pie', title = metadata, y = 'Reads/Cluster', labels=df.index, legend=False)
    ax_ratio.yaxis.set_visible(False)
    ax_ratio.get_figure().savefig(outprefix + '.readsPerCluster_pie.png', bbox_inches='tight')
    ax_bar = df.plot(kind='bar', title = metadata, y = 'Reads/Cluster', figsize = (12,6), ylabel = 'Reads/Cluster', legend=False)
    ax_bar.get_figure().savefig(outprefix + '.readsPerCluster_bar.png', bbox_inches='tight')


def build_count_df(labels, dpm_counts):
    df = pd.DataFrame.from_dict(labels, orient='index')
    df_counts = pd.DataFrame.from_dict(dpm_counts, orient='index')
    df = df.join(df_counts, lsuffix='cluster', rsuffix='reads')
    df.columns = ['Cluster Counts', 'Read Counts']
    df['Reads/Cluster'] = df['Read Counts']/df['Cluster Counts']
    df = df.sort_values(by='Reads/Cluster', ascending=False)
    return df

def assign_labels(clusterfile, min_oligos, proportion, max_size):
    labels = defaultdict(int)
    dpm = defaultdict(int)
    counter = 0
    with open(clusterfile, 'r') as clusters:
        for line in tqdm.tqdm(clusters):
            counter +=1 
            barcode, *reads = line.rstrip('\n').split('\t')
            assignment = get_single_label(reads, min_oligos, proportion, max_size)
            labels[assignment] += 1
            dpm_counts = line.count('DPM')
            dpm[assignment] += dpm_counts
    return labels, dpm, counter

def get_single_label(reads,min_oligos, proportion, max_size):
     """
     Identify the max represented bead label for a set of reads 
     Inputs:
         max_size(int) = maximum number of chromatin in cluster
     Returns:
         label(str) = name of protein associated with cluster
     """
     bead_reads = [read for read in reads if read.startswith('BPM')]
     if len(bead_reads) == 0:
        return 'none'
     cluster_size = len(reads) - len(bead_reads)
     if int(cluster_size) > int(max_size):
        return 'filtered'
     bead_labels = Counter([read.split(':')[0].split('_',1)[1] for read in bead_reads])
     candidate = bead_labels.most_common()[0]
     if candidate[1] < min_oligos:
         return 'uncertain'
     elif candidate[1]/sum(bead_labels.values()) < proportion:
         return 'ambiguous'
     else:
         return candidate[0]
     return 'malformed'


if __name__== "__main__":
    main()

