import pandas as pd
import numpy as np
import argparse
from collections import defaultdict
import seaborn as sns

def parse_args():

    parser = argparse.ArgumentParser(description='Profile the distribution of clusters and bead oligos for ChIP-n-DIP experiments')
    parser.add_argument('-c', '--clusters', dest='clusters', type=str, required=True,
                        help='Cluster file')
    parser.add_argument('-o', '--output', dest='output', type=str, required=True, help='Directory to save graphs to')
    opts = parser.parse_args()

    return opts

def generate_statistics(cluster_file):
    overall_stats =  list()
    bead_stats = list()
    chromatin_stats = list()
    with open(cluster_file, 'r') as f:
        for line in tqdm.tqdm(f):
            barcode, *reads = line.rstrip('\n').split('\t')
            beads = [read for read in reads if 'BPM' in read]
            chromatin_counts = len(reads) - len(beads)
            beads_counts = len(beads)
            overall_stats.append(len(reads))
            bead_stats.append(bead_counts)
            chromatin_stats.append(chromatin_counts)
    bins = [0,1,2,5,10,20,50,100,1000,100000]    
    results = {}
    results['Bead'] = np.histogram(bead_stats, bins)[0]/len(bead_stats)
    results['Chromatin'] = np.histogram(chromatin_stats, bins)[0]/len(chromatin_stats)
    results['Overall'] = np.histogram(overall_stats, bins)[0]/len(overall_stats)
    df = pd.DataFrame.from_dict(results).T
    df.columns = ['None', 'Singleton', '2-5', '6-10', '11-20', '21-50','51-100', '100-1000', '1000+']
    ax = df.plot.bar(stacked=True)
    ax.set(ylabel = 'Proportion of Clusters')
    ax.legend(loc='center left', bbox_to_anchor = [1.0,0.5])
    return df, ax.get_figure()


def max_representation_ecdf(clusterfile, ax):
    results = []
    with open(clusterfile, 'r') as clusters:
        for line in tqdm.tqdm(clusters):
            barcode, *reads = line.rstrip('\n').split('\t')
            bead_reads = [read for read in reads if read.startswith('BPM')]
            if len(bead_reads) == 0:
                results.append(0)
            else:
                bead_labels = Counter([read.split(':')[0].split('_',1)[1] for read in bead_reads])
                candidate = bead_labels.most_common()[0]
                results.append(candidate[1]/len(bead_reads))
    if ax == 'None':
        ax = plt.figure().subplots()
    ax = sns.ecdfplot(results, linewidth = 3, ax=ax)
    ax.set(xlabel = 'Max Bead Representation Proportion', ylabel='Proportion of Clusters')
    fig = ax.get_figure()
    return results, fig


def main():
    args = parse_arguments()
    out_max_rep = args.out + '.bead_ecdf.pdf'
    _ , fig_max_rep = maximum_represetnation_ecdf(args.clusters, 'None')
    out_cluster_zies = args.out + '.cluster_sizes.pdf'
    

for f in files:
       df1, df2 = count_statistics(f) 
       results1.append(df1)
       results2.append(df2)
    df1 = pd.concat(results1, axis=1).transpose()
    df2 = pd.concat(results2, axis=1).transpose()
    df1 = prepare_df(df1)
    df2 = prepare_df(df2)
    plot_profile(df1, df2, args.directory)


def count_statistics(filename):
    counts = []
    with open(filename, 'r') as f:
        for line in f:
            reads = line.split('\t')[1:]
            counts.append(len(reads))
    bins = pd.IntervalIndex.from_tuples([(0, 1), (1, 10), (10, 100), (100,1000), (1000,1000000)])
    name = filename.split('/')[-1]
    df = pd.DataFrame(counts,columns= [name])
    count_df = df[name].value_counts(bins =  bins).to_frame()
    df['Groups'] = pd.cut(df[name],bins)
    norm_df = df.groupby('Groups').sum()
    return count_df, norm_df

def plot_profile(df1, df2, output):
    plot = df1.plot(kind='bar',stacked=True)
    plot.legend(loc='center left',bbox_to_anchor=(1.0, 0.5))
    plot.set_ylabel = 'Cluster counts'
    fig = plot.get_figure()
    fig.savefig(output + '/cluster_sizes_counts.png', bbox_inches='tight')
    norm = df2.div(df2.sum(axis=1), axis=0)
    plot2 = norm.plot(kind='bar', stacked=True)
    plot2.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plot2.set_ylabel = 'Proportion'
    fig2 = plot2.get_figure()
    fig2.savefig(output + '/cluster_sizes.png', bbox_inches = 'tight')

def prepare_df(df):
    df.columns = df.columns.astype(str)
    rename_mapping = {"(0, 1]": "Single", "(1, 10]":"1-10", "(10, 100]":"10-100", "(100, 1000]":"100-1000", "(1000, 1000000]": ">1000"}
    df = df.rename(columns=rename_mapping)
    df = df[['Single', '1-10', '10-100', '100-1000', '>1000']]
    return df 


if __name__ == "__main__":
    main()
            












