import argparse
import pysam
import re
import os
from collections import defaultdict, Counter
from pathlib import Path
import itertools
import tqdm

def parse_args():

    parser = argparse.ArgumentParser(description='Add protein label to DNA bamfile')
    parser.add_argument('-i', '--input_bam',
                        dest='input_bam',
                        type=str,
                        required=True,
                        help='Aligned DNA Bamfile')
    parser.add_argument('-o', '--output_bam',
                        dest='output_bam',
                        type=str,
                        required=True,
                        help='Path to output bam with protein tag added')
    parser.add_argument('-c', '--clusters',
                        dest='clusters',
                        type=str,
                        required=True,
                        help='Clusters from which to assign protein tag')
    parser.add_argument('--num_tags',
                        dest='num_tags',
                        type=int,
                        required=True,
                        help='Number of tags in barcode')
    parser.add_argument('-d', '--dir',
                        dest='dir',
                        type=str,
                        action='store',
                        required=True,
                        help='Directory to write split bams')
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
                        required = False,
                        help = "The maximum cluster size to keep")

    opts = parser.parse_args()

    return opts


def main():
    args = parse_args()
    labels = assign_labels(args.clusters,args.min_oligos, args.proportion, args.max_size)
    label_bam_file(args.input_bam, args.output_bam, labels, args.num_tags)
    split_bam_by_RG(args.output_bam, args.dir)


def assign_labels(clusterfile, min_oligos, proportion, max_size):
    labels = {}
    with open(clusterfile, 'r') as clusters:
        for line in tqdm.tqdm(clusters):
            barcode, *reads = line.rstrip('\n').split('\t')
            multilabel = get_single_label(reads, min_oligos, proportion, max_size)
            labels[barcode] = multilabel
    return labels

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


def label_bam_file(input_bam, output_bam, labels, num_tags):
    count, duplicates, skipped = 0,0,0
    written = defaultdict(int)
    pattern = re.compile('::' + num_tags * '\[([a-zA-Z0-9_\-]+)\]')
    library = os.path.basename(input_bam).split('.')[0]
    header = construct_read_group_header(input_bam, labels)
    found = defaultdict(set)
    with pysam.AlignmentFile(input_bam, "rb") as in_bam, \
    pysam.AlignmentFile(output_bam, "wb", header=header) as out_bam:
        for read in in_bam.fetch(until_eof = True):
            count += 1
            if count % 1000000 == 0:
                print(count)
            name = read.query_name
            match = pattern.search(name)
            barcode = list(match.groups())[1:]
            barcode.append(library)
            full_barcode = '.'.join(barcode)
            position = read.reference_name + "_" + str(read.reference_start)
            if position in found[full_barcode]:
                duplicates +=1
            else:
                try:
                    readlabel = labels[full_barcode]
                    found[full_barcode].add(position)
                    read.set_tag('RG', readlabel, replace=True)
                    out_bam.write(read)
                    written[readlabel] +=1
                except KeyError:
                    skipped += 1

    print('Total reads:', count)
    print('Reads written:', written)
    print('Duplicate reads:', duplicates)
    print('Reads with an error not written out:', skipped)

def construct_read_group_header(input_bam, labels):
    proteins = set(labels.values())
    sample_name = input_bam.split('.',1)[0]
    with pysam.AlignmentFile(input_bam, 'rb') as input_file:
        bam_header = input_file.header.to_dict()
        read_group_dict = [{"ID":name, "SM":sample_name} for name in list(proteins)]
        bam_header["RG"] = read_group_dict
        return bam_header


def split_bam_by_RG(output_bam, output_dir):
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    format_string = output_dir + "/%*_%!.%."
    pysam.split("-f", format_string, output_bam)

if __name__== "__main__":
    main()

