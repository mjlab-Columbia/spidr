import pysam
import re
import gzip
import os
import sys
from collections import defaultdict, Counter
from tqdm import tqdm


class Position:
    """This class represents a genomic position, with type of nucleic acid (RNA or DNA)

    Methods:
    - to_string(): Returns a string representation of this position in the form
      "R/DPM(feature)_chrX:1000"
    """
    def __init__(self, read_type, feature, chromosome, start_coordinate, end_coordinate):
        self._type = read_type
        self._feature = feature
        self._chromosome = chromosome
        self._start_coordinate = start_coordinate
        self._end_coordinate = end_coordinate
    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        return (self._type == other._type and
                self._chromosome == other._chromosome and
                self._start_coordinate == other._start_coordinate)
    def __hash__(self):
        return hash((self._type, self._chromosome, 
                     self._start_coordinate, self._end_coordinate))
    def to_string(self):
        try:
            out = self._type + "[" + self._feature + "]" + "_" + \
                self._chromosome + ":" + \
                str(self._start_coordinate) + "-" + str(self._end_coordinate)
        except:
            print(self._type, self._feature, self._chromosome)
            print('Elements are not as expect!')
            sys.exit()
        return out
    def score(self):
        if self._type == 'RPM':
            if self._chromosome.startswith('chr'):
                return 2
            else:
                return 3
        elif self._type == 'BPM':
             return 1
        else:
             return 4
                
class Cluster:
    """This class represents a barcoding cluster as a collection of genomic
    positions.

    The underlying data structure is a set, so duplicate positions are
    discarded.

    Methods:
    - add_position(position): Adds a genomic position to this cluster

    - size(): Returns the number of reads or positions in this cluster

    - to_string(): Returns a string representation of this cluster as a
      tab-delimtited series of positions. See Position#to_string for how
      positions are represented as strings.

    - to_list(): Returns the Position class as a list (similar to_string())
    """

    def __init__(self):
        self._positions = set()
        self._label = ""
    def __iter__(self):
        return iter(self._positions)
    def change_label(self, label):
        self._label = label
    def get_label(self):
        return self._label
    def add_position(self, position):
        self._positions.add(position)
    def size(self, read_type=None):
        if read_type == None:
            return len(self._positions)
        else:
            return sum([1 if pos._type == read_type else 0 for pos in self._positions])
    def to_string(self):
        positions_sorted = sorted(list(self._positions), key = lambda x: x.score())
        position_strings = [position.to_string() for position in positions_sorted]
        return "\t".join(position_strings)
    def to_list(self):
        positions_sorted = sorted(list(self._positions), key = lambda x: x.score())
        position_strings = [position.to_string() for position in positions_sorted]
        return position_strings
    def count_type(self):
        count = Counter([position._type for position in self._positions])
        return count
    def get_coverage(self, read_type =None):
        if read_type == None:
            return Counter([pos._chromosome for pos in list(self._positions)])
        else:
            subset = [pos for pos in list(self._positions) if pos._type == read_type]
            return Counter([pos._chromosome for pos in subset])
    def contains_bead(self, name):
        labels = [pos._chromsome for pos in list(self._positions)]
        if name in labels:
            return True
        return False



class Clusters:
    """This class represents a collection of barcoding clusters.

    Methods:
    - get_cluster(barcode): Returns the cluster that corresponds to the given
      barcode. If the cluster does not exist, it is initialized (with zero
      positions), and this empty cluster is returned.

    - get_items(): iterate over clusters dictionary yielding keys and values

    - add_position(barcode, position): Adds the position to the cluster
      that corresponds with the given barcodes

    - to_strings(): Returns an iterator over the string representations of all
      of the contained clusters.

    - remove_cluster(barcode): Removes a cluster with the specified barcode

    - unique(): keep only unique cluster entries

    - make_lookup(): make a lookup table for converting cluster back into bam
    """
    def __init__(self):
        self._clusters = {}
    def __iter__(self):
        return iter(self._clusters.values())
    def __getitem__(self, barcode):
        return self._clusters[barcode]
    def get_cluster(self, barcode):
        if barcode not in self._clusters:
            self._clusters[barcode] = Cluster()
        return self._clusters[barcode]
    def add_cluster(self, barcode, cluster):
        self._clusters[barcode] = cluster
    def get_items(self):
        return self._clusters.items()
    def add_position(self, barcode, position):
        self.get_cluster(barcode).add_position(position)
    def to_strings(self):
        for barcode, cluster in self._clusters.items():
            yield barcode + '\t' + cluster.to_string()
    def to_labeled_strings(self):
        for barcode, cluster in self._clusters.items():
            yield cluster._label + '.' + barcode + '\t' + cluster.to_string()
    def remove_cluster(self, barcode):
        del self._clusters[barcode]
    def unique(self):
        for barcode, cluster in self._clusters.items():
            yield barcode + '\t' + cluster.unique()
    def make_lookup(self):
        lookup = defaultdict(set)
        for barcode, cluster in self._clusters.items():
            lookup[barcode].update(cluster.to_list())
        return lookup
    def make_stripped_lookup(self):
        lookup = defaultdict(set)
        for barcode, cluster in self._clusters.items():
            barcode_strip = barcode.split('.')[:-1]
            barcode_merge = '.'.join(barcode_strip)
            lookup[barcode_merge].update(cluster.to_list())
        return lookup


def get_clusters(filelist, num_tags):
    clusters = Clusters()
    pattern = re.compile('::' + num_tags * '\[([a-zA-Z0-9_\-]+)\]')
    dpm_counts = 0
    rpm_counts = 0
    bpm_counts = 0
    for sample in filelist:
        file_name = os.path.basename(sample)
        sample_name = file_name.split('.')[0]
        try:
            with pysam.AlignmentFile(sample, "rb") as f:
                for read in f.fetch(until_eof = True):
                    name = read.query_name
                    match = pattern.search(name)
                    barcode = list(match.groups())
                    read_type = barcode[0]
                    barcode_drop = barcode[1:]
                    if 'DPM' in read_type:
                        dpm_counts +=1
                        strand = '+' if not read.is_reverse else '-'
                        position = Position('DPM', strand, read.reference_name,
                                            read.reference_start, read.reference_end)
                    elif 'RPM' in read_type:
                        rpm_counts +=1
                        strand = '+' if not read.is_reverse else '-'
                        position = Position('RPM', strand, read.reference_name,
                                            read.reference_start, read.reference_end)
                    elif 'BEAD' in read_type:
                        bpm_counts +=1
                        UMI = read.reference_start
                        position = Position('BPM', '', read.reference_name, UMI, 0)
                    barcode_drop.append(sample_name)
                    barcode_str = ".".join(barcode_drop)
                    clusters.add_position(barcode_str, position)
        except ValueError:
            print('File provided has issues')
    print("Total BPM: ", bpm_counts)
    print("Total DPM: ", dpm_counts)
    print("Total RPM: ", rpm_counts)
    return clusters


def label_cluster_reads(reads, threshold, min_oligos, max_size):
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
     elif candidate[1]/sum(bead_labels.values()) < threshold:
         return 'ambiguous'
     else:  
         return candidate[0]

     return 'malformed'

def write_clusters_to_file(clusters, outfile, labeled=False):
    """Writes a Clusters object to a file"""

    count = 0
    with open(outfile, 'w') as f:

        if labeled:
            for cluster_string in clusters.to_labeled_strings():
                f.write(cluster_string)
                f.write("\n")
                count += 1
        else:
            for cluster_string in clusters.to_strings():
                f.write(cluster_string)
                f.write("\n")
                count += 1
    print('Number of clusters written: ',count)


def file_open(filename):
    """
    Open as normal or as gzip
    Faster using zcat?
    """
    #does file exist?
    f = open(filename,'rb')
    if (f.read(2) == b'\x1f\x8b'): #compressed alsways start with these two bytes
        f.seek(0) #return to start of file
        return gzip.GzipFile(fileobj=f, mode='rb')
    else:
        f.seek(0)
        return f

def parse_clusters(c_file, labeled=False):
    '''
    Parse cluster file

    Args:
        c_file(str): input path of cluster file
    '''

    total_reads = 0
    clusters = Clusters()
    pattern = re.compile('([a-zA-Z0-9]+)\[(.*)\]_(.+):([0-9]+)\-([0-9]+)')

    with file_open(c_file) as c:
        for line in tqdm(c):
            barcode, *reads = line.decode('utf-8').rstrip('\n').split('\t')
            if labeled:
                label, barcode = barcode.split('.', 1)
                clusters.get_cluster(barcode).change_label(label)
            for read in reads:
                total_reads += 1
                try:
                    match = pattern.search(read)
                    read_type, feature, chrom, start, end = match.groups()
                    position = Position(read_type, feature, chrom, start, end)
                    clusters.add_position(barcode, position)
                except:
                    print(read)
                    raise Exception('Pattern did not match above printed string')
    print('Total cluster reads:', total_reads)
    return(clusters)


def write_single_cluster(barcode, reads, out):
    positions_sorted = sorted(list(reads), key = lambda x: x.score())
    position_strings = [position.to_string() for position in positions_sorted]
    out_string =  "\t".join([barcode] + position_strings)
    out.write(out_string)
    out.write('\n')

def merge_clusters(in_file, out_file):
    current_barcode = ""
    current_reads = set()
    count = 0
    pattern = re.compile('([a-zA-Z0-9]+)\[(.*)\]_(.+):([0-9]+)\-([0-9]+)')
    with open(in_file, 'r') as in_clusters, \
    open(out_file, 'w') as out_clusters:
        for line in in_clusters:
             barcode, *reads = line.rstrip('\n').split('\t')
             if (barcode != current_barcode):
                  if current_barcode != "":
                      write_single_cluster(current_barcode, current_reads, out_clusters)
                      count += 1
                  current_barcode = barcode
                  current_reads = set()
             for read in reads:
                try:
                    match = pattern.search(read)
                    read_type, feature, chrom, start, end = match.groups()
                    position = Position(read_type, feature, chrom, start, end)
                    current_reads.add(position)
                except:
                    print(read)
                    raise Exception('Pattern did not match above printed string')
        write_single_cluster(current_barcode, current_reads, out_clusters)
        count +=1
    print("Total clusters written: ", count)
