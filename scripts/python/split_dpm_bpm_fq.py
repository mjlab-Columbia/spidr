import gzip
import os
import argparse
import re

def parse_args():

    parser = argparse.ArgumentParser(description='Split fastq based on dpm or bead barcodes')
    parser.add_argument('--r1', dest='read_1', type=str, required=True,
                        help='Fastq read 1')
    opts = parser.parse_args()

    return opts


def main():

    opts = parse_args()

    read_1_path = opts.read_1

    dpm_out_path = os.path.splitext(os.path.splitext(read_1_path)[0])[0] + '_dpm.fastq.gz'
    bpm_out_path =  os.path.splitext(os.path.splitext(read_1_path)[0])[0] + '_bpm.fastq.gz'
    other_out_path =  os.path.splitext(os.path.splitext(read_1_path)[0])[0] + '_other.fastq.gz'
    short_out_path =  os.path.splitext(os.path.splitext(read_1_path)[0])[0] + '_short.fastq.gz'


    dpm_count = 0
    bpm_count = 0
    other_count = 0
    incomplete = 0
    counter = 0

    pattern = re.compile('\[([a-zA-Z0-9_\-]+)\]')

    with file_open(read_1_path) as read_1, \
    gzip.open(dpm_out_path, 'wt') as dpm_out, \
    gzip.open(bpm_out_path, 'wt') as bpm_out, \
    gzip.open(other_out_path, 'wt') as other_out, \
    gzip.open(short_out_path, 'wt') as short_out:
        for qname, seq, thrd, qual in fastq_parse(read_1):
                counter +=1
                barcodes = pattern.findall(qname)
                if counter % 10000 == 0:
                    print(counter)
                if 'NOT_FOUND' in barcodes:
                    incomplete += 1
                    short_out.write(qname + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')
                elif 'DPM' in qname:
                    dpm_count += 1
                    dpm_out.write(qname + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')
                elif 'BEAD' in qname:
                    bpm_count += 1
                    bpm_out.write(qname + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')
                else:
                    other_count += 1
                    other_out.write(qname + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')
                    

    print('Reads without full barcode:', incomplete)
    print('DPM reads out:', dpm_count)
    print('BPM reads out:', bpm_count)
    print('Reads without Bead or DPM barcode:', other_count)



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


def fastq_parse(fp):
    """
    Parse fastq file.
    """
    linecount = 0
    name, seq, thrd, qual = [None] * 4
    for line in fp:

        linecount += 1
        if linecount % 4 == 1:
            try:
                name = line.decode('UTF-8').rstrip()
            except AttributeError:
                name = line.rstrip()
            assert name.startswith('@'),\
                   "ERROR: The 1st line in fastq element does not start with '@'.\n\
                   Please check FastQ file near line number %s" % (linecount)
        elif linecount % 4 == 2:
            try:
                seq = line.decode('UTF-8').rstrip()
            except AttributeError:
                seq = line.rstrip()
        elif linecount % 4 == 3:
            try:
                thrd = line.decode('UTF-8').rstrip()
            except AttributeError:
                thrd = line.rstrip()
            assert thrd.startswith('+'),\
                   "ERROR: The 3st line in fastq element does not start with '+'.\n\
                   Please check FastQ file near line number %s" % (linecount)
        elif linecount % 4 == 0:
            try:
                qual = line.decode('UTF-8').rstrip()
            except AttributeError:
                qual = line.rstrip()
            assert len(seq) == len(qual),\
                    "ERROR: The length of Sequence and Quality aren't equal.\n\
                    Please check FastQ file near line number %s" % (linecount)
            yield name, seq, thrd, qual,
            name, seq, thrd, qual = [None] * 4


if __name__ == "__main__":
    main()
