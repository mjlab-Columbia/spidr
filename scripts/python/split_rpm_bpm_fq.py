import gzip
import os
import argparse
import re
from click import command, option, Path


@command()
@option("--input", "-i", type=Path(exists=True), default=None,
        show_default=True, help="Path to barcoded read 1 gzipped, fastq file")
@option("--rpm_output", "-r", type=Path(), default=None, show_default=True, help="Output path for RPM reads")
@option("--bpm_output", "-b", type=Path(), default=None, show_default=True, help="Output path for BPM reads")
@option("--short_output", "-b", type=Path(), default=None, show_default=True,
        help="Output path for barcodes with incomplete reads")
def main(input, rpm_output, bpm_output, short_output):
    rpm_count = 0
    bpm_count = 0
    incomplete = 0
    counter = 0

    pattern = re.compile('\\[([a-zA-Z0-9_\\-]+)\\]')

    with file_open(input) as read_1, \
            gzip.open(rpm_output, 'wt') as rpm_out, \
            gzip.open(bpm_output, 'wt') as bpm_out, \
            gzip.open(short_output, 'wt') as short_out:
        for qname, seq, thrd, qual in fastq_parse(read_1):
            counter += 1
            barcodes = pattern.findall(qname)
            if counter % 10000 == 0:
                print(counter)
            if 'NOT_FOUND' in barcodes[1:]:
                incomplete += 1
                short_out.write(qname + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')
            elif 'BEAD' in qname:
                bpm_count += 1
                bpm_out.write(qname + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')
            else:
                rpm_count += 1
                rpm_out.write(qname + '\n' + seq + '\n' + thrd + '\n' + qual + '\n')

    print('Reads without full barcode:', incomplete)
    print('RPM reads out:', rpm_count)
    print('BPM reads out:', bpm_count)


def file_open(filename):
    """
    Open as normal or as gzip
    Faster using zcat?
    """
    # does file exist?
    f = open(filename, 'rb')
    if (f.read(2) == b'\x1f\x8b'):  # compressed alsways start with these two bytes
        f.seek(0)  # return to start of file
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
            assert name.startswith('@'), \
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
            assert thrd.startswith('+'), \
                "ERROR: The 3st line in fastq element does not start with '+'.\n\
                   Please check FastQ file near line number %s" % (linecount)
        elif linecount % 4 == 0:
            try:
                qual = line.decode('UTF-8').rstrip()
            except AttributeError:
                qual = line.rstrip()
            assert len(seq) == len(qual), \
                "ERROR: The length of Sequence and Quality aren't equal.\n\
                    Please check FastQ file near line number %s" % (linecount)
            yield name, seq, thrd, qual,
            name, seq, thrd, qual = [None] * 4


if __name__ == "__main__":
    main()
