#!/usr/bin/env python3

import json
import os
import re
from os.path import join, abspath
import argparse
from collections import defaultdict
import math

'''
Modified from https://github.com/crazyhottommy/pyflow-RNAseq/blob/master/fastq2json.py
'''

parser = argparse.ArgumentParser()
parser.add_argument("--fastq_dir", nargs='+', required=True,
                    help="Required. the FULL path to the fastq folder(s)")
parser.add_argument("--chunks", type=int, default=1,
                    help="Number of chunks to split the samples into (default: 1)")
args = parser.parse_args()

# default dictionary is quite useful!
FILES = defaultdict(lambda: defaultdict(list))

# build the dictionary with full path for each fastq.gz file
for folder in args.fastq_dir:
    for root, dirs, files in os.walk(folder):
        for f in files:
            if f.endswith("fastq.gz") or f.endswith("fq.gz") or f.endswith("fastq") or f.endswith("fq"):
                full_path = join(abspath(root), f)
                # R1 will be forward reads, R2 will be reverse reads
                m = re.search(r"(.+)_(R[12]).(fastq.gz|fq.gz|fastq|fq)", f)
                if m:
                    sample = m.group(1)
                    reads = m.group(2)
                    FILES[sample][reads].append(full_path)

all_samples = sorted(list(FILES.keys()))
total_samples = len(all_samples)

print()
print("total {} unique samples discovered".format(total_samples))
print("------------------------------------------")
for sample in all_samples:
    for read in FILES[sample]:
        print("{sample} {read} has {n} fastq".format(sample=sample, read=read, n=len(FILES[sample][read])))
print("------------------------------------------")

# --- Chunking Logic ---
if args.chunks < 1:
    raise ValueError("The number of chunks must be 1 or greater.")

if args.chunks == 1:
    # Original behavior if chunks is 1
    js = json.dumps(FILES, indent=4, sort_keys=True)
    with open('experiments.json', 'w') as f:
        f.write(js)
    print("Saved all samples to experiments.json")
else:
    # Divide samples into chunks cleanly
    # (Using a list comprehension to distribute remainders evenly if not perfectly divisible)
    chunk_size = total_samples / args.chunks
    chunks_list = [
        all_samples[int(round(chunk_size * i)): int(round(chunk_size * (i + 1)))]
        for i in range(args.chunks)
    ]
    
    # Write out each chunk to its respective file
    for index, chunk_samples in enumerate(chunks_list, start=1):
        # Skip creating empty files if there are more chunks requested than actual samples
        if not chunk_samples:
            continue
            
        chunk_dict = {sample: FILES[sample] for sample in chunk_samples}
        output_filename = "experiments_part{}.json".format(index)
        
        js = json.dumps(chunk_dict, indent=4, sort_keys=True)
        with open(output_filename, 'w') as f:
            f.write(js)
            
        print("Saved {} samples to {}".format(len(chunk_samples), output_filename))

print()
