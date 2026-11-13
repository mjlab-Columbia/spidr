#!/usr/bin/env python3
"""Abbreviated QC table: one column per library, built from workup/qc/summary/<lib>.qc_summary.tsv.

Usage: qc_summary_abbrev.py [-d SUMMARY_DIR] [-o OUT_PREFIX] [LIB ...]
With no LIB, every *.qc_summary.tsv in SUMMARY_DIR is used. Writes OUT_PREFIX.tsv and OUT_PREFIX.md."""
import argparse
import glob
import os
import sys

p = argparse.ArgumentParser()
p.add_argument('libs', nargs='*')
p.add_argument('-d', '--summary_dir', default='workup/qc/summary')
p.add_argument('-o', '--out_prefix', default='workup/qc/qc_summary_abbrev')
a = p.parse_args()

libs = a.libs or sorted(os.path.basename(f)[:-len('.qc_summary.tsv')]
                        for f in glob.glob(f'{a.summary_dir}/*.qc_summary.tsv'))
if not libs:
    sys.exit(f'no *.qc_summary.tsv found in {a.summary_dir}')

D = {}
for lib in libs:
    r = {}
    for l in open(f'{a.summary_dir}/{lib}.qc_summary.tsv'):
        k, *v = l.rstrip('\n').split('\t')
        if k != 'metric':
            r[k] = v[0] if v else 'NA'
    D[lib] = r


def num(x):
    try:
        return float(x.replace(',', ''))
    except ValueError:
        return None


def fmt(n):
    return f'{n:,.0f}'


def g(l, k): return num(D[l].get(k, '0')) or 0.0
def S(l, *ks): return sum(g(l, k) for k in ks)


conds = sorted({k[len('bam_reads_pooled_'):] for l in libs for k in D[l] if k.startswith('bam_reads_pooled_')})
A = [
    ('Reads in', None, None),
    ('', 'raw R1 reads', lambda l: g(l, 'raw_fastq_reads')),
    ('', 'after adapter trimming (Trim Galore)', lambda l: g(l, 'filtered_fastq_reads')),
    ('', 'fully barcoded (R1 category)', lambda l: g(l, 'pre_alignment_fully_barcoded_reads_r1')),
    ('', 'bead (BPM) reads', lambda l: g(l, 'total_bpm_reads')),
    ('', 'cDNA (RPM) reads', lambda l: g(l, 'total_rpm_reads')),
    ('RNA adapter trimming (cutadapt)', None, None),
    ('', 'retained', lambda l: g(l, 'rpm_reads_post_rna_adapter_trim')),
    ('', 'too short', lambda l: g(l, 'rpm_reads_too_short')),
    ('', 'too long', lambda l: g(l, 'rpm_reads_too_long')),
    ('', 'untrimmed', lambda l: g(l, 'rpm_reads_untrimmed')),
    ('Alignment', None, None),
    ('', 'bowtie2 unique reads (of retained)', lambda l: g(l, 'bowtie2_unique_mapped_reads')),
    ('', 'STAR input reads', lambda l: g(l, 'star_input_reads')),
    ('', 'STAR unique', lambda l: g(l, 'star_uniquely_mapped_reads')),
    ('', 'STAR too many loci', lambda l: g(l, 'star_mapped_too_many_loci')),
    ('', 'STAR unmapped (mismatches + too short + other)',
     lambda l: S(l, 'star_unmapped_too_many_mismatches', 'star_unmapped_too_short', 'star_unmapped_other')),
    ('', 'STAR chimeric', lambda l: g(l, 'star_chimeric_reads')),
    ('', 'unique cDNA reads after dedup', lambda l: g(l, 'post_alignment_unique_rpm_reads')),
    ('Clusters and condition BAMs', None, None),
    ('', 'reads assigned to clusters', lambda l: g(l, 'barcoded_reads_assigned_to_clusters')),
    ('', 'BPMs in matching BPM+RPM clusters', lambda l: g(l, 'bpms_in_matching_clusters')),
    ('', 'threshold split: input reads', lambda l: g(l, 'thresh_split_input_reads')),
]
for c in conds:
    A += [
        ('', f'{c}: reads assigned to condition BAMs', lambda l, c=c: g(l, f'barcoded_reads_assigned_to_bams_{c}')),
        ('', f'{c}:   BAM reads: pooled', lambda l, c=c: g(l, f'bam_reads_pooled_{c}')),
        ('', f'{c}:   BAM reads: ambiguous', lambda l, c=c: g(l, f'bam_reads_ambiguous_{c}')),
        ('', f'{c}:   BAM reads: no cluster / none', lambda l, c=c: g(l, f'bam_reads_none_{c}')),
        ('', f'{c}:   BAM reads: filtered', lambda l, c=c: g(l, f'bam_reads_filtered_{c}')),
        ('', f'{c}: threshold split: duplicates', lambda l, c=c: g(l, f'thresh_split_duplicates_{c}')),
        ('', f'{c}: threshold split: no matching cluster', lambda l, c=c: g(l, f'thresh_split_no_matching_cluster_{c}')),
    ]
rates = [('bowtie2 overall alignment rate (%)', 'bowtie2_overall_alignment_rate_pct'),
         ('STAR unique rate (%)', 'star_uniquely_mapped_pct'), ('bead duplication rate', 'bpm_duplication_rate')]

os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
with open(a.out_prefix + '.tsv', 'w') as o, open(a.out_prefix + '.md', 'w') as m:
    o.write('\t'.join(['section', 'metric'] + libs) + '\n')
    m.write('| metric | ' + ' | '.join(libs) + ' |\n|---|' + '---:|' * len(libs) + '\n')
    sec = ''
    for s, lab, fn in A:
        if lab is None:
            sec = s
            m.write(f'| **{s}** |' + ' |' * len(libs) + '\n')
            continue
        vals = [fmt(fn(l)) for l in libs]
        if all(v == '0' for v in vals):
            continue
        o.write(f'{sec}\t{lab}\t' + '\t'.join(vals) + '\n')
        m.write(f'| {lab} | ' + ' | '.join(vals) + ' |\n')
    for lab, k in rates:
        vals = [D[l].get(k, 'NA') for l in libs]
        o.write('Rates\t' + lab + '\t' + '\t'.join(vals) + '\n')
        m.write(f'| {lab} | ' + ' | '.join(vals) + ' |\n')
print(f'{len(libs)} libraries -> {a.out_prefix}.tsv/.md')
