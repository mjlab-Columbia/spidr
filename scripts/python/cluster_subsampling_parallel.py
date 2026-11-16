"""
Cluster subsampling split into independent, cluster-parallel jobs.

This is the scatter/gather version of cluster_subsampling.py. The original script walks through cumulative sets of
cluster chunk files (chunks 0..0, 0..1, 0..2, ...) in a single process. Each cumulative set is independent of the
others, so here every sequencing depth is its own job:

  depth  --depth K <chunk files 0..K-1>
      Load the first K per-chunk cluster files (the subset of clusters), merge the barcodes that occur in more than
      one chunk (a barcode seen in several chunks is ONE cluster whose reads are pooled), classify every merged
      cluster as BPM_ONLY / RPM_ONLY / BPM_and_RPM and count clusters and reads per category. Writes one JSON file.

  aggregate  <depth json files>
      Gather the per-depth JSON files, sort them by depth, and write the same two tables as cluster_subsampling.py
      ({prefix}.{condition}.cluster_stats.txt and {prefix}.{condition}.read_stats.txt) plus the optional line plot,
      once per condition.

Conditions: clusters are split by the condition in their ROUND1 barcode (ROUND1_<CONDITION>_<well>, e.g.
ROUND1_BLANK_G11 -> BLANK) and subsampled independently; there is always at least one condition.

The numbers produced for depth K and a condition are identical to row K-1 of cluster_subsampling.py run on that
condition's clusters from the same chunk files.

Extra count: BPM_and_RPM clusters that also pass the threshold_tag_and_split filters (--min_oligos, --proportion,
--max_size) are reported separately as bpm_and_rpm_pass (not included in "summed"). This requires keeping the bead id
of every BPM read, not just the BPM count, so it is slower and uses more memory.

Memory note: tags are reduced to per-cluster counts (BPM tags, RPM tags, other tags) while a chunk is read, instead of
holding a Python list of every tag per barcode. Merging chunks is still a real merge by barcode (a groupby-sum over the
barcodes of all K chunks); only the bookkeeping per cluster is smaller.
"""
import argparse
import gzip
import json
import os
import re
import sys
from collections import Counter

import pandas as pd
from tqdm import tqdm


def open_auto(path):
    # Decide on magic bytes rather than extension, like the other cluster scripts
    with open(path, 'rb') as f:
        is_gzip = f.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt', encoding='utf-8', newline='') if is_gzip else open(
        path, 'rt', encoding='utf-8', newline='')


ROUND1_PATTERN = re.compile(r"^ROUND1_(.+)_[^_]+$")


def round1_condition(barcode):
    """Condition of a cluster = the middle part of its ROUND1 barcode (ROUND1_<CONDITION>_<well>, e.g. ROUND1_BLANK_G11
    -> BLANK). Barcode fields are '.'-separated; the ROUND1 field may sit at any index."""
    for field in barcode.split("."):
        m = ROUND1_PATTERN.match(field)
        if m:
            return m.group(1)
    sys.exit(f"No ROUND1 barcode found in cluster barcode: {barcode}")


def load_chunk(path):
    """Return (df, bead_counts) for one chunk.

    df: one row per cluster in this chunk, columns barcode, n_bpm, n_rpm, n_other.
    bead_counts: dict barcode -> Counter({bead_id: number of BPM reads with that bead id}), only for clusters that
        have at least one BPM read. Needed for the thresh_and_split-style filter, which depends on WHICH bead ids
        are present and not only on how many BPM reads there are.

    The cluster barcode drops the trailing library field, exactly like cluster_subsampling.merge_extra_fields.
    A tag's category is its first character (B = BPM[...], R = RPM[...]); every tag counts towards the read total.
    """
    barcodes, n_bpm, n_rpm, n_other = [], [], [], []
    bead_counts = {}
    with open_auto(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            barcode = ".".join(parts[0].split(".")[:-1])
            barcodes.append(barcode)
            b = r = 0
            beads = None
            for tag in parts[1:]:
                c = tag[0]
                if c == "B":
                    b += 1
                    # Bead id of a BPM tag, same parse as threshold_tag_and_split.get_single_label:
                    # 'BPM[]_RBFOX2:14122321-0' -> 'RBFOX2'
                    if beads is None:
                        beads = Counter()
                    beads[tag.split(":")[0].split("_", 1)[1]] += 1
                elif c == "R":
                    r += 1
            if beads is not None:
                bead_counts[barcode] = beads
            n_bpm.append(b)
            n_rpm.append(r)
            n_other.append(len(parts) - 1 - b - r)
    df = pd.DataFrame({"barcode": barcodes, "n_bpm": n_bpm, "n_rpm": n_rpm, "n_other": n_other})
    return df, bead_counts


def merge_chunks(frames, bead_count_dicts):
    """Merge clusters across chunks: identical barcodes become one cluster with all of their tags pooled.

    Returns the merged per-cluster count table plus two extra columns describing the bead ids of each cluster:
    top_bead_count (reads of the most common bead id) and total_bead_count (all BPM reads); both 0 without BPM.
    """
    allclusters = pd.concat(frames, axis=0, ignore_index=True)
    merged = allclusters.groupby("barcode", sort=False)[["n_bpm", "n_rpm", "n_other"]].sum()

    # Pool the per-bead-id counters of barcodes that occur in several chunks
    pooled = {}
    for d in bead_count_dicts:
        for barcode, beads in d.items():
            if barcode in pooled:
                pooled[barcode].update(beads)  # Counter.update adds counts
            else:
                pooled[barcode] = beads
    # Reduce each cluster's counter to the two numbers the filter needs and align to the merged table
    top = pd.Series({bc: max(c.values()) for bc, c in pooled.items()}, dtype="int64")
    merged["top_bead_count"] = top.reindex(merged.index, fill_value=0).to_numpy()
    merged["total_bead_count"] = merged["n_bpm"]  # every BPM read has exactly one bead id
    return merged


def passes_thresholds(merged, min_oligos, proportion, max_size):
    """Boolean mask of clusters that would be assigned to a single bead by threshold_tag_and_split.get_single_label.

    Same rules, in the same order:
      1. at least one BPM read (otherwise label 'none')
      2. cluster size = reads that are not BPM; more than max_size -> 'filtered'
      3. most common bead id has fewer than min_oligos reads -> 'uncertain'
      4. most common bead id is a smaller fraction of all BPM reads than proportion -> 'ambiguous'
    Surviving clusters get a bead label. (This is only the mask; combine it with has-RPM for BPM_and_RPM clusters.)
    """
    has_bpm = merged["n_bpm"] > 0
    cluster_size = merged["n_rpm"] + merged["n_other"]  # reads minus BPM reads
    not_filtered = cluster_size <= max_size
    enough_oligos = merged["top_bead_count"] >= min_oligos
    # clip avoids 0/0 for clusters without BPM (they are excluded by has_bpm anyway)
    majority = (merged["top_bead_count"] / merged["total_bead_count"].clip(lower=1)) >= proportion
    return has_bpm & not_filtered & enough_oligos & majority


def compute_stats(merged, pass_mask):
    """Cluster- and read-level category proportions for merged clusters (same definitions as cluster_subsampling.py)."""
    has_b = merged["n_bpm"] > 0
    has_r = merged["n_rpm"] > 0
    reads = merged["n_bpm"] + merged["n_rpm"] + merged["n_other"]

    masks = {
        "bpm_only": has_b & ~has_r,
        "rpm_only": ~has_b & has_r,
        "bpm_and_rpm": has_b & has_r,
    }
    # Separate extra count (NOT part of "summed"): BPM_and_RPM clusters that also pass the thresh_and_split filters
    pass_mask_bpm_rpm = has_b & has_r & pass_mask
    total_clusters = merged.shape[0]
    total_reads = int(reads.sum())

    cluster_row = {k: int(m.sum()) / total_clusters for k, m in masks.items()}
    # explicit left-to-right addition (builtin sum() uses compensated summation on Python >= 3.12 and would differ by 1 ulp)
    cluster_row["summed"] = cluster_row["bpm_only"] + cluster_row["rpm_only"] + cluster_row["bpm_and_rpm"]
    cluster_row["bpm_and_rpm_pass"] = int(pass_mask_bpm_rpm.sum()) / total_clusters

    abs_reads = {k: int(reads[m].sum()) for k, m in masks.items()}
    read_row = {k: abs_reads[k] / total_reads for k in masks}
    read_row["summed"] = (abs_reads["bpm_only"] + abs_reads["rpm_only"] + abs_reads["bpm_and_rpm"]) / total_reads
    for k in masks:
        read_row[f"{k}_abs"] = abs_reads[k]
    read_row["bpm_and_rpm_pass"] = int(reads[pass_mask_bpm_rpm].sum()) / total_reads
    read_row["bpm_and_rpm_pass_abs"] = int(reads[pass_mask_bpm_rpm].sum())
    read_row["total_reads"] = total_reads
    return cluster_row, read_row


def run_depth(args):
    files = sorted(args.input)
    if args.depth is not None and args.depth != len(files):
        sys.exit(f"--depth {args.depth} does not match the {len(files)} input chunk files given")
    loaded = [load_chunk(f) for f in tqdm(files, desc="Loading chunks")]
    merged = merge_chunks([df for df, _ in loaded], [beads for _, beads in loaded])
    merged["passes"] = passes_thresholds(merged, args.min_oligos, args.proportion, args.max_size)
    # Subsample independently per condition: each condition's clusters are summarised on their own
    conditions = merged.index.map(round1_condition)
    out = {"depth": len(files), "input_files": files, "conditions": {}}
    for cond in sorted(set(conditions)):
        sub = merged[conditions == cond]
        cluster_row, read_row = compute_stats(sub, sub["passes"])
        out["conditions"][cond] = {"n_clusters": int(sub.shape[0]), "cluster_row": cluster_row, "read_row": read_row}
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(out, fh, indent=2)
    summary = ", ".join(f"{c}={v['n_clusters']}" for c, v in out["conditions"].items())
    print(f"depth {len(files)}: merged clusters per condition: {summary} -> {args.output}")


def run_aggregate(args):
    results = []
    for p in args.input:
        with open(p) as fh:
            results.append(json.load(fh))
    results.sort(key=lambda r: r["depth"])
    depths = [r["depth"] for r in results]
    if len(set(depths)) != len(depths):
        sys.exit(f"Duplicate depths among inputs: {depths}")
    if depths != list(range(1, len(depths) + 1)):
        sys.exit(f"Depths are not a contiguous series 1..N: {depths}")

    all_conditions = sorted({c for r in results for c in r["conditions"]})
    if args.conditions:
        missing = set(args.conditions) - set(all_conditions)
        if missing:
            sys.exit(f"Expected conditions with no clusters in the depth files: {sorted(missing)}")
        extra = set(all_conditions) - set(args.conditions)
        if extra:
            print(f"Warning: conditions found in data but not expected: {sorted(extra)}", file=sys.stderr)
        all_conditions = list(args.conditions) + sorted(extra)

    for cond in all_conditions:
        # a condition with no clusters at a shallow depth is skipped in that depth's JSON; that is an error for a
        # contiguous series, so require every depth to contain it
        absent = [r["depth"] for r in results if cond not in r["conditions"]]
        if absent:
            sys.exit(f"Condition {cond} has no clusters at depths {absent}")
        prefix = f"{args.output}.{cond}"
        cluster_df = pd.DataFrame([r["conditions"][cond]["cluster_row"] for r in results])
        read_df = pd.DataFrame([r["conditions"][cond]["read_row"] for r in results])
        cluster_df.to_csv(prefix + ".cluster_stats.txt", index=False)
        read_df.to_csv(prefix + ".read_stats.txt", index=False)
        print(f"[{cond}] stats tables written to {prefix}.{{cluster,read}}_stats.txt ({len(results)} depths)")
        if args.plot_data:
            plot_read_stats(read_df, prefix + ".subsampling_line_plot.pdf", cond)
            print(f"[{cond}] plot written to {prefix}.subsampling_line_plot.pdf")


def plot_read_stats(read_stats_df, path, condition=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
    for col, label, color in (("bpm_only", "BPM only", "blue"),
                              ("rpm_only", "RPM only", "orange"),
                              ("bpm_and_rpm", "BPM and RPM", "green")):
        axs[0].plot(read_stats_df.index, read_stats_df[col], label=label, color=color, marker='o')
        axs[1].plot(read_stats_df.index, read_stats_df[col + "_abs"], label=label, color=color, marker='o')
    # Dashed line: BPM+RPM clusters that also pass the thresh_and_split filters (subset of the green line)
    axs[0].plot(read_stats_df.index, read_stats_df["bpm_and_rpm_pass"], label="BPM and RPM (pass filters)",
                color="green", linestyle="--", marker="x")
    axs[1].plot(read_stats_df.index, read_stats_df["bpm_and_rpm_pass_abs"], label="BPM and RPM (pass filters)",
                color="green", linestyle="--", marker="x")
    axs[0].set_xlabel('Number of Cumulative Chunks Aggregated')
    axs[0].set_ylabel('Proportion')
    axs[0].set_title('Cluster Category Proportion (By Total Reads)')
    axs[1].set_xlabel('Number of Cumulative Chunks Aggregated')
    axs[1].set_ylabel('Read Count')
    axs[1].set_title('Absolute Read Counts by Cluster Category')
    axs[1].yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    for ax in axs:
        ax.legend()
        ax.grid(True, linestyle='--', alpha=0.5)
    if condition:
        fig.suptitle(f'Condition: {condition}')
    plt.tight_layout()
    fig.savefig(path, format="pdf", dpi=300)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    d = sub.add_parser("depth", help="Compute stats for one cumulative depth (first K chunk files)")
    d.add_argument("--input", "-i", nargs="+", required=True, help="The K cluster chunk files making up this depth")
    d.add_argument("--depth", "-d", type=int, default=None,
                   help="Number of chunks expected (checked against --input; optional sanity check)")
    d.add_argument("--output", "-o", required=True, help="Output JSON file for this depth")
    d.add_argument("--min_oligos", type=int, required=True,
                   help="thresh_and_split min_oligos: minimum reads of the top bead id")
    d.add_argument("--proportion", type=float, required=True,
                   help="thresh_and_split proportion: minimum fraction of BPM reads from the top bead id")
    d.add_argument("--max_size", type=int, required=True,
                   help="thresh_and_split max_size: clusters with MORE non-BPM reads than this are filtered")
    d.set_defaults(func=run_depth)

    a = sub.add_parser("aggregate", help="Gather per-depth JSON files into the stats tables and plot")
    a.add_argument("--input", "-i", nargs="+", required=True, help="Per-depth JSON files from the depth step")
    a.add_argument("--output", "-o", required=True,
                   help="Prefix; writes {prefix}.{condition}.cluster_stats.txt etc. for every condition")
    a.add_argument("--conditions", "-c", nargs="+", default=None,
                   help="Expected conditions (ROUND1 middle tags); error if one has no clusters")
    a.add_argument("--plot_data", "-p", action="store_true", help="Also write {prefix}.subsampling_line_plot.pdf")
    a.set_defaults(func=run_aggregate)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
