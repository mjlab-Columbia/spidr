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
      ({prefix}.cluster_stats.txt and {prefix}.read_stats.txt) plus the optional line plot.

The numbers produced for depth K are identical to row K-1 of cluster_subsampling.py run on the same chunk files.

Memory note: tags are reduced to per-cluster counts (BPM tags, RPM tags, other tags) while a chunk is read, instead of
holding a Python list of every tag per barcode. Merging chunks is still a real merge by barcode (a groupby-sum over the
barcodes of all K chunks); only the bookkeeping per cluster is smaller.
"""
import argparse
import gzip
import json
import os
import sys

import pandas as pd
from tqdm import tqdm


def open_auto(path):
    # Decide on magic bytes rather than extension, like the other cluster scripts
    with open(path, 'rb') as f:
        is_gzip = f.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rt', encoding='utf-8', newline='') if is_gzip else open(
        path, 'rt', encoding='utf-8', newline='')


def load_chunk(path):
    """Return a DataFrame (one row per cluster in this chunk) with columns barcode, n_bpm, n_rpm, n_other.

    The cluster barcode drops the trailing library field, exactly like cluster_subsampling.merge_extra_fields.
    A tag's category is its first character (B = BPM[...], R = RPM[...]); every tag counts towards the read total.
    """
    barcodes, n_bpm, n_rpm, n_other = [], [], [], []
    with open_auto(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            barcodes.append(".".join(parts[0].split(".")[:-1]))
            b = r = 0
            for tag in parts[1:]:
                c = tag[0]
                if c == "B":
                    b += 1
                elif c == "R":
                    r += 1
            n_bpm.append(b)
            n_rpm.append(r)
            n_other.append(len(parts) - 1 - b - r)
    return pd.DataFrame({"barcode": barcodes, "n_bpm": n_bpm, "n_rpm": n_rpm, "n_other": n_other})


def merge_chunks(frames):
    """Merge clusters across chunks: identical barcodes become one cluster with all of their tags pooled."""
    allclusters = pd.concat(frames, axis=0, ignore_index=True)
    return allclusters.groupby("barcode", sort=False)[["n_bpm", "n_rpm", "n_other"]].sum()


def compute_stats(merged):
    """Cluster- and read-level category proportions for merged clusters (same definitions as cluster_subsampling.py)."""
    has_b = merged["n_bpm"] > 0
    has_r = merged["n_rpm"] > 0
    reads = merged["n_bpm"] + merged["n_rpm"] + merged["n_other"]

    masks = {
        "bpm_only": has_b & ~has_r,
        "rpm_only": ~has_b & has_r,
        "bpm_and_rpm": has_b & has_r,
    }
    total_clusters = merged.shape[0]
    total_reads = int(reads.sum())

    cluster_row = {k: int(m.sum()) / total_clusters for k, m in masks.items()}
    # explicit left-to-right addition (builtin sum() uses compensated summation on Python >= 3.12 and would differ by 1 ulp)
    cluster_row["summed"] = cluster_row["bpm_only"] + cluster_row["rpm_only"] + cluster_row["bpm_and_rpm"]

    abs_reads = {k: int(reads[m].sum()) for k, m in masks.items()}
    read_row = {k: abs_reads[k] / total_reads for k in masks}
    read_row["summed"] = (abs_reads["bpm_only"] + abs_reads["rpm_only"] + abs_reads["bpm_and_rpm"]) / total_reads
    for k in masks:
        read_row[f"{k}_abs"] = abs_reads[k]
    read_row["total_reads"] = total_reads
    return cluster_row, read_row


def run_depth(args):
    files = sorted(args.input)
    if args.depth is not None and args.depth != len(files):
        sys.exit(f"--depth {args.depth} does not match the {len(files)} input chunk files given")
    frames = [load_chunk(f) for f in tqdm(files, desc="Loading chunks")]
    merged = merge_chunks(frames)
    cluster_row, read_row = compute_stats(merged)
    out = {
        "depth": len(files),
        "n_clusters": int(merged.shape[0]),
        "input_files": files,
        "cluster_row": cluster_row,
        "read_row": read_row,
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"depth {len(files)}: {out['n_clusters']} merged clusters -> {args.output}")


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

    cluster_df = pd.DataFrame([r["cluster_row"] for r in results])
    read_df = pd.DataFrame([r["read_row"] for r in results])
    cluster_df.to_csv(args.output + ".cluster_stats.txt", index=False)
    read_df.to_csv(args.output + ".read_stats.txt", index=False)
    print(f"Stats tables written to {args.output}.{{cluster,read}}_stats.txt ({len(results)} depths)")

    if args.plot_data:
        plot_read_stats(read_df, args.output + ".subsampling_line_plot.pdf")
        print(f"Plot written to {args.output}.subsampling_line_plot.pdf")


def plot_read_stats(read_stats_df, path):
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
    d.set_defaults(func=run_depth)

    a = sub.add_parser("aggregate", help="Gather per-depth JSON files into the stats tables and plot")
    a.add_argument("--input", "-i", nargs="+", required=True, help="Per-depth JSON files from the depth step")
    a.add_argument("--output", "-o", required=True, help="Prefix for {prefix}.cluster_stats.txt etc.")
    a.add_argument("--plot_data", "-p", action="store_true", help="Also write {prefix}.subsampling_line_plot.pdf")
    a.set_defaults(func=run_aggregate)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
