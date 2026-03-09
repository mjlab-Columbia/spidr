import pandas as pd
import gzip
import argparse
import sys
import json
import os
from tqdm import tqdm
from collections import Counter
import matplotlib.pyplot as plt


def find_cluster_category(tag_array):
    if ("B" in tag_array) and ("R" in tag_array):
        return "BPM_and_RPM"
    elif ("B" in tag_array) and ("R" not in tag_array):
        return "BPM_ONLY"
    elif ("B" not in tag_array) and ("R" in tag_array):
        return "RPM_ONLY"
    else:
        return None


def merge_extra_fields(filepath):
    rows = []

    if filepath.endswith(".gz"):
        with gzip.open(filepath, 'rt', encoding='utf-8', newline='') as f:
            for line in f:
                parts = line.strip().split("\t")
                barcode = ".".join(parts[0].split(".")[:-1])
                tags = [tag[0] for tag in parts[1:]]
                rows.append((barcode, tags))
    else:
        with open(filepath, 'rt', encoding='utf-8', newline='') as f:
            for line in f:
                parts = line.strip().split("\t")
                barcode = ".".join(parts[0].split(".")[:-1])
                tags = [tag[0] for tag in parts[1:]]
                rows.append((barcode, tags))

    return pd.DataFrame(rows, columns=["barcode", "tags"])


def load_cache(cache_path):
    if cache_path and os.path.exists(cache_path):
        with open(cache_path, 'r') as f:
            return json.load(f)
    return {"processed_chunks": [], "input_files": []}


def save_cache(cache_path, cache_data):
    with open(cache_path, 'w') as f:
        json.dump(cache_data, f, indent=2)


def append_stat_row(filepath, row_dict, write_header):
    df = pd.DataFrame([row_dict])
    df.to_csv(filepath, mode='a', index=False, header=write_header)


def compute_chunk_stats(grouped_df):
    def nonspecific_counter(x): return sum(x.values())

    total_count = grouped_df.shape[0]
    bpm_only_cluster = grouped_df[grouped_df["cluster_category"] == "BPM_ONLY"].shape[0] / total_count
    rpm_only_cluster = grouped_df[grouped_df["cluster_category"] == "RPM_ONLY"].shape[0] / total_count
    bpm_and_rpm_cluster = grouped_df[grouped_df["cluster_category"] == "BPM_and_RPM"].shape[0] / total_count

    total_reads = grouped_df["tag_count"].apply(nonspecific_counter).sum()
    bpm_only_read = grouped_df[grouped_df["cluster_category"] == "BPM_ONLY"]["tag_count"].apply(
        nonspecific_counter).sum() / total_reads
    rpm_only_read = grouped_df[grouped_df["cluster_category"] == "RPM_ONLY"]["tag_count"].apply(
        nonspecific_counter).sum() / total_reads
    bpm_and_rpm_read = grouped_df[grouped_df["cluster_category"] == "BPM_and_RPM"]["tag_count"].apply(
        nonspecific_counter).sum() / total_reads

    cluster_row = {
        "bpm_only": bpm_only_cluster,
        "rpm_only": rpm_only_cluster,
        "bpm_and_rpm": bpm_and_rpm_cluster,
        "summed": bpm_only_cluster + rpm_only_cluster + bpm_and_rpm_cluster,
    }
    read_row = {
        "bpm_only": bpm_only_read,
        "rpm_only": rpm_only_read,
        "bpm_and_rpm": bpm_and_rpm_read,
        "summed": bpm_only_read + rpm_only_read + bpm_and_rpm_read,
    }
    return cluster_row, read_row


def main():
    parser = argparse.ArgumentParser(description="Compute proportions of cluster categories from cluster files.")
    parser.add_argument(
        "--input",
        "-i",
        nargs='+',
        required=True,
        help="Input files (provide one or more paths; separate with spaces)")
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="Prefix for output tables and optional PDF plot of data",
        default=None)
    parser.add_argument(
        "--plot_data",
        "-p",
        action='store_true',
        required=False,
        help="Whether to save an output plot showing the tables as line plots (optional)",
        default=None)
    parser.add_argument(
        "--cache",
        "-c",
        required=False,
        default=None,
        help="Path to a JSON cache file. If specified, processed chunks are tracked and skipped on re-runs, "
             "and stats are written incrementally after each chunk.")
    args = parser.parse_args()

    cluster_files = sorted(args.input)
    if not cluster_files:
        print("No input files provided.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(cluster_files)} files.")

    # Load cache if specified
    cache_data = load_cache(args.cache)
    processed_chunks = set(cache_data.get("processed_chunks", []))

    if args.cache:
        cached_input_files = cache_data.get("input_files", [])
        if cached_input_files and cached_input_files != cluster_files:
            print(
                "WARNING: Input files differ from those recorded in the cache. "
                "Cache may be invalid. Proceeding anyway.",
                file=sys.stderr)
        if processed_chunks:
            print(f"Cache found: skipping {len(processed_chunks)} already-processed chunk(s) "
                  f"({sorted(processed_chunks)}).")
        cache_data["input_files"] = cluster_files

    cluster_stats_path = args.output + ".cluster_stats.txt"
    read_stats_path = args.output + ".read_stats.txt"

    # Determine whether output files already have content (from prior cached runs)
    output_files_exist = (
        args.cache
        and bool(processed_chunks)
        and os.path.exists(cluster_stats_path)
        and os.path.exists(read_stats_path)
    )

    # Load each cluster file separately
    separate_dataframes = []
    for file in tqdm(cluster_files, total=len(cluster_files), desc="Loading individual clusterfiles"):
        separate_dataframes.append(merge_extra_fields(file))

    # Determine which chunks still need processing
    all_chunks = set(range(len(separate_dataframes)))
    pending_chunks = sorted(all_chunks - processed_chunks)

    if not pending_chunks:
        print("All chunks already processed. Nothing to do.")
    else:
        print(f"Processing {len(pending_chunks)} chunk(s): {pending_chunks}")

    cluster_stats = []
    read_stats = []

    # Build cumulative dataframe incrementally; only compute stats for pending chunks
    cumulative_df = None
    progress_bar = tqdm(range(len(separate_dataframes)),
                        total=len(separate_dataframes),
                        desc="Processing cumulative chunks")

    for n in progress_bar:
        if cumulative_df is None:
            cumulative_df = separate_dataframes[n].copy()
        else:
            cumulative_df = pd.concat([cumulative_df, separate_dataframes[n]], axis=0)

        if n in processed_chunks:
            continue

        chunk_name = f"chunk_0_to_{n}"
        progress_bar.set_postfix(chunk=chunk_name)

        # When there's only one chunk, there's nothing to aggregate across barcodes
        if n == 0:
            grouped_df = cumulative_df.copy()
        else:
            grouped_series = cumulative_df.groupby("barcode").tags.agg(
                lambda x: [item for sublist in x for item in sublist])
            grouped_df = pd.DataFrame(grouped_series)

        grouped_df["cluster_category"] = grouped_df["tags"].apply(find_cluster_category)
        grouped_df["tag_count"] = grouped_df["tags"].apply(lambda x: Counter(x))

        cluster_row, read_row = compute_chunk_stats(grouped_df)
        cluster_stats.append(cluster_row)
        read_stats.append(read_row)

        if args.cache:
            write_header = not output_files_exist and (n == pending_chunks[0])
            append_stat_row(cluster_stats_path, cluster_row, write_header)
            append_stat_row(read_stats_path, read_row, write_header)

            processed_chunks.add(n)
            cache_data["processed_chunks"] = sorted(processed_chunks)
            save_cache(args.cache, cache_data)

    # Write output (only when not using incremental cache writes)
    if not args.cache:
        cluster_stats_df = pd.DataFrame(cluster_stats)
        read_stats_df = pd.DataFrame(read_stats)
        cluster_stats_df.to_csv(cluster_stats_path, index=False)
        read_stats_df.to_csv(read_stats_path, index=False)

    print(f"Stats table written to {args.output}")

    if args.plot_data:
        cluster_stats_df = pd.read_csv(cluster_stats_path)
        read_stats_df = pd.read_csv(read_stats_path)

        # Panel plot: left for cluster stats, right for read stats
        fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharex=True)

        # Left panel: cluster stats
        axs[0].plot(cluster_stats_df.index, cluster_stats_df['bpm_only'], label='BPM only', color='blue', marker='o')
        axs[0].plot(cluster_stats_df.index, cluster_stats_df['rpm_only'], label='RPM only', color='orange', marker='o')
        axs[0].plot(
            cluster_stats_df.index,
            cluster_stats_df['bpm_and_rpm'],
            label='BPM and RPM',
            color='green',
            marker='o')
        axs[0].set_xlabel('Number of Cumulative Chunks Aggregated')
        axs[0].set_ylabel('Proportion')
        axs[0].set_title('Cluster Category Proportion (By Cluster)')
        axs[0].legend()
        axs[0].grid(True, linestyle='--', alpha=0.5)

        # Right panel: read stats
        axs[1].plot(read_stats_df.index, read_stats_df['bpm_only'], label='BPM only', color='blue', marker='o')
        axs[1].plot(read_stats_df.index, read_stats_df['rpm_only'], label='RPM only', color='orange', marker='o')
        axs[1].plot(read_stats_df.index, read_stats_df['bpm_and_rpm'], label='BPM and RPM', color='green', marker='o')
        axs[1].set_xlabel('Number of Cumulative Chunks Aggregated')
        axs[1].set_ylabel('Proportion')
        axs[1].set_title('Cluster Category Proportion (By Total Reads)')
        axs[1].legend()
        axs[1].grid(True, linestyle='--', alpha=0.5)

        plt.tight_layout()
        print(f"\nPlot written to table written to {args.output + '.subsampling_line_plot.pdf'}")
        fig.savefig(args.output + ".subsampling_line_plot.pdf", format="pdf", dpi=300)


if __name__ == "__main__":
    main()
