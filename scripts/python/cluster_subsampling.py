import pandas as pd
import gzip
import argparse
import sys
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
    args = parser.parse_args()

    cluster_files = sorted(args.input)
    if not cluster_files:
        print("No input files provided.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(cluster_files)} files.")

    # Load each cluster file separately first
    separate_dataframes = []
    for file in tqdm(cluster_files, total=len(cluster_files), desc="Loading individual clusterfiles"):
        separate_dataframes.append(merge_extra_fields(file))

    # Concatenate the first 0...n dataframes from (1, n)
    concatenated_dataframes = {}
    progress_bar = tqdm(range(len(separate_dataframes)),
                        total=len(separate_dataframes),
                        desc="Concatenating dataframes")
    for n in progress_bar:
        key = f"chunk_0_to_{n}"
        concatenated_dataframes[key] = pd.concat(separate_dataframes[:n + 1], axis=0)

    cluster_stats = []
    read_stats = []

    def nonspecific_counter(x): return sum(x.values())

    for name, df in concatenated_dataframes.items():
        print(f"Processing {name}")

        # When there's only one chunk, there's nothing to aggregate
        if int(name.split("_")[-1]) == 0:
            grouped_df = df.copy()
        else:
            grouped_series = df.groupby("barcode").tags.agg(lambda x: [item for sublist in x for item in sublist])
            grouped_df = pd.DataFrame(grouped_series)

        grouped_df["cluster_category"] = grouped_df["tags"].apply(find_cluster_category)
        grouped_df["tag_count"] = grouped_df["tags"].apply(lambda x: Counter(x))

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

        cluster_stats.append({
            "bpm_only": bpm_only_cluster,
            "rpm_only": rpm_only_cluster,
            "bpm_and_rpm": bpm_and_rpm_cluster,
            "summed": bpm_only_cluster + rpm_only_cluster + bpm_and_rpm_cluster,
        })

        read_stats.append({
            "bpm_only": bpm_only_read,
            "rpm_only": rpm_only_read,
            "bpm_and_rpm": bpm_and_rpm_read,
            "summed": bpm_only_read + rpm_only_read + bpm_and_rpm_read,
        })

    # Optionally write summary table
    cluster_stats_df = pd.DataFrame(cluster_stats)
    read_stats_df = pd.DataFrame(read_stats)

    if args.plot_data:
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

    cluster_stats_df.to_csv(args.output + ".cluster_stats.txt", index=False)
    read_stats_df.to_csv(args.output + ".read_stats.txt", index=False)
    print(f"Stats table written to {args.output}")


if __name__ == "__main__":
    main()
