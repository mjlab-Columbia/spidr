from collections import Counter, defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import os
import gzip
import re

mpl.use("Agg")

ROUND1_PATTERN = re.compile(r"^ROUND1_(.+)_[^_]+$")

"""
Generate maximum representation ecdfs for bead type representation within clusters.
"""


def main():
    args = parse_arguments()

    # One pass over the clusters; per-condition (ROUND1) and overall values for both plots
    props, counts = collect_max_representation(args.clusterfile)
    clustername = os.path.basename(args.clusterfile).replace(".clusters.gz", "")

    ax = plot_ecdfs(props, clustername)
    ax.set(xlabel="Maximum Bead Representation Proportion", ylabel="Proportion of Beads")
    ax.get_figure().savefig(args.output_ecdf, bbox_inches="tight")
    plt.close(ax.get_figure())

    ax = plot_ecdfs(counts, clustername)
    ax.set(xlabel="Num Oligos for Most Common Type", ylabel="Proportion of Beads")
    ax.set(xlim=(0, int(args.xlim)))
    ax.get_figure().savefig(args.output_counts, bbox_inches="tight")
    plt.close(ax.get_figure())


def round1_condition(barcode):
    """Condition = middle part of the ROUND1 barcode field (ROUND1_<CONDITION>_<well>, e.g. ROUND1_BLANK_G11 -> BLANK).
    Same rule as cluster_subsampling_parallel.py."""
    for field in barcode.split("."):
        m = ROUND1_PATTERN.match(field)
        if m:
            return m.group(1)
    raise ValueError(f"No ROUND1 barcode found in cluster barcode: {barcode}")


def collect_max_representation(clusterfile):
    """
    For each non-singleton cluster, the max bead-type proportion and count among its BPM reads.

    Returns:
        (props, counts): dicts {condition: [values]}, each also holding an "All" key with every cluster
    """
    props = defaultdict(list)
    counts = defaultdict(list)
    with open_auto(clusterfile) as clusters:
        for line in clusters:
            barcode, *reads = line.rstrip("\n").split("\t")
            bead_reads = [read for read in reads if read.startswith("BPM")]
            if len(bead_reads) > 1:  # ignore singleton clusters
                bead_labels = Counter(
                    [read.split(":")[0].split("_", 1)[1] for read in bead_reads]
                )
                top = bead_labels.most_common()[0][1]
                cond = round1_condition(barcode)
                for key in (cond, "All"):
                    props[key].append(top / len(bead_reads))
                    counts[key].append(top)
    return props, counts


def plot_ecdfs(values_by_condition, clustername):
    """One ECDF line per condition (sorted), plus the all-cluster line in black dashes."""
    fig, ax = plt.subplots()
    for cond in sorted(k for k in values_by_condition if k != "All"):
        sns.ecdfplot(values_by_condition[cond], linewidth=3, ax=ax,
                     label=f"{cond} (n={len(values_by_condition[cond])})")
    if "All" in values_by_condition:
        sns.ecdfplot(values_by_condition["All"], linewidth=2, linestyle="--", color="black", ax=ax,
                     label=f"All (n={len(values_by_condition['All'])})")
    ax.set_title(clustername)
    ax.legend(title="ROUND1 condition")
    return ax


def open_auto(filename, mode='rt', *args, **kwargs):
    with open(filename, 'rb') as f:
        is_gzip = f.read(2) == b'\x1f\x8b'

    opener = gzip.open if is_gzip else open
    return opener(filename, mode, *args, **kwargs)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate the maximum representation ecdf plots to check for bead type uniqueness within clusters."
    )
    parser.add_argument(
        "--clusterfile",
        metavar="FILE",
        action="store",
        required=True,
        help="Path to a single clusters file",
    )
    parser.add_argument(
        "--output_ecdf",
        metavar="FILE",
        action="store",
        required=True,
        help="Output path for the maximum representation ECDF plot",
    )
    parser.add_argument(
        "--output_counts",
        metavar="FILE",
        action="store",
        required=True,
        help="Output path for the maximum representation counts plot",
    )
    parser.add_argument(
        "--xlim",
        action="store",
        required=True,
        help="The maximum x value on the counts plot",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
