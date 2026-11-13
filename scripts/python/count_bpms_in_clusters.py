import click
from collections import Counter


def get_single_label(reads, min_oligos, proportion, max_size):
    """
    Same labelling as threshold_tag_and_split.get_single_label (kept here so this script does not need pysam).
    Returns the target name, or one of none / filtered / uncertain / ambiguous.
    """
    bead_reads = [read for read in reads if read.startswith("BPM")]
    if len(bead_reads) == 0:
        return "none"
    cluster_size = len(reads) - len(bead_reads)
    if cluster_size > max_size:
        return "filtered"
    bead_labels = Counter([read.split(":")[0].split("_", 1)[1] for read in bead_reads])
    candidate = bead_labels.most_common()[0]
    if candidate[1] < min_oligos:
        return "uncertain"
    elif candidate[1] / sum(bead_labels.values()) < proportion:
        return "ambiguous"
    return candidate[0]


NON_TARGET_LABELS = {"none", "filtered", "uncertain", "ambiguous"}


@click.command()
@click.option("--clusters", "-c", type=click.Path(exists=True), required=True, help="*.complete.clusters file")
@click.option("--output", "-o", type=click.Path(), required=True)
@click.option("--min_oligos", type=int, required=True)
@click.option("--proportion", type=float, required=True)
@click.option("--max_size", type=int, required=True)
def main(clusters, output, min_oligos, proportion, max_size):
    """
    Count BPM reads in clusters, by cluster type:
      bpms_in_clusters          all BPMs in all clusters
      bpms_in_bpm_only_clusters BPMs in clusters with no RPM (non-BPM) reads
      bpms_in_matching_clusters BPMs in clusters with >=1 RPM read that pass the min_oligos / proportion / max_size
                                cutoffs, i.e. whose RPM reads are written to a target BAM by thresh_and_split_condition
                                (not to the ambiguous / uncertain / filtered / none BAMs)
    Output: one "category<TAB>count" line per category.
    """
    counts = Counter()
    with open(clusters) as handle:
        for line in handle:
            _, *reads = line.rstrip("\n").split("\t")
            n_bpm = sum(read.startswith("BPM") for read in reads)
            counts["bpms_in_clusters"] += n_bpm
            if len(reads) == n_bpm:
                counts["bpms_in_bpm_only_clusters"] += n_bpm
            elif get_single_label(reads, min_oligos, proportion, max_size) not in NON_TARGET_LABELS:
                counts["bpms_in_matching_clusters"] += n_bpm
    with open(output, "w") as out:
        for key in ("bpms_in_clusters", "bpms_in_bpm_only_clusters", "bpms_in_matching_clusters"):
            out.write(f"{key}\t{counts[key]}\n")


if __name__ == "__main__":
    main()
