#!/usr/bin/env python3
"""Aggregate high-level SPIDR QC metrics into a TSV per library (metrics as rows)."""

import json
import re
from pathlib import Path

import click
import pandas as pd


COUNT_METRIC_PREFIXES = (
    "barcoded_reads_assigned_to_bams_",
    "thresh_split_total_reads_",
)

COUNT_METRICS = {
    "raw_fastq_reads",
    "fastp_reads_before_filtering",
    "fastp_reads_after_filtering",
    "filtered_fastq_reads",
    "pre_alignment_fully_barcoded_reads_r1",
    "total_bpm_reads",
    "total_rpm_reads",
    "bowtie2_total_reads",
    "star_input_reads",
    "star_uniquely_mapped_reads",
    "star_mapped_multiple_loci",
    "star_mapped_too_many_loci",
    "star_unmapped_too_many_mismatches",
    "star_unmapped_too_short",
    "star_unmapped_other",
    "star_chimeric_reads",
    "post_alignment_unique_rpm_reads",
    "barcoded_reads_assigned_to_clusters",
}

PERCENTAGE_METRICS = {
    "bowtie2_overall_alignment_rate_pct",
    "star_uniquely_mapped_pct",
}

RATE_METRICS = {
    "fastp_duplication_rate",
    "fastp_q20_rate_after_filtering",
    "fastp_q30_rate_after_filtering",
    "fastp_gc_content_after_filtering",
    "bpm_duplication_rate",
}

STAR_COUNT_KEYS = {
    "star_input_reads": "Number of input reads",
    "star_uniquely_mapped_reads": "Uniquely mapped reads number",
    "star_mapped_multiple_loci": "Number of reads mapped to multiple loci",
    "star_mapped_too_many_loci": "Number of reads mapped to too many loci",
    "star_unmapped_too_many_mismatches": "Number of reads unmapped: too many mismatches",
    "star_unmapped_too_short": "Number of reads unmapped: too short",
    "star_unmapped_other": "Number of reads unmapped: other",
    "star_chimeric_reads": "Number of chimeric reads",
}


def is_count_metric(metric_name: str) -> bool:
    return metric_name in COUNT_METRICS or metric_name.startswith(COUNT_METRIC_PREFIXES)


def format_count(value) -> str:
    count = int(round(float(value)))
    if count < 0:
        raise ValueError(f"Read count must be a natural number, got {count}")
    if count >= 1000:
        return f"{count:,}"
    return str(count)


def format_display_value(metric_name: str, value) -> str:
    if value is None:
        return ""

    if is_count_metric(metric_name):
        return format_count(value)

    if metric_name in PERCENTAGE_METRICS:
        return f"{float(value):.2f}"

    if metric_name in RATE_METRICS:
        return f"{float(value):.4f}"

    if isinstance(value, float) and value.is_integer():
        return format_count(value)

    return str(value)


def read_integer(path: Path) -> int:
    return int(path.read_text().strip())


def read_float(path: Path) -> float:
    return float(path.read_text().strip())


def parse_fastp_json(path: Path) -> dict:
    with path.open() as handle:
        data = json.load(handle)

    read1_before = data["read1_before_filtering"]
    read1_after = data["read1_after_filtering"]
    duplication = data.get("duplication", {})

    def base_rates(read_section: dict) -> dict:
        total_bases = read_section["total_bases"]
        return {
            "q20_rate": read_section["q20_bases"] / total_bases,
            "q30_rate": read_section["q30_bases"] / total_bases,
        }

    after_rates = base_rates(read1_after)

    return {
        # fastp summary.total_reads counts R1 + R2; use read1_* for R1-only counts
        # that match raw_fastq_reads and downstream R1-based QC metrics.
        "fastp_reads_before_filtering": read1_before["total_reads"],
        "fastp_reads_after_filtering": read1_after["total_reads"],
        "fastp_duplication_rate": duplication.get("rate"),
        "fastp_q20_rate_after_filtering": after_rates["q20_rate"],
        "fastp_q30_rate_after_filtering": after_rates["q30_rate"],
        "fastp_gc_content_after_filtering": data["summary"]["after_filtering"]["gc_content"],
    }


def parse_bowtie2_qc(path: Path) -> dict:
    total_reads = 0
    weighted_aligned = 0.0
    chunk_reads = None

    with path.open() as handle:
        for line in handle:
            reads_match = re.match(r"^(\d+) reads;", line.strip())
            if reads_match:
                chunk_reads = int(reads_match.group(1))
                total_reads += chunk_reads
                continue

            rate_match = re.match(r"^([\d.]+)% overall alignment rate", line.strip())
            if rate_match and chunk_reads is not None:
                weighted_aligned += chunk_reads * float(rate_match.group(1)) / 100.0
                chunk_reads = None

    alignment_rate = (weighted_aligned / total_reads * 100.0) if total_reads else None
    return {
        "bowtie2_total_reads": total_reads,
        "bowtie2_overall_alignment_rate_pct": alignment_rate,
    }


def parse_star_log(path: Path) -> dict:
    metrics = {}
    with path.open() as handle:
        for line in handle:
            if "|" not in line:
                continue
            key, value = line.rsplit("|", 1)
            metrics[key.strip()] = value.strip()
    return metrics


def aggregate_star_logs(star_logs: list[Path]) -> dict:
    totals = {column: 0 for column in STAR_COUNT_KEYS}
    for log_path in star_logs:
        chunk_metrics = parse_star_log(log_path)
        for column, star_key in STAR_COUNT_KEYS.items():
            if star_key not in chunk_metrics:
                raise KeyError(f"Missing '{star_key}' in STAR log {log_path}")
            totals[column] += int(chunk_metrics[star_key])

    input_reads = totals["star_input_reads"]
    unique_reads = totals["star_uniquely_mapped_reads"]
    unique_pct = (unique_reads / input_reads * 100.0) if input_reads else None

    return {
        **totals,
        "star_uniquely_mapped_pct": unique_pct,
    }


def parse_condition_assignment(path: Path) -> tuple[str, int]:
    stem = path.name.removesuffix(".barcoded_reads_assigned_to_bams.txt")
    condition = stem.rsplit(".", 1)[1]
    return condition, read_integer(path)


def parse_thresh_split_log(path: Path) -> tuple[str, int]:
    marker = ".thresh_and_split_condition."
    if marker not in path.name:
        raise ValueError(f"Unexpected thresh-and-split log name: {path.name}")
    condition = path.name.split(marker, 1)[1].removesuffix(".log")
    with path.open() as handle:
        first_line = handle.readline().strip()
    match = re.match(r"Total reads:\s*(\d+)", first_line)
    if not match:
        raise ValueError(f"Could not parse total reads from {path}")
    return condition, int(match.group(1))


@click.command()
@click.option("--experiment", required=True, help="Library / aliquot identifier")
@click.option("--output", "-o", required=True, type=click.Path(), help="Output TSV path")
@click.option("--raw-reads", required=True, type=click.Path(exists=True))
@click.option("--fastp-json", required=True, type=click.Path(exists=True))
@click.option("--filtered-reads", required=True, type=click.Path(exists=True))
@click.option("--total-bpm-reads", required=True, type=click.Path(exists=True))
@click.option("--total-rpm-reads", required=True, type=click.Path(exists=True))
@click.option("--bowtie2-qc-log", required=True, type=click.Path(exists=True))
@click.option("--post-alignment-barcoded-count", required=True, type=click.Path(exists=True))
@click.option("--bpm-duplication-rate", required=True, type=click.Path(exists=True))
@click.option("--barcoded-reads-in-clusters", required=True, type=click.Path(exists=True))
@click.option(
    "--pre-alignment-barcoded-count",
    multiple=True,
    type=click.Path(exists=True),
    help="Per-chunk fully barcoded R1 count files",
)
@click.option(
    "--star-log",
    "star_logs",
    multiple=True,
    type=click.Path(exists=True),
    help="STAR Log.final.out files (one per chunk)",
)
@click.option(
    "--barcoded-in-bams",
    multiple=True,
    type=click.Path(exists=True),
    help="Per-condition barcoded read assignment count files",
)
@click.option(
    "--thresh-split-log",
    "thresh_split_logs",
    multiple=True,
    type=click.Path(exists=True),
    help="Per-condition threshold-and-split QC logs",
)
def main(
    experiment,
    output,
    raw_reads,
    fastp_json,
    filtered_reads,
    total_bpm_reads,
    total_rpm_reads,
    bowtie2_qc_log,
    post_alignment_barcoded_count,
    bpm_duplication_rate,
    barcoded_reads_in_clusters,
    pre_alignment_barcoded_count,
    star_logs,
    barcoded_in_bams,
    thresh_split_logs,
):
    metrics = {}

    metrics["raw_fastq_reads"] = read_integer(Path(raw_reads))
    metrics.update(parse_fastp_json(Path(fastp_json)))
    metrics["filtered_fastq_reads"] = read_integer(Path(filtered_reads))
    metrics["pre_alignment_fully_barcoded_reads_r1"] = sum(
        read_integer(Path(path)) for path in pre_alignment_barcoded_count
    )
    metrics["total_bpm_reads"] = read_integer(Path(total_bpm_reads))
    metrics["total_rpm_reads"] = read_integer(Path(total_rpm_reads))
    metrics.update(parse_bowtie2_qc(Path(bowtie2_qc_log)))
    metrics.update(aggregate_star_logs([Path(path) for path in star_logs]))
    metrics["post_alignment_unique_rpm_reads"] = read_integer(Path(post_alignment_barcoded_count))
    metrics["bpm_duplication_rate"] = read_float(Path(bpm_duplication_rate))
    metrics["barcoded_reads_assigned_to_clusters"] = read_integer(Path(barcoded_reads_in_clusters))

    for path in barcoded_in_bams:
        condition, count = parse_condition_assignment(Path(path))
        metrics[f"barcoded_reads_assigned_to_bams_{condition}"] = count

    for path in thresh_split_logs:
        condition, count = parse_thresh_split_log(Path(path))
        metrics[f"thresh_split_total_reads_{condition}"] = count

    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    formatted_values = [
        format_display_value(metric_name, value)
        for metric_name, value in metrics.items()
    ]

    summary_df = pd.DataFrame({"metric": list(metrics.keys()), experiment: formatted_values})
    summary_df.to_csv(output_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
