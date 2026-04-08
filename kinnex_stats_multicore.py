#!/usr/bin/env python3
"""
kinnex_stats_multicore.py

Compute per-sample FLNC stats across one or more BAMs and write a CSV.
Multiple BAMs for the same sample are streamed and pooled before stats
are computed — no intermediate merged files are written.

Usage:
    # Single BAM:
    python kinnex_stats_multicore.py --bams sample.flnc.bam --sample MySample --out stats.csv

    # Multiple BAMs for one sample:
    python kinnex_stats_multicore.py --bams s1_c1.bam s1_c2.bam --sample MySample --out stats.csv

    # Legacy FOFN mode (single-sample, kept for back-compat):
    python kinnex_stats_multicore.py --fofn sample.fofn --sample MySample --out stats.csv

Requirements:
    pip install pysam

Authored by ChatGPT per prompts by K Munson 10/9/2025
Updated for multi-BAM-per-sample support 2025-10

"""

import argparse
import csv
import math
import os
import sys
import statistics
import pysam
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from typing import List, Dict, Any, Tuple

DEFAULT_WORKERS = 4
OUT_CSV = "stats.csv"


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def phred_from_rq(rq_raw) -> float:
    """Convert PacBio rq value to PHRED score."""
    try:
        rq = float(rq_raw)
    except Exception:
        raise ValueError(f"Cannot convert rq value: {rq_raw}")
    if rq > 1.0 and rq <= 1e6:
        rq = rq / 1e6
    if rq >= 1.0:
        rq = 1.0 - 1e-12
    if rq <= 0.0:
        return 0.0
    return -10 * math.log10(1 - rq)


def n50_from_list(values: List[float]) -> float:
    """Return N50 from a list of numbers."""
    if not values:
        return 0.0
    sorted_vals = sorted(values, reverse=True)
    half_total = sum(sorted_vals) / 2.0
    cumsum = 0.0
    for v in sorted_vals:
        cumsum += v
        if cumsum >= half_total:
            return v
    return sorted_vals[-1]


def compute_stats(sample: str, lengths: List[int], qualities: List[float]) -> Dict[str, Any]:
    """Compute summary stats from pooled length and quality lists."""
    result = {
        "Sample": sample,
        "n": 0,
        "n_bases": 0,
        "AvgLen": "NA", "MedLen": "NA", "N50Len": "NA",
        "AvgQ":   "NA", "MedQ":   "NA", "N50Q":   "NA",
        "error": None,
    }
    n = len(lengths)
    result["n"] = n
    if n == 0:
        return result

    result["n_bases"] = sum(lengths)
    result["AvgLen"] = round(result["n_bases"] / n, 2)
    result["MedLen"] = round(statistics.median(lengths), 2)
    result["N50Len"] = int(n50_from_list(lengths))

    if qualities:
        result["AvgQ"] = round(sum(qualities) / len(qualities), 2)
        result["MedQ"] = round(statistics.median(qualities), 2)
        result["N50Q"] = round(n50_from_list(qualities), 2)

    return result


# ---------------------------------------------------------------------------
# Per-BAM worker — returns raw lists, not summary stats
# ---------------------------------------------------------------------------

def collect_reads(bam_path: str, min_length: int = 0) -> Tuple[List[int], List[float], str | None]:
    """
    Stream one BAM and return (lengths, qualities, error_or_None).
    Returns raw lists so the caller can pool across multiple BAMs before
    computing any statistics — avoids large intermediate files.
    """
    lengths: List[int] = []
    qualities: List[float] = []

    if not os.path.exists(bam_path):
        return lengths, qualities, f"File not found: {bam_path}"

    try:
        with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
            for read in bam:
                # Skip secondary and supplementary; keep primary + unmapped
                if read.is_secondary or read.is_supplementary:
                    continue
                qlen = read.query_length
                if qlen is None:
                    seq = read.query_sequence
                    qlen = len(seq) if seq is not None else 0
                if min_length and qlen < min_length:
                    continue
                lengths.append(int(qlen))
                if read.has_tag("rq"):
                    try:
                        qualities.append(phred_from_rq(read.get_tag("rq")))
                    except Exception:
                        continue
    except Exception as exc:
        return lengths, qualities, str(exc)

    return lengths, qualities, None


# ---------------------------------------------------------------------------
# Multi-BAM sample processing
# ---------------------------------------------------------------------------

def process_sample(
    sample: str,
    bam_paths: List[str],
    min_length: int = 0,
    workers: int = DEFAULT_WORKERS,
    progress: bool = False,
) -> Dict[str, Any]:
    """
    Process all BAMs for a sample in parallel, pool the raw reads,
    then compute stats once. No intermediate files are written.
    """
    all_lengths: List[int] = []
    all_qualities: List[float] = []
    errors = []

    worker_fn = partial(collect_reads, min_length=min_length)
    n_bams = len(bam_paths)
    workers_to_use = min(workers, n_bams)

    if progress:
        print(f"[{sample}] Processing {n_bams} BAM(s) with {workers_to_use} worker(s)...",
              file=sys.stderr)

    with ProcessPoolExecutor(max_workers=workers_to_use) as exe:
        futures = {exe.submit(worker_fn, p): p for p in bam_paths}
        for fut in as_completed(futures):
            bam_path = futures[fut]
            lengths, qualities, error = fut.result()
            all_lengths.extend(lengths)
            all_qualities.extend(qualities)
            if error:
                errors.append(f"{os.path.basename(bam_path)}: {error}")
            if progress:
                print(
                    f"[{sample}] {os.path.basename(bam_path)} "
                    f"reads={len(lengths)}" + (f" ERROR={error}" if error else ""),
                    file=sys.stderr,
                )

    result = compute_stats(sample, all_lengths, all_qualities)
    if errors:
        result["error"] = "; ".join(errors)
    return result


# ---------------------------------------------------------------------------
# FOFN helper (back-compat)
# ---------------------------------------------------------------------------

def read_fofn(fofn_path: str) -> List[str]:
    paths = []
    with open(fofn_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                paths.append(line)
    return paths


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compute per-sample FLNC stats across one or more BAMs."
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--bams", nargs="+", metavar="BAM",
        help="One or more BAM paths for this sample",
    )
    input_group.add_argument(
        "--fofn", metavar="FOFN",
        help="File of BAM paths (one per line) — back-compat mode",
    )

    parser.add_argument(
        "--sample", required=True,
        help="Sample name written to the output CSV",
    )
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS,
                        help=f"Parallel BAM workers (default {DEFAULT_WORKERS})")
    parser.add_argument("--out", default=OUT_CSV,
                        help=f"Output CSV (default {OUT_CSV})")
    parser.add_argument("--min-length", type=int, default=0,
                        help="Skip reads shorter than this (default 0 = no filter)")
    parser.add_argument("--progress", action="store_true",
                        help="Print per-BAM progress to stderr")
    args = parser.parse_args()

    if args.bams:
        bam_paths = args.bams
    else:
        try:
            bam_paths = read_fofn(args.fofn)
        except Exception as e:
            print(f"Error reading FOFN: {e}", file=sys.stderr)
            sys.exit(1)

    if not bam_paths:
        print("No BAM paths found.", file=sys.stderr)
        sys.exit(1)

    result = process_sample(
        sample=args.sample,
        bam_paths=bam_paths,
        min_length=args.min_length,
        workers=args.workers,
        progress=args.progress,
    )

    if args.progress:
        print(
            f"[{args.sample}] total reads={result['n']}"
            + (f" ERRORS={result['error']}" if result.get("error") else ""),
            file=sys.stderr,
        )

    fieldnames = ["Sample", "n", "n_bases", "AvgLen", "MedLen", "N50Len", "AvgQ", "MedQ", "N50Q"]
    try:
        with open(args.out, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow({k: result.get(k, "NA") for k in fieldnames})
    except Exception as e:
        print(f"Error writing CSV: {e}", file=sys.stderr)
        sys.exit(1)

    if args.progress:
        print(f"Done. Wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
