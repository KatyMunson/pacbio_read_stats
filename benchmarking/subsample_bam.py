#!/usr/bin/env python3
"""
subsample_bam.py

Draw N primary reads from a PacBio FLNC BAM and write a smaller BAM,
preserving the original header, rq tags, read lengths, and BGZF encoding.

Used by run_benchmark.py to create scaled-down test inputs, or standalone:

    python subsample_bam.py --bam /data/sample.flnc.bam --n 100000 --out sub_100k.bam

Requirements: pysam
"""

import argparse
import os
import random
import sys

import pysam


def reservoir_sample(bam_path: str, n: int, seed: int = 42) -> tuple:
    """
    Reservoir-sample n primary reads from a BAM in a single streaming pass.
    Skips secondary and supplementary alignments, matching collect_reads behaviour.
    Returns (list_of_reads, header_dict).
    """
    rng = random.Random(seed)
    reservoir = []
    count = 0

    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        header = bam.header.to_dict()
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
            count += 1
            if len(reservoir) < n:
                reservoir.append(read.to_string())
            else:
                j = rng.randint(0, count - 1)
                if j < n:
                    reservoir[j] = read.to_string()

    if count < n:
        print(
            f"Warning: BAM contains {count} primary reads, fewer than requested {n}. "
            "All reads will be written.",
            file=sys.stderr,
        )

    return reservoir, header, count


def write_subsampled_bam(read_strings: list, header: dict, out_path: str) -> None:
    """Write sampled reads (as SAM strings) to a new BAM using the original header."""
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    # Re-open source header via a temporary in-memory AlignmentHeader
    pysam_header = pysam.AlignmentHeader.from_dict(header)
    with pysam.AlignmentFile(out_path, "wb", header=pysam_header) as out_bam:
        for s in read_strings:
            read = pysam.AlignedSegment.fromstring(s, pysam_header)
            out_bam.write(read)


def subsample(bam_path: str, n_reads: int, out_path: str, seed: int = 42) -> str:
    """
    Subsample n_reads primary reads from bam_path and write to out_path.
    Skips generation if out_path already exists (cache-friendly for repeated runs).
    Returns out_path.
    """
    if os.path.exists(out_path):
        return out_path

    reads, header, total = reservoir_sample(bam_path, n_reads, seed=seed)
    write_subsampled_bam(reads, header, out_path)
    actual = min(n_reads, total)
    print(f"  Wrote {actual:,} reads → {out_path}", file=sys.stderr)
    return out_path


def main():
    parser = argparse.ArgumentParser(
        description="Subsample N primary reads from a PacBio FLNC BAM."
    )
    parser.add_argument("--bam", required=True, help="Input BAM path")
    parser.add_argument("--n", type=int, required=True, help="Number of reads to sample")
    parser.add_argument("--out", required=True, help="Output BAM path")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    args = parser.parse_args()

    subsample(args.bam, args.n, args.out, seed=args.seed)


if __name__ == "__main__":
    main()
