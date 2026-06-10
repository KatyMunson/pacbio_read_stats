#!/usr/bin/env python3
"""
subsample_bam.py

Draw N reads from a real PacBio FLNC BAM and write a smaller BAM,
preserving the original header, rq tags, read lengths, and BGZF encoding.

Used by run_benchmark.py to create scaled-down test inputs, or standalone:

    python subsample_bam.py --bam /data/sample.flnc.bam --n 100000 --out sub_100k.bam

The --n count refers to primary reads only:
- For UBAMs (the common case) every record is primary, so --n equals total records written.
- For aligned BAMs, secondary and supplementary alignments are skipped so that --n counts
  distinct reads — consistent with how kinnex_stats_multicore.py counts reads.

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
    For UBAMs (the common case) every record is primary, so this filter is a no-op
    with negligible overhead (two bit-flag comparisons per read vs. BGZF decompression
    cost). For aligned BAMs the filter ensures n counts distinct reads, not alignment
    records, and the coordinate ordering of a sorted BAM does not affect sampling
    uniformity — reservoir sampling guarantees each primary read has equal selection
    probability regardless of stream order.

    Returns (list_of_read_strings, header_dict, total_primary_count).
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
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
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
        description="Subsample N primary reads from a PacBio FLNC BAM.",
        epilog=(
            "Note: --n counts primary reads only. For UBAMs (the common case) this "
            "equals total records. For aligned BAMs, secondary/supplementary alignments "
            "are excluded — consistent with how kinnex_stats_multicore.py counts reads."
        ),
    )
    parser.add_argument("--bam", required=True, help="Input BAM path")
    parser.add_argument("--n", type=int, required=True, help="Number of primary reads to sample")
    parser.add_argument("--out", required=True, help="Output BAM path")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    args = parser.parse_args()

    subsample(args.bam, args.n, args.out, seed=args.seed)


if __name__ == "__main__":
    main()
