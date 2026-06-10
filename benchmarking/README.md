# Benchmarking the Kinnex Stats Pipeline

This folder contains tools for measuring real-world performance of
`kinnex_stats_multicore.py` and the Snakemake pipeline. Use the results to
set `threads`, `mem`, and `hrs` in `config.yaml` before submitting cluster jobs.

Benchmarks use **real subsampled BAMs** — subsets of your actual data — so
measurements reflect true read length distributions, BGZF compression ratios,
and I/O patterns.

---

## Contents

| File | Purpose |
|------|---------|
| `subsample_bam.py` | Draw N primary reads from a real BAM → smaller BAM (reservoir sampling, one pass, no index required) |
| `run_benchmark.py` | Run the pipeline across a parameter matrix; collect wall time, peak RSS, and throughput |

---

## Prerequisites

- Python 3.11 and pysam 0.22 (already in the project environment)
- Snakemake (only required for `snakemake` mode)
- At least one real PacBio Kinnex FLNC BAM file

Activate the project environment before running:

```bash
conda activate kinnex_stats
# or: module load python/3.11 pysam/0.22
```

---

## A note on UBAM vs. aligned BAM input

The majority use case for this pipeline is **UBAM (unaligned BAM)** — files
that come directly from PacBio CCS processing with no alignment step. In a
UBAM every record is a primary read, so `subsample_bam.py`'s filter for
secondary/supplementary alignments is a no-op (two bit-flag checks per read
that are always false — negligible overhead compared to BGZF decompression).

If you supply an **aligned BAM**, the subsampler skips secondary and
supplementary alignment records and counts only primary alignments toward
`--n`. This matches how `kinnex_stats_multicore.py` counts reads. The
coordinate ordering of a sorted aligned BAM does not bias the sample —
reservoir sampling guarantees each primary read has equal probability of
selection regardless of stream order.

In both cases, the `n_reported` value in the benchmark results CSV should
equal `n_reads_subsample × n_bams`.

---

## Quick Start

### Subsample a BAM manually

```bash
python subsample_bam.py \
    --bam /path/to/sample.flnc.bam \
    --n 50000 \
    --out /tmp/sub_50k.bam
```

### Run a small benchmark (~2–5 minutes)

```bash
cd benchmarking/

python run_benchmark.py standalone \
    --bam /path/to/sample.flnc.bam \
    --reads-list 10000 100000 \
    --bams-list 1 3 \
    --workers-list 1 4 \
    --results quick_results.csv
```

### Run the full parameter matrix (~30–90 minutes)

```bash
python run_benchmark.py standalone \
    --bam /path/to/sample.flnc.bam \
    --results full_results.csv
```

Default matrix: 3 read-count levels × 3 BAM-count levels × 4 worker counts = **36 combinations**.

### Benchmark the full Snakemake pipeline

```bash
python run_benchmark.py snakemake \
    --manifest /path/to/manifest.tsv \
    --cores 8 \
    --results snakemake_results.csv
```

This copies the Snakefile and config into a temporary directory, runs
`snakemake --cores N`, and reports per-sample wall time and throughput by
parsing the Snakemake log alongside the per-sample output CSVs.

---

## Disk Space

Subsampled BAMs are stored in a temporary directory and deleted after the
benchmark (use `--keep-bams` to retain them). Approximate sizes:

| Reads per BAM | Approx. BAM size |
|---------------|-----------------|
| 10,000        | ~30–80 MB        |
| 100,000       | ~300–800 MB      |
| 1,000,000     | ~3–8 GB          |

For the full matrix with 6 BAMs at 1M reads, you may need **18–48 GB** of
scratch space. Use `--out-dir` to direct scratch output to a partition with
sufficient space. Skip the 1M-read tier with `--reads-list 10000 100000` if
space is limited.

> **Note on subsampling time:** reservoir sampling streams the entire source
> BAM in a single pass regardless of how many reads are requested. Subsampling
> time is proportional to source BAM size, not `--n`. Use `--keep-bams` and
> `--out-dir` to cache subsampled BAMs across runs.

---

## How to Interpret Results

### Worker scaling

`kinnex_stats_multicore.py` is **I/O bound** — each BAM is streamed
sequentially, and multiple BAMs are read in parallel (one worker per BAM).
Speedup saturates once `workers ≥ n_bams`. Increasing workers beyond
`n_bams` has no effect because all files are already being read in parallel.
On NFS or shared HDD, contention between concurrent reads may reduce the
benefit of parallelism; local SSDs approach near-linear scaling.

### Memory scaling

Memory scales with **total reads across all BAMs**, not file size. The
`process_sample` function pools all read lengths and qualities into Python
lists in the parent process before computing statistics. Based on measured
results, expect roughly **29 MB per 100k reads** (~290 MB/M reads) for a
single BAM. With multiple BAMs the pooled total scales linearly with total
read count. Setting `mem` in `config.yaml` too low causes SGE to kill jobs;
the benchmark output's `peak_rss_mb` column gives the empirical value to use
as a floor.

### N50 computation cost

The O(n log n) sort inside `n50_from_list` is visible at high read counts
(≥1M) and becomes a meaningful fraction of total wall time above ~10M reads.
At 24M reads it added roughly 50s on top of BAM streaming time, which lowers
the apparent `reads_per_sec` figure even though streaming speed is unchanged.
The `reads_per_sec` column is therefore an end-to-end throughput rate, not a
pure I/O rate.

### Snakemake per-sample timing

The `snakemake` mode reports one row per sample with `wall_time_s` parsed
from the Snakemake log (second-resolution timestamps). Because samples run
in parallel up to `--cores`, the largest sample's `wall_time_s` determines
the overall walltime limit needed for cluster submissions. The `snakemake_total`
row shows the true end-to-end elapsed time including scheduler overhead
(typically 5–30 seconds).

### Setting `config.yaml` resources

Use benchmark results as follows:

| `config.yaml` parameter | How to set |
|--------------------------|------------|
| `threads` | Set to the **maximum number of BAMs any single sample has** — this matches the parallelism ceiling. Adding more threads beyond this wastes slot reservations. |
| `mem` (GB per thread) | Round `peak_rss_mb / threads` up to the next whole GB, then add 20% headroom for OS overhead. |
| `hrs` | Use the per-sample `wall_time_s` from the `snakemake` mode for your largest sample, then multiply by 1.5 as a safety buffer. |

---

## Measured Results

*Measured on a lab HPC cluster (qlogin on a compute node), source BAMs on NFS
and subsampled benchmark on local hard drive. Results will vary with disk type, CPU, and concurrent
I/O load. NFS latency typically reduces throughput vs. local SSD.*

### Standalone script benchmark

Subsampled from a single real Kinnex FLNC UBAM; multiple-BAM runs use
independently seeded subsamples to simulate distinct SMRT cells.

| reads/BAM | n_bams | workers | wall_time_s | peak_rss_mb | reads_per_sec | Mbases_per_sec |
|-----------|--------|---------|-------------|-------------|---------------|----------------|
| 10,000    | 1      | 1       | 0.80        | 19.5        | 12,453        | 23.2           |
| 10,000    | 1      | 4       | 0.30        | 18.6        | 33,223        | 62.0           |
| 10,000    | 3      | 1       | 3.31        | 18.4        | 9,061         | 16.9           |
| 10,000    | 3      | 4       | 0.30        | 20.8        | 99,668        | 185.4          |
| 100,000   | 1      | 1       | 1.50        | 29.0        | 66,489        | 123.2          |
| 100,000   | 1      | 4       | 1.30        | 28.3        | 76,746        | 142.2          |
| 100,000   | 3      | 1       | 3.61        | 50.9        | 83,195        | 154.2          |
| 100,000   | 3      | 4       | 1.50        | 51.0        | 199,601       | 370.0          |

Key observations:
- **Workers help most when `n_bams > 1`:** 3 BAMs × 4 workers is ~11× faster than 3 BAMs × 1 worker at 10k reads
- **Per-BAM startup cost** (~0.8s fixed overhead for pysam open + process fork) dominates at small read counts
- **Memory scales with total reads:** 100k × 3 BAMs → 51 MB ≈ 3× the single-BAM value (170 MB/M reads)

### Full pipeline benchmark (Snakemake)

Run against 4 real 1KG Kinnex samples with `--cores 8` (all 4 samples ran in parallel).

| sample      | n_reads    | wall_time_s | reads_per_sec | Mbases_per_sec |
|-------------|------------|-------------|---------------|----------------|
| CHM13-W-0   | 23,713,188 | 317         | 74,805        | 167.4          |
| GM07037-W-0 | 10,127,275 | 113         | 89,622        | 166.1          |
| GM11930-J-0 | 16,127,891 | 177         | 91,118        | 161.2          |
| GM11933-J-0 | 10,084,807 | 107         | 94,251        | 158.4          |
| **Total**   | 60,053,161 | 512 (wall)  |               |                |

Key observations:
- CHM13 has lower throughput (74k vs. 89–94k reads/s) likely due to more BAMs per sample
  causing more parallel NFS I/O contention when competing with the other three concurrent jobs
- All 4 samples ran simultaneously; total wall time (512s) ≈ CHM13's individual time (317s)
  plus ~200s Snakemake overhead and job scheduling on the shared node
- **For `config.yaml` on this cluster:** `hrs: 6` is very conservative; CHM13 at 317s
  with 1.5× buffer suggests `hrs: 1` is sufficient. `mem` should be validated against your
  largest sample's read count (see memory scaling above)

### Sanity-check bounds

If your measurements fall outside these ranges, check for disk contention,
NFS latency, or resource limits:

| Scenario | Expected wall time |
|----------|--------------------|
| 10k reads, 1 BAM, 1 worker | < 10 seconds |
| 100k reads, 1 BAM, 1 worker | < 30 seconds |
| 1M reads, 1 BAM, 1 worker | < 2 minutes |
| 10M reads, 1 BAM, 1 worker | < 10 minutes |
| 25M reads, 1 BAM, 1 worker | < 15 minutes |
| Memory per 1M total reads | ~170–290 MB RSS |
| Speedup at n_bams=3, workers=3 vs workers=1 | 3×–11× (NFS: lower end; local SSD: higher end) |
