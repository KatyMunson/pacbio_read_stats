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

This copies the Snakefile and config into a temporary directory and runs
`snakemake --cores N`, measuring end-to-end wall time including scheduler
overhead.

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

---

## How to Interpret Results

### Worker scaling

`kinnex_stats_multicore.py` is **I/O bound** — each BAM is streamed
sequentially, and multiple BAMs are read in parallel (one worker per BAM).
Speedup saturates once `workers ≥ n_bams`. Increasing workers beyond
`n_bams` has no effect because all files are already being read in parallel.
On a shared HDD, contention between concurrent reads may reduce the benefit
of parallelism; SSDs approach near-linear scaling.

### Memory scaling

Memory scales with **total reads across all BAMs**, not file size. The
`process_sample` function pools all read lengths and qualities into Python
lists in the parent process before computing statistics. Expect roughly
100–200 MB RSS per 1 million reads. Setting `mem` in `config.yaml` too low
causes SGE to kill jobs; the benchmark output's `peak_rss_mb` column gives
the empirical value to use as a floor.

### N50 computation cost

The O(n log n) sort inside `n50_from_list` is visible at high read counts
(≥1M). This CPU cost is bounded and predictable; it adds a few seconds even
after BAM I/O completes.

### Snakemake overhead

Mode B adds Snakemake rule-scheduling overhead on top of per-sample compute
time, typically 5–30 seconds for the scheduler itself. The `wall_time_per_sample_s`
column in the Snakemake results CSV is the most useful number for estimating
walltime limits for cluster submissions.

### Setting `config.yaml` resources

Use benchmark results as follows:

| `config.yaml` parameter | How to set |
|--------------------------|------------|
| `threads` | Set to the **maximum number of BAMs any single sample has** — this matches the parallelism ceiling. Adding more threads beyond this wastes slot reservations. |
| `mem` (GB per thread) | Round `peak_rss_mb / threads` up to the next whole GB, then add 20% headroom for OS overhead. |
| `hrs` | Use the standalone benchmark wall time for your largest sample (most reads × most BAMs), then multiply by 1.5 as a safety buffer. |

---

## Example Results

Run the benchmark with your data and paste the output table here.

```bash
python run_benchmark.py standalone \
    --bam /path/to/sample.flnc.bam \
    --results full_results.csv
cat full_results.csv
```

Replace the placeholder table below with your measured values.

| reads_per_bam | n_bams | workers | wall_time_s | peak_rss_mb | reads_per_sec | mbases_per_sec |
|---------------|--------|---------|-------------|-------------|---------------|----------------|
| 10,000        | 1      | 1       | TBD         | TBD         | TBD           | TBD            |
| 10,000        | 1      | 4       | TBD         | TBD         | TBD           | TBD            |
| 10,000        | 3      | 1       | TBD         | TBD         | TBD           | TBD            |
| 10,000        | 3      | 4       | TBD         | TBD         | TBD           | TBD            |
| 100,000       | 1      | 1       | TBD         | TBD         | TBD           | TBD            |
| 100,000       | 3      | 1       | TBD         | TBD         | TBD           | TBD            |
| 100,000       | 3      | 4       | TBD         | TBD         | TBD           | TBD            |
| 100,000       | 6      | 4       | TBD         | TBD         | TBD           | TBD            |
| 1,000,000     | 1      | 1       | TBD         | TBD         | TBD           | TBD            |
| 1,000,000     | 3      | 4       | TBD         | TBD         | TBD           | TBD            |
| 1,000,000     | 6      | 4       | TBD         | TBD         | TBD           | TBD            |
| 1,000,000     | 6      | 8       | TBD         | TBD         | TBD           | TBD            |

*Measured on: [hardware description — CPU model, RAM, disk type (HDD/SSD/NFS)].
Results will vary with disk type, CPU speed, and concurrent I/O load.*

### Sanity-check bounds

If your measurements fall outside these ranges, check for disk contention,
NFS latency, or resource limits:

| Scenario | Expected wall time |
|----------|--------------------|
| 10k reads, 1 BAM, 1 worker | < 30 seconds |
| 100k reads, 1 BAM, 1 worker | < 5 minutes |
| 1M reads, 1 BAM, 1 worker | < 30 minutes |
| Memory per 1M total reads | < 500 MB RSS |
| Speedup at n_bams=6, workers=6 vs workers=1 | 2×–6× (HDD: lower end; SSD: higher end) |
