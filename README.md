# pacbio_read_stats

Snakemake pipeline for computing per-sample read statistics from PacBio Kinnex FLNC BAMs. Each sample may span one or more BAMs (e.g. multiple SMRT cells); all BAMs are streamed and pooled before statistics are computed — no merged BAMs or intermediate FASTQ files are written to disk.

## Output

A single merged CSV at `results/kinnex_stats.csv` (path configurable) with one row per sample:

| Column | Description |
|--------|-------------|
| `Sample` | Sample name from manifest |
| `n` | Total read count |
| `n_bases` | Total bases |
| `AvgLen` | Mean read length |
| `MedLen` | Median read length |
| `N50Len` | Read length N50 |
| `AvgQ` | Mean read quality (PHRED, from `rq` tag) |
| `MedQ` | Median read quality |
| `N50Q` | Read quality N50 |

Per-sample CSVs are retained in `results/per_sample/` so that adding new samples to the manifest only requires processing the new samples — existing results are reused and the final table is regenerated.

## Requirements

- Snakemake
- Python 3.11
- pysam 0.22

Load via modules (`envmodules` directives are already in the Snakefile):
```
module load python/3.11 pysam/0.22
```
Or use conda/mamba with the provided environment file:
```
mamba env create -f envs/kinnex_stats.yaml
conda activate kinnex_stats
```

## Setup

### 1. Edit the manifest

Create `manifest.tsv` (tab-separated, no header). Comment lines starting with `#` are ignored.

```
# sample    semicolon-delimited BAM paths
sample1     /path/to/sample1.flnc.bam
sample2     /path/to/sample2_cell1.flnc.bam;/path/to/sample2_cell2.flnc.bam
```

### 2. Edit `config.yaml`

```yaml
manifest: "manifest.tsv"          # path to manifest
stats_script: "kinnex_stats_multicore.py"
out_csv: "results/kinnex_stats.csv"
min_length: 0                      # skip reads shorter than this (0 = no filter)

resources:
  kinnex_stats:
    threads: 3    # set to max BAMs-per-sample in your manifest
    mem: 6        # GB per thread
    hrs: 6
  merge_stats:
    threads: 1
    mem: 4
    hrs: 1
```

Set `threads` under `kinnex_stats` to the maximum number of BAMs any single sample has — this controls how many BAMs are read in parallel within a job.

## Running

### Local

```bash
snakemake --cores 8
```

### Cluster (SGE example)

```bash
snakemake \
    --profile profile/sge \
    --jobs 50
```

The `resources.mem` and `resources.hrs` values in `config.yaml` map to SGE memory and walltime limits. Each sample runs as an independent cluster job, so all samples can run in parallel.

A minimal SGE profile (`profile/sge/config.yaml`) looks like:

```yaml
cluster: >
  qsub -cwd -V
  -pe smp {threads}
  -l h_vmem={resources.mem}G,h_rt={resources.hrs}:00:00
jobs: 50
latency-wait: 60
rerun-incomplete: true
keep-going: true
```

### Dry run

```bash
snakemake -n
```

## Adding samples

Edit `manifest.tsv` to add new rows, then rerun. Snakemake will only process the new samples; existing per-sample CSVs are reused and the final table is regenerated.

## Performance & Benchmarking

The `benchmarking/` folder contains tools to measure runtime on your own data
and choose appropriate `config.yaml` resource settings.

**Quick rules of thumb** (refine with the benchmarking toolkit):

| `config.yaml` parameter | Guidance |
|--------------------------|---------|
| `threads` | Set to the maximum number of BAMs any single sample has — parallelism saturates there |
| `mem` (GB per thread) | ~170–290 MB per million reads ÷ threads, rounded up to next GB, with 20% headroom |
| `hrs` | Measure your largest sample with the benchmark, then multiply by 1.5 |

See [`benchmarking/README.md`](benchmarking/README.md) for full instructions,
a parameter-sweep benchmark script, and an example results table.

## File structure

```
.
├── Snakefile
├── config.yaml
├── manifest.tsv
├── kinnex_stats_multicore.py   # per-sample stats script
├── envs/
│   └── kinnex_stats.yaml       # conda environment
├── benchmarking/
│   ├── subsample_bam.py        # subsample N reads from a real BAM
│   ├── run_benchmark.py        # timing + memory benchmark harness
│   └── README.md               # benchmarking guide and results table
└── results/
    ├── kinnex_stats.csv         # final merged output
    ├── per_sample/
    │   └── {sample}.stats.csv  # per-sample intermediates (kept)
    └── logs/
        ├── {sample}/
        │   └── kinnex_stats.log
        └── merge_stats.log
```
