#!/usr/bin/env python3
"""
run_benchmark.py

Benchmark the kinnex_stats_multicore.py pipeline using real (subsampled) BAM data.

Two subcommands:

  standalone   Subsample a real BAM at multiple sizes and run kinnex_stats_multicore.py
               directly, sweeping over n_bams and workers. Best for per-script performance.

  snakemake    Run the full Snakemake pipeline against a real manifest and measure
               end-to-end wall time including scheduler overhead.

Examples:

  # Quick sanity check (~2 min):
  python run_benchmark.py standalone \\
      --bam /data/sample.flnc.bam \\
      --reads-list 10000 100000 \\
      --bams-list 1 3 \\
      --workers-list 1 4 \\
      --results quick_results.csv

  # Full matrix (~30-60 min depending on hardware):
  python run_benchmark.py standalone --bam /data/sample.flnc.bam --results full_results.csv

  # Snakemake end-to-end:
  python run_benchmark.py snakemake \\
      --manifest /data/manifest.tsv \\
      --cores 8 \\
      --results snakemake_results.csv
"""

import argparse
import csv
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from itertools import product
from pathlib import Path

# Path to kinnex_stats_multicore.py relative to this file
_HERE = Path(__file__).resolve().parent
_DEFAULT_SCRIPT = _HERE.parent / "kinnex_stats_multicore.py"
_DEFAULT_SNAKEFILE = _HERE.parent / "Snakefile"
_DEFAULT_CONFIG = _HERE.parent / "config.yaml"

DEFAULT_READS_LIST = [10_000, 100_000, 1_000_000]
DEFAULT_BAMS_LIST = [1, 3, 6]
DEFAULT_WORKERS_LIST = [1, 2, 4, 8]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _poll_rss_kb(pid: int, interval: float, result: list, stop: threading.Event) -> None:
    """
    Background thread: poll /proc/<pid>/status for VmRSS every `interval` seconds.
    Records the peak value in result[0]. Stops when the process exits or stop is set.

    VmRSS covers the main process RSS. The ProcessPoolExecutor workers are separate
    processes; their memory is not included here, but the dominant cost — the pooled
    all_lengths and all_qualities lists — lives in the main process after workers return.
    """
    peak_kb = 0
    status_path = f"/proc/{pid}/status"
    while not stop.is_set():
        try:
            with open(status_path) as f:
                for line in f:
                    if line.startswith("VmRSS:"):
                        kb = int(line.split()[1])
                        if kb > peak_kb:
                            peak_kb = kb
                        break
        except FileNotFoundError:
            break  # process has exited
        time.sleep(interval)
    result.append(peak_kb)


def _measure_peak_rss_kb(proc: subprocess.Popen, poll_interval: float = 0.1) -> int:
    """Poll a running Popen process for peak RSS; return KB. Falls back to 0 on non-Linux."""
    if platform.system() != "Linux":
        proc.wait()
        return 0
    result: list = []
    stop = threading.Event()
    t = threading.Thread(target=_poll_rss_kb, args=(proc.pid, poll_interval, result, stop),
                         daemon=True)
    t.start()
    proc.wait()
    stop.set()
    t.join()
    return result[0] if result else 0


def compute_throughput(n_reads: int, n_bases: int, wall_time_s: float) -> dict:
    if wall_time_s <= 0:
        return {"reads_per_sec": None, "mbases_per_sec": None}
    return {
        "reads_per_sec": round(n_reads / wall_time_s, 1),
        "mbases_per_sec": round(n_bases / 1e6 / wall_time_s, 3),
    }


def parse_result_csv(csv_path: str) -> dict:
    """Return the first data row from a kinnex_stats output CSV as a dict."""
    try:
        with open(csv_path, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                return row
    except Exception:
        pass
    return {}


def print_summary_table(rows: list) -> None:
    if not rows:
        return
    cols = ["n_reads_subsample", "n_bams", "workers", "wall_time_s",
            "peak_rss_mb", "reads_per_sec", "mbases_per_sec", "status"]
    widths = {c: max(len(c), max(len(str(r.get(c, ""))) for r in rows)) for c in cols}
    header = "  ".join(c.ljust(widths[c]) for c in cols)
    sep = "  ".join("-" * widths[c] for c in cols)
    print(header)
    print(sep)
    for r in rows:
        print("  ".join(str(r.get(c, "")).ljust(widths[c]) for c in cols))


def write_results_csv(rows: list, path: str) -> None:
    if not rows:
        return
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nResults written to: {path}")


# ---------------------------------------------------------------------------
# Standalone mode
# ---------------------------------------------------------------------------

def build_matrix(reads_list, bams_list, workers_list) -> list:
    return [
        {"n_reads": r, "n_bams": b, "workers": w}
        for r, b, w in product(reads_list, bams_list, workers_list)
    ]


def ensure_subsampled_bam(source_bam: str, n_reads: int, out_dir: str, seed: int = 42) -> str:
    """Return path to a subsampled BAM, generating it if needed."""
    sys.path.insert(0, str(_HERE))
    from subsample_bam import subsample
    out_path = os.path.join(out_dir, f"sub_{n_reads}_seed{seed}.bam")
    return subsample(source_bam, n_reads, out_path, seed=seed)


def run_one_standalone(
    script_path: str,
    bam_paths: list,
    sample: str,
    workers: int,
    out_csv: str,
) -> dict:
    """Run kinnex_stats_multicore.py as a subprocess and return timing/memory metrics."""
    cmd = [
        sys.executable, str(script_path),
        "--bams", *bam_paths,
        "--sample", sample,
        "--workers", str(workers),
        "--out", out_csv,
    ]
    t0 = time.perf_counter()
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        peak_rss_kb = _measure_peak_rss_kb(proc)  # waits for proc to finish
        wall_time = time.perf_counter() - t0
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            raise subprocess.CalledProcessError(proc.returncode, cmd, stdout, stderr)
        return {
            "wall_time_s": round(wall_time, 3),
            "peak_rss_kb": peak_rss_kb,
            "peak_rss_mb": round(peak_rss_kb / 1024, 1),
            "stderr": stderr,
            "status": "ok",
        }
    except subprocess.CalledProcessError as e:
        wall_time = time.perf_counter() - t0
        stderr = e.stderr or ""
        return {
            "wall_time_s": round(wall_time, 3),
            "peak_rss_kb": 0,
            "peak_rss_mb": 0.0,
            "stderr": stderr,
            "status": f"error: {stderr[:120].strip()}",
        }


def run_standalone(args) -> None:
    script_path = args.script or _DEFAULT_SCRIPT
    if not os.path.exists(script_path):
        sys.exit(f"Script not found: {script_path}")
    if not os.path.exists(args.bam):
        sys.exit(f"BAM not found: {args.bam}")

    out_dir = args.out_dir or tempfile.mkdtemp(prefix="bench_bams_")
    bam_cache_dir = os.path.join(out_dir, "bams")
    run_dir = os.path.join(out_dir, "runs")
    os.makedirs(bam_cache_dir, exist_ok=True)
    os.makedirs(run_dir, exist_ok=True)

    matrix = build_matrix(args.reads_list, args.bams_list, args.workers_list)
    total = len(matrix) * args.repeats
    print(f"Standalone benchmark: {total} run(s) across {len(matrix)} parameter combination(s)")
    benchmark_start = time.perf_counter()
    if args.dry_run:
        for m in matrix:
            print(f"  reads={m['n_reads']:>9,}  n_bams={m['n_bams']}  workers={m['workers']}")
        return

    rows = []
    run_idx = 0
    for m in matrix:
        n_reads, n_bams, workers = m["n_reads"], m["n_bams"], m["workers"]

        # Subsample the source BAM once per (n_reads, seed) combination.
        # For n_bams > 1 use different seeds to simulate distinct SMRT cells.
        print(f"\n[{run_idx+1}/{total}] reads={n_reads:,}  n_bams={n_bams}  workers={workers}")
        bam_paths = []
        for i in range(n_bams):
            bam_paths.append(
                ensure_subsampled_bam(args.bam, n_reads, bam_cache_dir, seed=42 + i)
            )

        for rep in range(1, args.repeats + 1):
            run_idx += 1
            out_csv = os.path.join(run_dir, f"r{n_reads}_b{n_bams}_w{workers}_rep{rep}.csv")
            metrics = run_one_standalone(
                script_path, bam_paths,
                sample=f"bench_r{n_reads}_b{n_bams}",
                workers=workers,
                out_csv=out_csv,
            )
            result_row = parse_result_csv(out_csv)
            n_reported = int(result_row.get("n", 0)) if result_row else 0
            n_bases = int(result_row.get("n_bases", 0)) if result_row else 0
            total_reads = n_reads * n_bams

            throughput = compute_throughput(n_reported, n_bases, metrics["wall_time_s"])
            row = {
                "mode": "standalone",
                "n_reads_subsample": n_reads,
                "n_bams": n_bams,
                "total_reads": total_reads,
                "workers": workers,
                "repeat": rep,
                "wall_time_s": metrics["wall_time_s"],
                "peak_rss_kb": metrics["peak_rss_kb"],
                "peak_rss_mb": metrics["peak_rss_mb"],
                "reads_per_sec": throughput["reads_per_sec"],
                "mbases_per_sec": throughput["mbases_per_sec"],
                "n_reported": n_reported,
                "status": metrics["status"],
            }
            rows.append(row)
            print(
                f"  rep {rep}: {metrics['wall_time_s']:.2f}s  "
                f"{metrics['peak_rss_mb']:.0f} MB RSS  "
                f"{throughput['reads_per_sec']:,.0f} reads/s  "
                f"[{metrics['status']}]"
            )
            if n_reported != total_reads and metrics["status"] == "ok":
                print(
                    f"  Warning: expected {total_reads:,} reads, "
                    f"got {n_reported:,} in output CSV"
                )

    if not args.keep_bams and args.out_dir is None:
        shutil.rmtree(bam_cache_dir, ignore_errors=True)

    total_elapsed = time.perf_counter() - benchmark_start
    print(f"\n--- Summary (total elapsed: {total_elapsed:.1f}s) ---")
    print_summary_table(rows)
    write_results_csv(rows, args.results)


# ---------------------------------------------------------------------------
# Snakemake mode
# ---------------------------------------------------------------------------

def _parse_snakemake_log(log_path: str) -> dict:
    """
    Parse a Snakemake log file and return {sample: wall_time_s} for kinnex_stats jobs.

    Snakemake logs timestamps before each rule block and before each "Finished job N."
    line. We match job start (rule kinnex_stats + wildcards: sample=NAME) to finish
    (Finished job N.) via the jobid field.
    """
    import re
    from datetime import datetime

    # Snakemake timestamp format: [Mon Jan  1 00:00:00 2026]
    ts_re = re.compile(r'^\[(\w{3} \w{3}\s+\d+ \d+:\d+:\d+ \d{4})\]')
    finished_re = re.compile(r'^Finished job (\d+)\.')

    job_starts: dict = {}   # jobid -> (sample, datetime)
    per_sample: dict = {}   # sample -> wall_time_s

    current_ts = None
    in_kinnex = False
    cur_jobid = None
    cur_sample = None

    with open(log_path) as f:
        for raw in f:
            line = raw.strip()
            m = ts_re.match(line)
            if m:
                current_ts = datetime.strptime(m.group(1), "%a %b %d %H:%M:%S %Y")
                in_kinnex = False
                cur_jobid = None
                cur_sample = None
                continue

            if line == "rule kinnex_stats:":
                in_kinnex = True
                continue

            if in_kinnex:
                if line.startswith("jobid:"):
                    cur_jobid = int(line.split(":", 1)[1].strip())
                elif line.startswith("wildcards:") and "sample=" in line:
                    # wildcards may have multiple fields: sample=X, other=Y
                    for part in line.split("wildcards:", 1)[1].split(","):
                        part = part.strip()
                        if part.startswith("sample="):
                            cur_sample = part.split("=", 1)[1].strip()
                    if cur_jobid is not None and cur_sample is not None:
                        job_starts[cur_jobid] = (cur_sample, current_ts)
                continue

            m2 = finished_re.match(line)
            if m2 and current_ts is not None:
                jobid = int(m2.group(1))
                if jobid in job_starts:
                    sample, start = job_starts.pop(jobid)
                    per_sample[sample] = round((current_ts - start).total_seconds(), 1)

    return per_sample


def _find_snakemake_log(work_dir: str) -> str | None:
    """Return path to the most recent Snakemake log file in work_dir."""
    log_dir = os.path.join(work_dir, ".snakemake", "log")
    if not os.path.isdir(log_dir):
        return None
    logs = sorted(
        (os.path.join(log_dir, f) for f in os.listdir(log_dir) if f.endswith(".snakemake.log")),
        key=os.path.getmtime,
    )
    return logs[-1] if logs else None


def run_snakemake_benchmark(args) -> None:
    if not os.path.exists(args.manifest):
        sys.exit(f"Manifest not found: {args.manifest}")
    if not os.path.exists(_DEFAULT_SNAKEFILE):
        sys.exit(f"Snakefile not found: {_DEFAULT_SNAKEFILE}")

    work_dir = tempfile.mkdtemp(prefix="bench_snakemake_")
    print(f"Snakemake benchmark working directory: {work_dir}")

    shutil.copy(_DEFAULT_SNAKEFILE, work_dir)
    shutil.copy(_DEFAULT_CONFIG, work_dir)
    shutil.copy(_DEFAULT_SCRIPT, work_dir)
    shutil.copy(args.manifest, os.path.join(work_dir, "manifest.tsv"))
    envs_src = _HERE.parent / "envs"
    if envs_src.exists():
        shutil.copytree(str(envs_src), os.path.join(work_dir, "envs"))

    cmd = [
        "snakemake",
        "--cores", str(args.cores),
        "--config", "manifest=manifest.tsv",
        "--nolock",
    ]
    if args.use_conda:
        cmd.append("--use-conda")

    print(f"Running: {' '.join(cmd)}")
    t0 = time.perf_counter()
    try:
        result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)
        total_wall = time.perf_counter() - t0
        success = result.returncode == 0
    except FileNotFoundError:
        sys.exit("snakemake not found. Install it or activate the project conda environment.")

    overall_status = "ok" if success else f"error (rc={result.returncode})"
    print(f"\nSnakemake finished in {total_wall:.1f}s  [{overall_status}]")
    if not success:
        print("--- stderr (last 40 lines) ---")
        for line in result.stderr.splitlines()[-40:]:
            print(" ", line)

    # --- Per-sample timing from Snakemake log ---
    per_sample_times: dict = {}
    log_path = _find_snakemake_log(work_dir)
    if log_path:
        try:
            per_sample_times = _parse_snakemake_log(log_path)
        except Exception as e:
            print(f"  Warning: could not parse Snakemake log for per-sample times: {e}")

    # --- Per-sample read stats from output CSVs ---
    per_sample_dir = os.path.join(work_dir, "results", "per_sample")
    rows = []
    if os.path.isdir(per_sample_dir):
        for fname in sorted(os.listdir(per_sample_dir)):
            if not fname.endswith(".stats.csv"):
                continue
            sample = fname.replace(".stats.csv", "")
            stats = parse_result_csv(os.path.join(per_sample_dir, fname))
            n_reads = int(stats.get("n", 0)) if stats else 0
            n_bases = int(stats.get("n_bases", 0)) if stats else 0
            wall_time = per_sample_times.get(sample)
            throughput = compute_throughput(n_reads, n_bases, wall_time or 0)
            rows.append({
                "mode": "snakemake",
                "sample": sample,
                "n_reads": n_reads,
                "n_bases": n_bases,
                "wall_time_s": wall_time if wall_time is not None else "NA",
                "reads_per_sec": throughput["reads_per_sec"] if wall_time else "NA",
                "mbases_per_sec": throughput["mbases_per_sec"] if wall_time else "NA",
                "status": "ok" if success else overall_status,
            })

    # Summary row
    n_samples = len(rows)
    rows.append({
        "mode": "snakemake_total",
        "sample": f"{n_samples} samples",
        "n_reads": sum(int(r["n_reads"]) for r in rows if r["n_reads"] != "NA"),
        "n_bases": sum(int(r["n_bases"]) for r in rows if r["n_bases"] != "NA"),
        "wall_time_s": round(total_wall, 3),
        "reads_per_sec": "",
        "mbases_per_sec": "",
        "status": overall_status,
    })

    # Print summary table
    if rows:
        cols = ["sample", "n_reads", "wall_time_s", "reads_per_sec", "mbases_per_sec", "status"]
        widths = {c: max(len(c), max(len(str(r.get(c, ""))) for r in rows)) for c in cols}
        print("  ".join(c.ljust(widths[c]) for c in cols))
        print("  ".join("-" * widths[c] for c in cols))
        for r in rows:
            print("  ".join(str(r.get(c, "")).ljust(widths[c]) for c in cols))

    write_results_csv(rows, args.results)

    if not args.keep_work_dir:
        shutil.rmtree(work_dir, ignore_errors=True)
    else:
        print(f"Work directory retained at: {work_dir}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark kinnex_stats pipeline using real (subsampled) BAM data.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    sub = parser.add_subparsers(dest="subcommand", required=True)

    # ---- standalone subcommand ----
    sa = sub.add_parser("standalone", help="Benchmark kinnex_stats_multicore.py directly")
    sa.add_argument("--bam", required=True, help="Source BAM to subsample from")
    sa.add_argument(
        "--reads-list", type=int, nargs="+", default=DEFAULT_READS_LIST,
        metavar="N", dest="reads_list",
        help=f"Subsample sizes (default: {DEFAULT_READS_LIST})",
    )
    sa.add_argument(
        "--bams-list", type=int, nargs="+", default=DEFAULT_BAMS_LIST,
        metavar="N", dest="bams_list",
        help=f"Number of BAMs per sample (default: {DEFAULT_BAMS_LIST})",
    )
    sa.add_argument(
        "--workers-list", type=int, nargs="+", default=DEFAULT_WORKERS_LIST,
        metavar="N", dest="workers_list",
        help=f"Worker counts to test (default: {DEFAULT_WORKERS_LIST})",
    )
    sa.add_argument("--repeats", type=int, default=1, help="Repeats per combination (default: 1)")
    sa.add_argument("--results", default="benchmark_results.csv", help="Output CSV path")
    sa.add_argument("--out-dir", default=None, dest="out_dir",
                    help="Scratch directory for BAMs and per-run CSVs (default: system tmp)")
    sa.add_argument("--keep-bams", action="store_true", dest="keep_bams",
                    help="Don't delete subsampled BAMs after benchmark")
    sa.add_argument("--script", default=None,
                    help="Path to kinnex_stats_multicore.py (default: ../kinnex_stats_multicore.py)")
    sa.add_argument("--dry-run", action="store_true", dest="dry_run",
                    help="Print parameter matrix without running")

    # ---- snakemake subcommand ----
    sm = sub.add_parser("snakemake", help="Benchmark the full Snakemake pipeline end-to-end")
    sm.add_argument("--manifest", required=True, help="Manifest TSV pointing to real BAMs")
    sm.add_argument("--cores", type=int, default=8, help="--cores for snakemake (default: 8)")
    sm.add_argument("--results", default="snakemake_results.csv", help="Output CSV path")
    sm.add_argument("--use-conda", action="store_true", dest="use_conda",
                    help="Pass --use-conda to snakemake")
    sm.add_argument("--keep-work-dir", action="store_true", dest="keep_work_dir",
                    help="Don't delete the temporary Snakemake working directory")

    return parser.parse_args()


def main():
    args = parse_args()
    if args.subcommand == "standalone":
        run_standalone(args)
    else:
        run_snakemake_benchmark(args)


if __name__ == "__main__":
    main()
