# =============================================================================
# Snakefile — Kinnex FLNC Stats (multi-BAM-per-sample)
#
# Each sample may have 1–N BAMs (e.g. multiple SMRT cells). All BAMs for a
# sample are streamed and pooled in a single SGE job before stats are computed.
# No merged BAMs or intermediate FASTQ files are written.
#
# Manifest format (TSV, no header):
#   sample_name <TAB> /path/bam1.bam;/path/bam2.bam[;...]
#
# Usage: snakemake -s Snakefile --configfile config.yaml --cores <N>
# See README.md for full documentation and cluster submission instructions.
#
# Author:  KM
# Created: 2025-10
# =============================================================================

import os

configfile: "config.yaml"

# ---------------------------------------------------------------------------
# Wildcard constraints
# ---------------------------------------------------------------------------
wildcard_constraints:
    sample="[^/]+",

# ---------------------------------------------------------------------------
# Manifest parsing  (TSV: sample <tab> bam1;bam2;...)
# ---------------------------------------------------------------------------
SAMPLES = []
BAMS    = {}   # sample -> list of bam paths

with open(config["manifest"]) as fh:
    for lineno, line in enumerate(fh, 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            raise ValueError(
                f"Manifest line {lineno} has {len(fields)} field(s); "
                "expected at least 2 (sample, semicolon-delimited BAM paths)"
            )
        sample = fields[0]
        bams   = [b.strip() for b in fields[1].split(";") if b.strip()]
        if not bams:
            raise ValueError(f"Manifest line {lineno}: no BAM paths found for '{sample}'")
        if sample in BAMS:
            raise ValueError(f"Duplicate sample '{sample}' on line {lineno}")
        SAMPLES.append(sample)
        BAMS[sample] = bams

if not SAMPLES:
    raise ValueError("No samples found in manifest!")

# ---------------------------------------------------------------------------
# Target
# ---------------------------------------------------------------------------
rule all:
    input:
        config["out_csv"]


# ---------------------------------------------------------------------------
# Per-sample stats  (streams all BAMs, no merged intermediates)
# ---------------------------------------------------------------------------
rule kinnex_stats:
    input:
        bams=lambda wc: BAMS[wc.sample],
        script=config["stats_script"],
    output:
        csv="results/per_sample/{sample}.stats.csv",
    log:
        "results/logs/{sample}/kinnex_stats.log",
    params:
        min_length=config.get("min_length", 0),
        bam_args=lambda wc: " ".join(BAMS[wc.sample]),
    threads: config["resources"]["kinnex_stats"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["kinnex_stats"]["mem"] * attempt,
        hrs=config["resources"]["kinnex_stats"]["hrs"],
    conda:
        "envs/kinnex_stats.yaml"
    envmodules:
        "python/3.11",
        "pysam/0.22",
    shell:
        """
        python {input.script} \
            --bams {params.bam_args} \
            --sample {wildcards.sample} \
            --workers {threads} \
            --min-length {params.min_length} \
            --out {output.csv} \
            --progress \
            2> {log}
        """


# ---------------------------------------------------------------------------
# Merge per-sample CSVs into one final table
# ---------------------------------------------------------------------------
rule merge_stats:
    input:
        csvs=expand("results/per_sample/{sample}.stats.csv", sample=SAMPLES),
    output:
        csv=config["out_csv"],
    log:
        "results/logs/merge_stats.log",
    threads: config["resources"]["merge_stats"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["merge_stats"]["mem"] * attempt,
        hrs=config["resources"]["merge_stats"]["hrs"],
    run:
        import csv, os

        fieldnames = ["Sample", "n", "n_bases", "AvgLen", "MedLen", "N50Len", "AvgQ", "MedQ", "N50Q"]
        all_rows = []

        for csv_path in input.csvs:
            with open(csv_path, newline="") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    all_rows.append(row)

        all_rows.sort(key=lambda r: r["Sample"])

        os.makedirs(os.path.dirname(output.csv) or ".", exist_ok=True)
        with open(output.csv, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_rows)

        with open(log[0], "w") as fh:
            fh.write(f"Merged {len(all_rows)} rows from {len(input.csvs)} files.\n")
