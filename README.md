# Parallel Metagenomic Profiling Pipeline

A SLURM-based bash pipeline that runs read QC plus multiple taxonomic and antimicrobial-resistance (AMR) profilers on paired-end shotgun metagenomic samples. Per-sample SLURM jobs are generated from a single template script and submitted in parallel.

## Overview

Two SLURM jobs are run per sample:

**FastQC job** (`FastQC.sh`): runs FastQC on R1 and R2 independently.

**Main workflow job** (`fullWorkflow.sh`): runs the following tools sequentially within one SLURM job:

1. **Kraken2** — k-mer–based taxonomic classification
2. **CARD (via DIAMOND blastx)** — antimicrobial resistance gene screening
3. **Centrifuger** — taxonomic classification (with quantification step)
4. **Sylph** — containment-based species profiling, with GTDB taxonomy conversion
5. **MetaPhlAn** — marker-gene–based taxonomic profiling

Each step is **idempotent**: if the expected output already exists, the step is skipped. This means the script can be safely re-run after partial failures.

## Requirements

### SLURM resources

Main workflow (per sample):
- 16 CPU cores
- 128 GB RAM (Kraken2 PlusPF DB requires ~74 GB)
- 6 hour wall time
- Partitions: `hsph`, `sapphire`, `shared`

FastQC (per sample):
- 1 CPU, 4 GB RAM, 9 hour wall time
- Partitions: `shared`, `sapphire`

### Conda environments
| Environment | Used for |
|---|---|
| `/n/netscratch/hhealy_lab/avdarling_conda_envs/fastqc_env` | FastQC |
| `/n/holylabs/hhealy_lab/Lab/avdarling_conda_envs/kraken2_latest` | Kraken2 |
| `diamond_env` | CARD (DIAMOND blastx) |
| `/n/holylabs/hhealy_lab/Lab/avdarling_conda_envs/centrifuger_env` | Centrifuger |
| `/n/holylabs/hhealy_lab/Lab/avdarling_conda_envs/sylph_env` | Sylph + sylph-tax |
| `metaphlan_env` | MetaPhlAn |

### Reference databases
- Kraken2: `/n/holylabs/hhealy_lab/Lab/databases/kraken2_pluspf`
- CARD (DIAMOND): `/n/home11/avdarling/databases/CARD/card_db`
- Centrifuger: `/n/holylabs/hhealy_lab/Lab/databases/centrifuger_db/cfr_hpv+gbsarscov2`
- Sylph: `/n/holylabs/hhealy_lab/Lab/databases/gtdb-r220-c200-dbv1.syldb` (taxonomy: `GTDB_r220`)
- MetaPhlAn: bundled with the conda env (`metaphlan_databases`)

## Input

Paired-end gzipped FASTQs located under `/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/Ginkgo_rpip`:
- `${SAMPLE}_R1.fastq.gz`
- `${SAMPLE}_R2.fastq.gz`

Both the FastQC and main workflow scripts verify gzip integrity with `gunzip -t` before running. The `SAMPLE` variable is hard-coded to `BANANA` in the template and replaced per job by the `makeScripts*` generator (see "Parallel submission" below).

## Output

`${base}` below is the sample name (the `_R1.fastq.gz` suffix stripped from the input filename).

### FastQC
`/n/netscratch/hhealy_lab/avdarling/fastqc/`
- `${base}_R1_fastqc.html`, `${base}_R1_fastqc.zip`
- `${base}_R2_fastqc.html`, `${base}_R2_fastqc.zip`

### Kraken2
`/n/netscratch/hhealy_lab/avdarling/kraken_out/kraken_output_ct0_5_min_hit_3/`
- `${base}.kraken2` — per-read classifications
- `${base}.k2report` — Kraken2 report

Run with `--confidence 0.5` and `--minimum-hit-groups 3` (these settings are encoded in the directory name so multiple parameter sweeps can coexist).

### CARD (DIAMOND blastx)
`/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/CARD_output/`
- `${base}.card` — DIAMOND tab-separated hits

R1 and R2 are concatenated into a temporary combined FASTQ, aligned with `--max-target-seqs 1 --evalue 1e-10 --id 80`, and the temp file is removed afterward.

### Centrifuger
`/n/netscratch/hhealy_lab/avdarling/kraken_out/centrifuger_output/`
- `${base}.centrifuger` — raw per-read assignments
- `${base}.centrifuger_quant.tsv` — quantification (`--output-format 1`)

### Sylph
`/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/sylph_output/`
- `${base}_results.tsv` — Sylph profile (containment ANI, abundance)

`/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/sylph_output/sylph_taxonomy_out/`
- `${base}*.sylphmpa` — MetaPhlAn-style taxonomy file from `sylph-tax taxprof` against GTDB r220

### MetaPhlAn
`/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/metaphlan_output_parallel_nonanadama_rerun/`
- `profiles/${base}_profile.txt` — taxonomic profile (all levels, `--tax_lev a`)
- `bowtie2/${base}.bowtie2.bz2` — read-to-marker mapping
- `sams/${base}.sam.bz2` — SAM alignment

## Logs

SLURM stdout/stderr land in `/n/home11/avdarling/slurm/`:
- Main workflow: `%j.kraken_bracken.output` / `%j.kraken_bracken.err`
- FastQC: `%j.FastQC.out.output` / `%j.FastQC.out.err`

(`%j` is the SLURM job ID.)

## Parallel submission

Per-sample scripts are generated from `BANANA`-templated base scripts and submitted in bulk. Two pairs of helper scripts live in `/n/home11/avdarling/scripts/`:

### Main workflow

- **Template:** `fullWorkflow.sh` — the script described above with `SAMPLE="BANANA"`.
- **Generator:** `makeScriptsFullWorkflow.sh` — `cd`s into the input directory, finds every `*_R1.fastq.gz`, derives the sample name, and uses `sed 's@BANANA@<sample>@'` to write a per-sample copy of `fullWorkflow.sh` into `/n/home11/avdarling/scripts/workflow_scripts_per_sample/`.
- **Submitter:** `submitAllScriptsForFullWorkflow.sh` — loops over every `*.sh` in `workflow_scripts_per_sample/` and `sbatch`es it.

### FastQC

- **Template:** `FastQC.sh` — single-sample FastQC script with `SAMPLE="BANANA"`.
- **Generator:** `makeScriptsforFastQC.sh` — same pattern; finds `*_R1.fastq.gz`, `sed`-substitutes the sample name, and writes per-sample scripts to `/n/home11/avdarling/scripts/FastQCscripts/${sample}.fastqc.sh`.
- **Submitter:** submitted analogously (`for f in *.sh; do sbatch "$f"; done` from inside `FastQCscripts/`).

### Typical run order

```
bash makeScriptsforFastQC.sh           # generate per-sample FastQC scripts
bash submitAllScriptsForFastQC.sh      # (or equivalent) sbatch them
bash makeScriptsFullWorkflow.sh        # generate per-sample workflow scripts
bash submitAllScriptsForFullWorkflow.sh
```

## Notes & caveats

- All steps depend on the input FASTQs being valid; the early `gunzip -t` check in both templates guards against corrupt uploads.
- Outputs are split between `holylabs` (persistent project space) and `netscratch` (high-throughput scratch). Anything on netscratch should be copied off before the scratch retention window expires.
- The CARD step concatenates R1+R2 rather than running them separately; this is fine for translated alignment but means read-pairing information is not preserved in `${base}.card`.
- The main workflow's SLURM log is named `%j.kraken_bracken.*` for historical reasons even though Bracken is no longer part of the pipeline. Consider renaming to something tool-agnostic.
