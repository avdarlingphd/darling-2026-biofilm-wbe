# Parallel Metagenomic Profiling Pipeline

A SLURM-based bash pipeline that runs multiple taxonomic and antimicrobial-resistance (AMR) profilers in series on a single paired-end shotgun metagenomic sample. Designed to be submitted in parallel (one job per sample) by templating the `SAMPLE` variable.

## Overview

For each sample, the script runs the following tools sequentially within a single SLURM job:

1. **Kraken2 + Bracken** — k-mer–based taxonomic classification and abundance re-estimation
2. **CARD (via DIAMOND blastx)** — antimicrobial resistance gene screening
3. **Centrifuger** — taxonomic classification (with quantification step)
4. **Sylph** — containment-based species profiling, with GTDB taxonomy conversion
5. **MetaPhlAn** — marker-gene–based taxonomic profiling

Each step is **idempotent**: if the expected output already exists, the step is skipped. This means the script can be safely re-run after partial failures.

## Requirements

### SLURM resources (per sample)
- 16 CPU cores
- 128 GB RAM (Kraken2 PlusPF DB requires ~74 GB)
- 6 hour wall time
- Partitions: `hsph`, `sapphire`, `shared`

### Conda environments
| Environment | Used for |
|---|---|
| `/n/holylabs/hhealy_lab/Lab/avdarling_conda_envs/kraken2_latest` | Kraken2 |
| `bracken_env` | Bracken |
| `diamond_env` | CARD (DIAMOND blastx) |
| `/n/holylabs/hhealy_lab/Lab/avdarling_conda_envs/centrifuger_env` | Centrifuger |
| `/n/holylabs/hhealy_lab/Lab/avdarling_conda_envs/sylph_env` | Sylph + sylph-tax |
| `metaphlan_env` | MetaPhlAn |

### Reference databases
- Kraken2/Bracken: `/n/holylabs/hhealy_lab/Lab/databases/kraken2_pluspf`
- CARD (DIAMOND): `/n/home11/avdarling/databases/CARD/card_db`
- Centrifuger: `/n/holylabs/hhealy_lab/Lab/databases/centrifuger_db/cfr_hpv+gbsarscov2`
- Sylph: `/n/holylabs/hhealy_lab/Lab/databases/gtdb-r220-c200-dbv1.syldb` (taxonomy: `GTDB_r220`)
- MetaPhlAn: bundled with the conda env (`metaphlan_databases`)

## Input

Paired-end gzipped FASTQs located under `/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/Ginkgo_rpip`:
- `${SAMPLE}_R1.fastq.gz`
- `${SAMPLE}_R2.fastq.gz`

The script verifies gzip integrity with `gunzip -t` before proceeding. The `SAMPLE` variable is replaced per job (here templated as `BANANA`) when submitting in parallel.

## Output

`${base}` below is the sample name with `_R1.fastq.gz` stripped.

### Kraken2
`/n/netscratch/hhealy_lab/avdarling/kraken_out/kraken_output_ct0_5_min_hit_3/`
- `${base}.kraken2` — per-read classifications
- `${base}.k2report` — Kraken2 report

Run with `--confidence 0.5` and `--minimum-hit-groups 3`.

### Bracken
`/n/netscratch/hhealy_lab/avdarling/kraken_out/kraken_output_ct0_5_min_hit_3/bracken_output_ct0_5_min_hit_3/`
- `${base}.bracken_S` — species
- `${base}.bracken_G` — genus
- `${base}.bracken_F` — family
- `${base}.bracken_O` — order
- `${base}.bracken_C` — class
- `${base}.bracken_P` — phylum

Run with read length `-r 150` and threshold `-t 10`.

### CARD (DIAMOND blastx)
`/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/CARD_output/`
- `${base}.card` — DIAMOND tab-separated hits

R1 and R2 are concatenated to a temporary combined FASTQ, aligned with `--max-target-seqs 1 --evalue 1e-10 --id 80`, and the temp file is removed afterward.

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

SLURM stdout/stderr land in:
- `/n/home11/avdarling/slurm/%j.kraken_bracken.output`
- `/n/home11/avdarling/slurm/%j.kraken_bracken.err`

(`%j` is the SLURM job ID. The log filename is shared across all steps despite the name; consider renaming to something tool-agnostic if helpful.)

## Running in parallel

The script template hard-codes `SAMPLE="BANANA"`; the parallel submission wrapper substitutes the actual sample name into a per-sample copy of the script before `sbatch`-ing it. Each sample gets its own job, its own log files, and writes to disjoint output filenames keyed on `${base}`.

## Notes & caveats

- All steps after the first depend on the input FASTQs being valid; the early `gunzip -t` check guards against corrupt uploads.
- Outputs are split between `holylabs` (persistent project space) and `netscratch` (high-throughput scratch). Anything on netscratch should be copied off before the scratch retention window expires.
- The CARD step concatenates R1+R2 rather than running them separately; this is fine for translated alignment but means read-pairing information is not preserved in `${base}.card`.
- Kraken2 confidence and min-hit-groups settings (`0.5` / `3`) are encoded in the output directory name (`ct0_5_min_hit_3`) so multiple parameter sweeps can coexist.
