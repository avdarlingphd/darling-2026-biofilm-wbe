# Parallel Metagenomic Profiling Pipeline

A SLURM-based bash pipeline that runs multiple taxonomic and antimicrobial-resistance (AMR) profilers on paired-end shotgun metagenomic samples. Per-sample SLURM jobs are generated from a single template script and submitted in parallel.

## Overview

For each sample, the main workflow job (`fullWorkflow.sh`) runs the following tools sequentially within one SLURM job:

1. **Kraken2** — k-mer–based taxonomic classification
2. **CARD (via DIAMOND blastx)** — antimicrobial resistance gene screening
3. **Centrifuger** — taxonomic classification (with quantification step)
4. **Sylph** — containment-based species profiling, with GTDB taxonomy conversion
5. **MetaPhlAn** — marker-gene–based taxonomic profiling

Each step is **idempotent**: if the expected output already exists, the step is skipped. This means the script can be safely re-run after partial failures.

## Adapting to your environment

Paths and SLURM settings in this README reflect the original run on the Harvard FASRC Cannon cluster. To run this pipeline elsewhere, replace the following:

| What | This run | Replace with |
|---|---|---|
| Input FASTQ dir | `/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/Ginkgo_rpip` | your input directory |
| Script dir | `/n/home11/avdarling/scripts/` | wherever you keep your scripts |
| Per-sample script output | `/n/home11/avdarling/scripts/workflow_scripts_per_sample/` | a directory the generator can write to |
| SLURM log dir | `/n/home11/avdarling/slurm/` | wherever you want SLURM logs |
| Kraken2 DB | `/n/holylabs/hhealy_lab/Lab/databases/kraken2_pluspf` | path to a Kraken2 PlusPF (or other) DB |
| CARD DB | `/n/home11/avdarling/databases/CARD/card_db` | a CARD DIAMOND database |
| Centrifuger DB | `/n/holylabs/hhealy_lab/Lab/databases/centrifuger_db/cfr_hpv+gbsarscov2` | a Centrifuger index |
| Sylph DB | `/n/holylabs/hhealy_lab/Lab/databases/gtdb-r220-c200-dbv1.syldb` | a Sylph `.syldb` |
| Kraken2 output dir | `/n/netscratch/hhealy_lab/avdarling/kraken_out/...` | your scratch/output dir |
| Centrifuger output dir | `/n/netscratch/hhealy_lab/avdarling/kraken_out/centrifuger_output` | your scratch/output dir |
| CARD output dir | `/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/CARD_output` | your output dir |
| Sylph output dir | `/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/sylph_output` | your output dir |
| MetaPhlAn output dir | `/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/metaphlan_output_parallel_nonanadama_rerun` | your output dir |
| SLURM partitions | `hsph`, `sapphire`, `shared` | partitions available on your cluster |
| Conda envs | absolute paths under `/n/holylabs/...` and named envs (`diamond_env`, `metaphlan_env`) | your own conda envs with the appropriate tools installed |

Quick sanity checks: confirm your cluster has enough memory (~128 GB) and wall time (~6 hours) per sample, and that the conda envs listed below contain the tools at the versions you want.

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

The script verifies gzip integrity with `gunzip -t` before running. The `SAMPLE` variable is hard-coded to `BANANA` in the template and replaced per job by the `makeScriptsFullWorkflow.sh` generator (see "Parallel submission" below).

## Output

`${base}` below is the sample name (the `_R1.fastq.gz` suffix stripped from the input filename).

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
- `%j.kraken_bracken.output`
- `%j.kraken_bracken.err`

(`%j` is the SLURM job ID.)

## Parallel submission

Per-sample scripts are generated from a `BANANA`-templated base script and submitted in bulk. Helper scripts live in `/n/home11/avdarling/scripts/`:

- **Template:** `fullWorkflow.sh` — the script described above with `SAMPLE="BANANA"`.
- **Generator:** `makeScriptsFullWorkflow.sh` — `cd`s into the input directory, finds every `*_R1.fastq.gz`, derives the sample name, and uses `sed 's@BANANA@<sample>@'` to write a per-sample copy of `fullWorkflow.sh` into `/n/home11/avdarling/scripts/workflow_scripts_per_sample/`.
- **Submitter:** `submitAllScriptsForFullWorkflow.sh` — loops over every `*.sh` in `workflow_scripts_per_sample/` and `sbatch`es it.

### Run order

```
bash /n/home11/avdarling/scripts/makeScriptsFullWorkflow.sh
bash /n/home11/avdarling/scripts/submitAllScriptsForFullWorkflow.sh
```

## Notes & caveats

- All steps depend on the input FASTQs being valid; the early `gunzip -t` check guards against corrupt uploads.
- Outputs are split between `holylabs` (persistent project space) and `netscratch` (high-throughput scratch). Anything on netscratch should be copied off before the scratch retention window expires.
- The CARD step concatenates R1+R2 rather than running them separately; this is fine for translated alignment but means read-pairing information is not preserved in `${base}.card`.
- The SLURM log filename `%j.kraken_bracken.*` is historical — Bracken is no longer part of the pipeline. Consider renaming to something tool-agnostic.
