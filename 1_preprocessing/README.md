# FastQC — Per-Sample Read QC

Runs FastQC on R1 and R2 of every paired-end sample, in parallel, as one SLURM job per sample.

## Overview

For each sample, a per-sample SLURM script is generated from the `BANANA`-templated `FastQC.sh` and submitted independently. Each job:

1. Locates `${SAMPLE}_R1.fastq.gz` and `${SAMPLE}_R2.fastq.gz` under the input directory.
2. Validates gzip integrity with `gunzip -t`.
3. Runs FastQC on R1 and R2 separately, writing HTML + ZIP reports to the output directory.

The job is **idempotent**: if both `${base}_R1_fastqc.html` and `${base}_R1_fastqc.zip` already exist for a given read, that read is skipped (R1 and R2 are checked independently).

## Adapting to your environment

Paths and SLURM settings in this README reflect the original run on the Harvard FASRC Cannon cluster. To run this elsewhere, replace the following:

| What | This run | Replace with |
|---|---|---|
| Input FASTQ dir | `/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/Ginkgo_rpip` | your input directory |
| Output dir | `/n/netscratch/hhealy_lab/avdarling/fastqc` | your output directory |
| Conda env | `/n/netscratch/hhealy_lab/avdarling_conda_envs/fastqc_env` | a conda env with `fastqc` installed |
| Script dir | `/n/home11/avdarling/scripts/` | wherever you keep your scripts |
| Per-sample script output | `/n/home11/avdarling/scripts/FastQCscripts/` | a directory the generator can write to |
| SLURM log dir | `/n/home11/avdarling/slurm/` | wherever you want SLURM logs |
| SLURM partitions | `shared`, `sapphire` | partitions available on your cluster |

## Requirements

### SLURM resources (per sample)
- 1 CPU core
- 4 GB RAM
- 9 hour wall time
- Partitions: `shared`, `sapphire`

### Conda environment
`/n/netscratch/hhealy_lab/avdarling_conda_envs/fastqc_env`

## Input

Paired-end gzipped FASTQs located under:
`/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/Ginkgo_rpip`

Expected naming convention:
- `${SAMPLE}_R1.fastq.gz`
- `${SAMPLE}_R2.fastq.gz`

`SAMPLE` is hard-coded to `BANANA` in the template and replaced per job by the generator script.

## Output

`/n/netscratch/hhealy_lab/avdarling/fastqc/`

For each sample, four files (where `${base}` = `${SAMPLE}`):
- `${base}_R1_fastqc.html` — interactive QC report for R1
- `${base}_R1_fastqc.zip` — raw FastQC data for R1
- `${base}_R2_fastqc.html` — interactive QC report for R2
- `${base}_R2_fastqc.zip` — raw FastQC data for R2

## Logs

SLURM stdout/stderr land in `/n/home11/avdarling/slurm/`:
- `%j.FastQC.out.output`
- `%j.FastQC.out.err`

(`%j` is the SLURM job ID.)

## Scripts

All scripts live in `/n/home11/avdarling/scripts/`:

- **`FastQC.sh`** — single-sample template with `SAMPLE="BANANA"`. Not run directly; copied + edited per sample.
- **`makeScriptsforFastQC.sh`** — `cd`s into the input directory, finds every `*_R1.fastq.gz`, derives the sample name with `basename`, and uses `sed "s@BANANA@${sample}@g"` to write a per-sample copy of `FastQC.sh` to `/n/home11/avdarling/scripts/FastQCscripts/${sample}.fastqc.sh`.
- **`submitAllScriptsForFastQC.sh`** (or equivalent) — loops over every `*.sh` in `FastQCscripts/` and `sbatch`es it.

### Run order

```
bash /n/home11/avdarling/scripts/makeScriptsforFastQC.sh        # generate per-sample scripts
bash /n/home11/avdarling/scripts/submitAllScriptsForFastQC.sh   # sbatch each one
```

After all jobs finish, the `${OUTDIR}` will contain HTML + ZIP reports for every R1 and R2 across all samples. These can be aggregated downstream with MultiQC if desired.

## Notes & caveats

- Output goes to `netscratch`, which has a retention window — copy the reports to persistent storage if you need them long-term.
- R1 and R2 are skipped independently, so a partial rerun (e.g. R2 finished but R1 failed) will only re-run the missing read.
- The conda env path lives on `netscratch`; if scratch is cleaned, the env will need to be recreated before this script can run.
