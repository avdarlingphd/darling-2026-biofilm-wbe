# RPIP Pathogen Genome Coverage Estimation

Pipeline for estimating reference genome coverage of RPIP panel-targeted pathogens across hospital wastewater metagenomic samples. For each sample × pathogen combination, Kraken2-classified reads are extracted, aligned to an NCBI reference genome, and coverage metrics are calculated.

---

## Configuration

Edit these variables at the top of each script before running.

### Paths — `kraken_coverage_array.sh`

| Variable | Description |
|----------|-------------|
| `COMBINATIONS` | Tab-separated file of sample × taxid pairs (one per line) |
| `KRAKEN` | Kraken2 output file per sample: `{DIR}/{SAMPLE}.kraken` |
| `REPORT` | Kraken2 report file per sample: `{DIR}/{SAMPLE}.kreport` |
| `R1` / `R2` | Paired FASTQ files per sample |
| `REF` | Reference genome per taxid: `{BASE}/taxid_{TAXID}/reference_genome.fna` |
| `OUTDIR` | Output directory per combination: `{BASE}/{SAMPLE}_{TAXID}/` |
| `PYTHON` | Python binary with biopython installed (used for KrakenTools) |
| `BWA` | BWA-MEM binary |
| `SAMTOOLS` | samtools binary |
| `MOSDEPTH` | mosdepth binary |
| `KRAKENTOOLS` | Path to `extract_kraken_reads.py` |

### SLURM settings — `kraken_coverage_array.sh`

| Setting | Current value | Description |
|---------|---------------|-------------|
| `--array` | `0-4540%20` | One task per combination; max 20 running at once |
| `-c` | `8` | CPUs per task |
| `--mem` | `100G` | Memory per task |
| `-t` | `3-00:00` | Wall time (3 days) |
| `-p` | `hsph,shared,sapphire` | SLURM partitions |

### Paths — `download_references.sh`

| Variable | Description |
|----------|-------------|
| `TAXIDS_FILE` | File of taxids to download (one per line) |
| `REF_BASE` | Base directory where reference genomes will be saved |

### Paths — `aggregate_results.sh`

| Variable | Description |
|----------|-------------|
| `RESULTS_BASE` | Base directory containing all per-combination output folders |
| `OUTPUT` | Filename for the aggregated output TSV |

---

## Overview

**Samples:** 209 hospital wastewater samples (Yale New Haven Hospital, 2024)
**Targets:** 101 RPIP panel pathogens
**Combinations:** 4,541 sample × taxid pairs
**Kraken2 settings:** confidence threshold 0.5, minimum hit groups 3

## Pipeline

```
1. download_references.sh     — Download NCBI reference genomes for all taxids
2. kraken_coverage_array.sh   — SLURM array: extract reads → align → calculate coverage
3. aggregate_results.sh       — Combine all per-combination TSVs into one table
```

### Step 1: Download Reference Genomes

Downloads one reference genome per taxid from NCBI using the `datasets` CLI.

```bash
sbatch scripts/download_references.sh
```

- Reads taxids from `data/taxids.txt` (one per line)
- Skips taxids with existing reference genomes
- Saves to: `reference_genomes/taxid_XXXXX/reference_genome.fna`

### Step 2: Run Coverage Array

SLURM array job — one task per sample × taxid combination.

```bash
sbatch scripts/kraken_coverage_array.sh
```

Each task:
1. Extracts Kraken2-classified reads for the target taxid using KrakenTools
2. Aligns extracted reads to the reference genome with BWA-MEM
3. Calculates coverage metrics with mosdepth
4. Writes `coverage_metrics.tsv` to the output directory

If `coverage_metrics.tsv` already exists for a combination, the task exits immediately — safe to resubmit without reprocessing completed work.

Output columns: `sample`, `taxid`, `extracted_r1_reads`, `extracted_r2_reads`, `mapped_reads`, `reference_length`, `breadth_1x`, `breadth_10x`, `mean_depth`, `median_depth`

### Step 3: Aggregate Results

```bash
bash scripts/aggregate_results.sh
```

Collects all per-combination `coverage_metrics.tsv` files into a single table: `all_coverage_metrics.tsv`

---

## Coverage Metrics

| Metric | Description |
|--------|-------------|
| `breadth_1x` | Fraction of reference genome with ≥1x read depth |
| `breadth_10x` | Fraction of reference genome with ≥10x read depth |
| `mean_depth` | Mean reads per base across the full reference genome (including uncovered bases) |
| `median_depth` | Median reads per base, computed from per-base depth file |

> **Note:** For sparse metagenomic data, median depth is often 0. `breadth_1x` is generally more informative for presence/absence determination.

---

## Quick Summary Stats

After the array job completes, compute summary statistics across all results:

```bash
cat /path/to/kraken_coverage/*/coverage_metrics.tsv | grep -v "^sample" | awk '
BEGIN {total=0; breadth1=0; breadth10=0; depth_sum=0; detected=0}
{
    total++
    if ($7 > 0) breadth1++
    if ($8 > 0) breadth10++
    if ($7 > 0) { depth_sum += $9; detected++ }
}
END {
    print "Total combinations: " total
    print ">=1x coverage: " breadth1 " (" int(breadth1/total*100) "%)"
    print ">=10x coverage: " breadth10 " (" int(breadth10/total*100) "%)"
    print "Mean depth (detected only): " depth_sum/detected
}'
```

Percentages reflect the fraction of sample × pathogen combinations with any coverage, not the fraction of genome bases covered.

---

## Dependencies

| Tool | Version | Use |
|------|---------|-----|
| Kraken2 | 2.1.7 | Metagenomic classification (pre-run) |
| KrakenTools | — | Extract taxid-classified reads |
| BWA-MEM | — | Read alignment |
| samtools | — | BAM processing |
| mosdepth | — | Coverage calculation |
| NCBI datasets CLI | — | Reference genome download |

## Conda Environments

```
seqtk_env    — Python / KrakenTools / biopython  (use base miniconda Python if seqtk_env is broken)
bwa_env      — BWA, samtools, mosdepth
kraken2_2.17 — NCBI datasets CLI
```

---

## Directory Structure

```
data/
  taxids.txt                                     # All 101 taxid targets
  new_taxids_to_download.txt                     # Taxids needing new reference download
  kraken_coverage_combinations_hospital.txt      # 4,541 sample × taxid pairs
  nonmunicipal_samples.txt                       # 209 hospital sample IDs

scripts/
  download_references.sh
  kraken_coverage_array.sh
  aggregate_results.sh

reference_genomes/
  taxid_XXXXX/
    reference_genome.fna

kraken_coverage/
  SAMPLE_TAXID/
    coverage_metrics.tsv
    mapped.sorted.bam
    mosdepth.mosdepth.summary.txt
    mosdepth.quantized.bed.gz
    mosdepth.per-base.bed.gz
    bwa.log

all_coverage_metrics.tsv                         # Final aggregated output
```
