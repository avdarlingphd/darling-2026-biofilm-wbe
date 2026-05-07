# GeCoCheck Analysis: YNHH Wastewater RPIP 2024

## Overview

Genome coverage analysis of wastewater samples collected as part of the Yale New Haven Hospital (YNHH) respiratory pathogen identification project (RPIP) 2024. Sequencing was performed using the Ginkgo RPIP panel. This pipeline uses GeCoCheck (coverage_pipeline.py) to validate pathogen detection by computing bowtie2-based genome coverage metrics for taxa identified by Kraken2.

---

## Samples

- **Total sequenced**: 250 samples
- **Non-municipal wastewater** (used for coverage analysis): 209 samples
- **Municipal wastewater** (excluded): 41 samples
- Sample IDs follow the format `GKWWBAA#####`
- Raw paired-end FASTQs (R1/R2) located in: `Ginkgo_rpip_fastqs/`
- Sample metadata: `Ginkgo_rpip_fastqs/sample_metadata_nonmunicipal.csv`

---

## Pipeline

### Kraken2 Classification
Reads were classified using Kraken2 with the following parameters:
- Confidence threshold: 0.5 (`ct0_5`)
- Minimum hit groups: 3 (`min_hit_3`)
- Output directory: `kraken_output_ct0_5_min_hit_3/`

### Genome Coverage (GeCoCheck v10)
GeCoCheck was run on the 209 non-municipal samples using Kraken2 kreport files as input for taxon nomination. Key parameters:

```bash
coverage_pipeline.py \
    --processors 8 \
    --sample_metadata sample_metadata_nonmunicipal.csv \
    --project_name Ginkgo_rpip \
    --fastq_dir gcc_merged_fastq/ \
    --kraken_kreport_dir kraken_output_ct0_5_min_hit_3/ \
    --kraken_outraw_dir kraken_output_ct0_5_min_hit_3/ \
    --output_dir Ginkgo_rpip_coverage_results_v10/ \
    --genome_dir gcc_genome_cache/ \
    --bowtie2_db_dir gcc_bowtie2_cache/ \
    --coverage_program Bowtie2 \
    --no_grouped_samples \
    --read_lim 50000 \
    --skip_cleanup
```

FASTQ preparation: Paired-end gzipped FASTQs were merged per sample prior to running GeCoCheck:
```bash
zcat ${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz > ${SAMPLE}.fastq
```

SLURM resources: 8 cores, 100GB RAM, 7-day wall time, `intermediate` partition.

---

## Key Outputs

| File | Description |
|------|-------------|
| `coverage_checker_output.tsv` | Per-sample, per-taxon coverage metrics (mean depth, genome fraction, bowtie2 stats) |
| `no_genome.txt` | Taxa that failed reference genome download |
| `checkpoint.txt` | Pipeline checkpoint tracking |
| `assemblies_using.csv` | Reference genome accessions used per taxon |
| `gcc_genome_cache/` | Downloaded reference genomes (60 genomes, ~253MB) |
| `gcc_bowtie2_cache/` | Bowtie2 indices (360 files, ~834MB) |

---

## Results Summary

### Coverage-Validated Detections
GeCoCheck processed **~60 taxa** that met the read threshold across 209 non-municipal samples. Coverage metrics reported include mean bowtie2 depth, genome fraction, and mean mapping quality per taxon per sample.

### Target Pathogens Not Processed by GeCoCheck
Several target pathogens of interest were detected by Bracken but were absent from the GeCoCheck output. This is expected behavior: GeCoCheck nominates taxa for genome download based on **direct Kraken2 kreport read counts**, not Bracken-redistributed estimates. The following pathogens had insufficient direct Kraken2 reads to clear the pipeline's read threshold:

| Pathogen | Median Direct Kraken Reads | Bracken Inflation Factor |
|---|---|---|
| *Stenotrophomonas maltophilia* | 115 | ~31× |
| *Aeromonas hydrophila* | 1,315 | ~34× |
| *Citrobacter freundii* | 647 | ~37× |
| *Vibrio cholerae* | 14 | ~263× |
| *Naegleria fowleri* | 24 | ~184× |
| *Mycobacterium avium* | 29 | ~2003× |
| *Mycobacterium tuberculosis* | 43 | ~521× |
| *Burkholderia mallei* | 32 | ~176× |

Bracken inflates read counts by redistributing genus-level reads probabilistically to the species level. These inflated estimates do not reflect direct sequencing evidence and cannot substitute for the kreport-based threshold used by GeCoCheck.

---

## Scripts

| Script | Purpose |
|--------|---------|
| `merge_fastq_nonmunicipal.sh` | SLURM array job to merge R1/R2 fastqs for 209 non-municipal samples |
| `run_gecocheck_v10.sh` | SLURM job to run GeCoCheck on 209 non-municipal samples |
| `move_gcc_files.sh` | Migrates genome/bowtie2 caches and key outputs from scratch to lab storage |
| `summarize_coverage.R` | R script to compute per-pathogen mean/median depth and genome fraction from output TSV |

---

## Notes

- GeCoCheck was previously run on all 250 samples (v8) using `ct0_5` kreports; v10 restricts to non-municipal samples and uses the more stringent `ct0_5_min_hit_3` kreports.
- Genome and bowtie2 caches from v8 were retained and reused in v10 to avoid redundant downloads.
- Pipeline run on Harvard FASRC cluster (netscratch for intermediate files, holylabs for long-term storage).
