# BLASTN Validation Pipeline for Kraken2 Metagenomics Classifications

This pipeline validates Kraken2 taxonomic classifications from paired-end metagenomic wastewater sequencing data by extracting Kraken2-classified reads and BLASTing them against the NCBI nucleotide (nt) database. Each sample/pathogen combination is classified as a true positive, false positive, uncertain, or uncultured-dominant hit.

## Overview

Kraken2 is a fast k-mer-based classifier but can produce false positives, particularly for low-abundance organisms or organisms with shared k-mer content. This pipeline provides orthogonal validation by taking reads that Kraken2 assigned to a target taxon and confirming their identity using BLASTN alignment against the full NCBI nt database.

## Pipeline Steps

### Step 1: Read Extraction
Reads classified by Kraken2 to each target pathogen taxon are extracted from paired-end FASTQ files using KrakenTools (`extract_kraken_reads.py`). The `--include-children` flag ensures reads assigned to any descendant taxon are captured. Extracted reads are converted from FASTQ to FASTA using seqtk.

**Scripts:**
- `extract_kraken_base.sh` — SLURM template for per-sample extraction (one job per sample, loops over all target taxids)
- `makeScripts_extract.py` — generates per-sample SLURM scripts and a `submitAll_extract.sh` from the template

```bash
python makeScripts_extract.py \
    --taxid_list sample_taxid_list_all_pathogens.tsv \
    --base_script extract_kraken_base.sh \
    --outdir extract_scripts/

bash extract_scripts/submitAll_extract.sh
```

### Step 2: BLASTN Alignment
Extracted reads are aligned against the NCBI nt database (v5 format, requires BLAST+ ≥ 2.10.0) using `blastn`. To keep runtime manageable, each read set is subsampled to a maximum of 500 reads before BLASTing. Up to 20 target sequences are returned per read (`-max_target_seqs 20`). Results are written as tab-separated files with columns: query ID, scientific name, % identity, alignment length, e-value, and sequence title.

Jobs gracefully skip if FASTA files are missing or BLAST results already exist, making it safe to resubmit the full job list after partial completion.

**Scripts:**
- `blast_validation_base.sh` — SLURM template for per-sample/taxid BLASTN jobs (one job per combination)
- `makeScripts_blast.py` — generates per-combination SLURM scripts and a `submitAll_blast.sh` from the template

```bash
python makeScripts_blast.py \
    --taxid_list sample_taxid_list_all_pathogens.tsv \
    --base_script blast_validation_base.sh \
    --outdir blast_scripts/

bash blast_scripts/submitAll_blast.sh
```

### Step 3: False Positive Detection
BLAST results are parsed and each sample/taxid combination is classified using `false_positive_detection.py`. For each read, all 20 BLAST hits are considered. The number of reads with any hit to each organism is counted and organisms are ranked by read count. The expected organism's rank among all detected organisms, the percentage of reads matching the expected organism, and the prevalence of uncultured/environmental hits are used to classify each combination.

```bash
python false_positive_detection.py \
    --blast_dir /path/to/blast_validation/ \
    --taxid_list sample_taxid_list_all_pathogens.tsv \
    --taxid_names taxid_to_name_all_pathogens.csv \
    --output blast_false_positive_report.tsv
```

## Classification Criteria

| Classification | Criteria |
|---|---|
| **TRUE_POSITIVE** | ≥80% of reads have any BLAST hit to the expected organism across all 20 hits; expected organism is ranked #1 or #2 by read count in both R1 and R2; and >2 reads match the expected organism in both R1 and R2 |
| **FALSE_POSITIVE** | <10% of reads match the expected organism |
| **UNCULTURED_DOMINANT** | >50% of best hits per read are uncultured/environmental sequences — takes priority over all other criteria |
| **UNCERTAIN** | All other cases — including high pct_match but low rank, or 10–79% match |
| **NO_DATA** | BLAST results exist but are empty |
| **BLAST_NOT_RUN** | No output directory found for this combination |

**Organism matching** uses case-insensitive substring matching. For species-level expected organisms (two+ words), genus-level matches alone are not accepted — e.g., "Naegleria sp." does not match "Naegleria fowleri". For genus-only expected organisms, genus-level matching is permitted.

## Output

`blast_false_positive_report.tsv` contains one row per sample/taxid combination with the following columns:

| Column | Description |
|---|---|
| `sample` | Sample ID |
| `taxid` | NCBI taxonomic ID of the target pathogen |
| `expected_organism` | Organism name for the taxid |
| `classification` | TRUE_POSITIVE, FALSE_POSITIVE, UNCERTAIN, UNCULTURED_DOMINANT, NO_DATA, or BLAST_NOT_RUN |
| `avg_pct_match` | Average % of reads matching expected organism across R1 and R2 |
| `avg_pct_uncultured` | Average % of best hits that are uncultured/environmental |
| `R1_n_reads` / `R2_n_reads` | Number of reads BLASTed per read set |
| `R1_n_match` / `R2_n_match` | Number of reads with any hit to the expected organism |
| `R1_pct_match` / `R2_pct_match` | % of reads matching expected organism |
| `R1_pct_uncultured` / `R2_pct_uncultured` | % of best hits that are uncultured |
| `R1_expected_rank` / `R2_expected_rank` | Rank of expected organism by read count (1 = most reads) |
| `R1_mean_pident` / `R2_mean_pident` | Mean % nucleotide identity across all hits |
| `R1_top_organisms` / `R2_top_organisms` | Top 5 organisms by read count with read counts in parentheses |

## Dependencies

- BLAST+ ≥ 2.10.0 (required for NCBI nt v5 database format)
- [KrakenTools](https://github.com/jenniferlu717/KrakenTools) (`extract_kraken_reads.py`)
- seqtk
- Python ≥ 3.8
- pandas

## Input Files

- **Kraken2 output**: `.kraken` and `.kreport` files per sample
- **Paired FASTQ files**: gzipped paired-end reads per sample
- **`sample_taxid_list_all_pathogens.tsv`**: two-column TSV (no header) with sample ID and taxid for each combination to validate
- **`taxid_to_name_all_pathogens.csv`**: two-column CSV with `taxid` and `organism_name` for all target pathogens

## SLURM Configuration

| Step | CPUs | Memory | Time limit |
|---|---|---|---|
| Read extraction | 1 | 16 GB | 4 hours |
| BLASTN | 8 | 32 GB | 24 hours |

Extraction jobs run on `hsph`, `shared`, or `sapphire` partitions. BLAST jobs run on `hsph`.
