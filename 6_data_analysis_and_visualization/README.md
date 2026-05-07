# Biofilm YNHH Project — Metagenomic & 16S Analysis

This repository contains the R analysis pipeline for the Healy Lab Biofilm Project at Yale New Haven Hospital (YNHH). The code integrates shotgun metagenomic data (Kraken2/Bracken) with 16S rRNA amplicon data to characterize the microbial communities and pathogen distributions across sink biofilms, sewer (endcap) biofilms, and hospital wastewater collected from five hospital sites.

## Overview

The pipeline ingests Kraken2/Bracken classifications, FastQC read counts, BLAST false-positive reports, GeCoCheck genome coverage output, sample metadata, and a 16S phyloseq object, then produces a wide range of figures and summary tables describing pathogen detection, community composition, diversity, and differential abundance across sample types and time points.

## Project Structure

```
Biofilm YNHH Project R/
├── Input Data/                 # Raw inputs (Kraken/Bracken TSVs, metadata, BLAST, 16S RDS, etc.)
├── Output Spreadsheets/        # Generated TSV/XLSX summaries
├── Biofilm Project Figures/    # Generated PNG/SVG figures
└── analysis.Rmd                # The R code in this repo
```

Set the working directory and the input data path at the top of the script:

```r
setwd("/Users/amd689/Documents/Healy Lab/Biofilm Project/Biofilm YNHH Project R")
data_path = "Input Data/"
```

## Input Files

The pipeline expects the following files in `Input Data/`:

- `kraken_ct_0_5_min_hit_3_all_samples.tsv` — Kraken2 reports across samples (confidence threshold 0.5, min-hit 3)
- `kraken_ct_0_15_min_hit_2_all_samples.tsv` — Alternate Kraken2 thresholds for stacked-bar plots
- `fastqc_read_counts.tsv` — Per-sample total sequence counts after QC, trimming, and dehosting
- `blast_false_positive_report.tsv` — BLAST validation results used to flag false-positive Kraken hits
- `coverage_checker_output.tsv` — GeCoCheck output (Bowtie2 genome coverage stats per taxon)
- `metadata_Gi_AD.csv` — Sample metadata (SampleID, Location, Type, Date, etc.)
- `16S_Metadata.xlsx` — 16S sample metadata with timepoint and sub-timepoint annotations
- `phyloseq_object_updated.rds` — Phyloseq object containing the 16S OTU table, taxonomy, and sample data
- `Species_Bracken_output.xlsx` — Bracken species-level outputs (generated/consumed by the pipeline)

## Sample Types

| Code     | Label             |
| -------- | ----------------- |
| `Drain`  | Sink Biofilm      |
| `Endcap` | Sewer Biofilm     |
| `H_WW`   | Hospital Wastewater |
| `Mu_WW`  | Municipal Wastewater (filtered out of most analyses) |
| `tap`    | Tap Water (16S only) |

## Sites

Five hospital locations are analyzed and renamed for publication:

| Code   | Display Name |
| ------ | ------------ |
| MARIO  | Site 1       |
| FIONA  | Site 2       |
| LUIGI  | Site 3       |
| SHREK  | Site 4       |
| OSCAR  | Site 5       |

## Pipeline Sections

The `.Rmd` is organized into the following major sections (mirrored by the headers in the file):

1. **Setup** — Working directory, libraries, custom palettes (`fortyfive_pal`, viridis-derived colors), themes, and helper functions including `pairwise.adonis` / `pairwise.adonis2` and a log1p transform.
2. **Relative abundance from Kraken2** — Computes per-sample species fractions using total QC'd reads as the denominator, with a 0.0001% relative-abundance floor.
3. **BLAST false-positive filtering** — Loads the BLAST report and removes flagged taxa from Kraken/Bracken tables.
4. **Genome coverage (GeCoCheck)** — Estimates mean depth per taxon using `(reads_mapped × 151bp) / genome_length` and exports a pathogen coverage summary.
5. **Metadata wrangling** — Cleans dates, recodes locations and timepoints, and merges metadata with abundance tables.
6. **Read-count summary tables** — Total reads annotated by Kraken vs total QC'd reads per sample.
7. **RPIP pathogen panel** — Defines the curated RPIP pathogen list (current and legacy versions, both Illumina nomenclature and Kraken/NCBI nomenclature) and computes detection-frequency and relative-abundance comparisons against non-RPIP species.
8. **Date matching with 16S** — Restricts the metagenomic data to dates with paired 16S samples and selects each non-biofilm sample closest in time to its corresponding biofilm sample (per Site × Timepoint).
9. **NMDS & PERMANOVA** — Bray–Curtis ordination at species/OTU level, faceted plots by sample type and site, and PERMANOVA with explained-variance bar charts and pairwise comparisons.
10. **Alpha diversity** — Rarefaction-based Shannon, Simpson, richness, and evenness, plotted with `ggstatsplot::ggbetweenstats`.
11. **Bray–Curtis dissimilarity by sample type** — Per-Site × Type within-group dissimilarity boxplots.
12. **Venn / Euler diagrams** — Pathogen overlap between Sink Biofilm, Sewer Biofilm, and Hospital Wastewater (per Site, averaged across Sites, and totals).
13. **Core sewer biofilm microbiome** — Stacked bars of the top 44 pathogens by mean relative abundance, prevalence × detection-threshold heatmaps using the `microbiome` package at multiple taxonomic ranks (Genus, Family, Order, Class, Phylum, Kingdom).
14. **16S analyses** — Rarefaction curves, NMDS/PERMANOVA on OTU-level data, stacked bars by Genus / Family / Order / Class / Phylum (faceted and time-collapsed), Venn/Euler overlap, and Bray–Curtis dissimilarity.
15. **Heatmaps** — Top-10 RPIP pathogens by median relative abundance in sewer biofilm, plotted as log10 relative-abundance heatmaps with sample-type facets.
16. **CV analysis** — Coefficient-of-variation categories (Highly Stable → Unstable) per pathogen × Site × Type, paired with relative-abundance boxplots.
17. **UpSet plots** — Pathogen sharing across sample types using `ggupset`, both at Genus and Species level, including a custom dot-and-bar layout for total detections per environment.
18. **Stacked bars (real relative abundance)** — Annotated-RPIP-read fractions split into top species, "Other Enriched Pathogens", and "Non Enriched Taxa".
19. **MaAsLin3 differential abundance** — Mixed-effects model with `Type` as a fixed effect and `Location` as a random effect, producing volcano plots, faceted boxplots of significant pathogens, and result tables.

## Key Outputs

Figures are written to `Biofilm Project Figures/` as both `.png` (300–600 dpi) and `.svg`. Summary tables are written to `Output Spreadsheets/` as `.tsv` and `.xlsx`. Notable outputs include:

- `Output Spreadsheets/kraken_relative_abundance.tsv`
- `Output Spreadsheets/Ginkgo_rpip_pathogen_coverage_summary.tsv`
- `Output Spreadsheets/Ginkgo_rpip_pathogen_persample.tsv`
- `Output Spreadsheets/Total Reads Annotated and Total Reads QCd.xlsx`
- `Output Spreadsheets/kraken_sample_taxid_list_fixed_with_0.0001%_rel_ab_threshold.xlsx`
- `Biofilm Project Figures/Kraken2 CT 0.5 ... .png/.svg` (NMDS, heatmaps, stacked bars, UpSet plots, Euler diagrams, etc.)
- `Biofilm Project Figures/16S ... .png/.svg` (16S equivalents for stacked bars, NMDS, core microbiome, Venn diagrams)
- `maaslin3_output/` (MaAsLin3 model fits, residuals, and results tables)

## Dependencies

### CRAN

`ggh4x`, `ggplot2`, `ggstatsplot`, `dplyr`, `tidyverse`, `readxl`, `vegan`, `OTUtable`, `MASS`, `sysfonts`, `showtext`, `tidytext`, `ggtext`, `RColorBrewer`, `writexl`, `forcats`, `viridis`, `patchwork`, `scales`, `eulerr`, `ggvenn`, `ggupset`, `ComplexUpset`, `ggsignif`, `ggrepel`, `gridGraphics`, `grid`, `svglite`, `cowplot`, `ggpubr`, `rstatix`, `reshape2`

### Bioconductor / GitHub

`phyloseq`, `Biostrings`, `genefilter`, `DESeq2`, `decontam` (Bioconductor)
`microbiome` (`devtools::install_github("microbiome/microbiome")`)
`pairwiseAdonis` (`remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")`)
`microshades`, `ape`, `maaslin3`

Install MaAsLin3 from Bioconductor or its GitHub repo before running section 19.

## Filtering Conventions

- **Species level only**: `taxonomy_lvl == "S"` (or `rank == "S"` for raw Kraken).
- **Relative-abundance floor**: 0.0001 (0.01%) of total QC'd reads when computing Kraken-derived relative abundance.
- **BLAST validation**: Records flagged as `FALSE_POSITIVE` or `UNCERTAIN` are optionally removed.
- **Date-matched subset**: Wastewater samples restricted to dates `2024-09-09`, `09-11`, `09-12`, `10-30`, `10-31`, `12-04`, `12-09`; for each Site × Timepoint × Type, the sample closest to the biofilm collection date is retained.
- **Timepoint windows**:
  - Pre: 2023-12-01 → 2024-04-30
  - Timepoint 1: 2024-08-01 → 2024-09-30
  - Timepoint 2: 2024-10-01 → 2024-11-30
  - Timepoint 3: 2024-12-01 → 2024-12-31

## Reproducing the Analysis

1. Place all input files in `Input Data/`.
2. Set the working directory in the first chunk to your local clone.
3. Run the chunks sequentially. Several chunks depend on objects built earlier (e.g., `bracken_merged_df`, `kraken_merged_df_filtered`, `sixteens_df`, `rpip_pathogen_species_list`).
4. Output directories (`Output Spreadsheets/`, `Biofilm Project Figures/`, `maaslin3_output/`) must exist or will be created by the relevant `ggsave` / `write_xlsx` calls.

## Notes

- Several chunks are wrapped in `eval=FALSE, include=FALSE` (e.g., theme overrides for the UpSet plot) — toggle them on for the corresponding figure styles.
- The script defines two RPIP pathogen lists (current and legacy) plus a Kraken/NCBI-nomenclature variant; the active panel is set near the bottom of the "Count Pathogens" chunk via `rpip_pathogen_species_list = rpip_pathogen_species_list_kraken_nomenclature_old_version`.
- Color palettes (`fortyfive_pal`, `posterCol*`, `viridisCol*`) are defined once and reused across figures for consistency.
- All figures are saved at high DPI (300–600) for publication.
