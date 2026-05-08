

#  SET WORKING DIRECTORY #######################################

```{r}
#Setting working directory for importing data from excel
setwd("/Users/amd689/Documents/Healy Lab/Biofilm Project/Biofilm YNHH Project R")
```


# SET NAME OF PATH TO DATA
```{r}
data_path = "Input Data/"

```
 

#LOAD LIBRARIES
```{r}
rm(list = ls())
#install.packages("ggh4x")
library(ggh4x); packageVersion("ggh4x")
library(ggplot2); packageVersion("ggplot2")
library(ggstatsplot)
library(dplyr); packageVersion("dplyr")
library(tidyverse); packageVersion("tidyverse")
library(readxl)
library(vegan); packageVersion("vegan")
library(OTUtable)
library(MASS)
library(sysfonts)
library(showtext)
library(tidytext)
library(ggtext)
library(RColorBrewer)

#BiocManager::install("decontam")
#BiocManager::install("DESeq2")
#BiocManager::install("genefilter")
#BiocManager::install("phyloseq")
#install.packages("permute")
#install.packages("vegan")
#install.packages("viridis")
#install.packages("patchwork")
#install.packages("ComplexUpset")
#install.packages("ggupset")
#install.packages("ggstatsplot")
#install.packages("ggtext")
#install.packages("forcats")
#install.packages("ggsignif")
library(ggsignif)
library(forcats)
library(ggtext)
library(ggstatsplot)
library(ggupset)
library(ComplexUpset)
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(genefilter)
library(vegan); packageVersion("vegan")
library(DESeq2); packageVersion("DESeq2")
library(ggplot2)
library(scales)
library(dplyr)
library(tidyverse)
library(viridis)
library(patchwork)
#theme_set(theme_bw())
```
#INSTALL FUNCTION FROM GITHUB
```{r eval=FALSE, include=FALSE}
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)

#function pulled from github

#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(devtools)

#devtools::install_github("gauravsk/ranacapa", force = TRUE)

# library(ranacapa)
# 
# data(iris)
# ranacapa::pairwise_adonis(iris[, 1:4], iris$Species)
# 
# ranacapa::pairwise_adonis(iris[, 1:4], iris$Species, reduce = 'setosa')
# 
# # similarity euclidean from vegdist and holm correction
# pairwise_adonis(x = iris[, 1:4], factors = iris$Species,
# sim_method = 'euclidian', p_adjust_m = 'holm')
# 
# #similarity manhattan from daisy and bonferroni correction
# pairwise_adonis(x = iris[, 1:4], factors = iris$Species,
# # sim_method = 'manhattan', p_adjust_m = 'bonferroni')
# 
# 
# install.packages("spaa")
# library(spaa)
# 
# devtools::install_github("GuillemSalazar/EcolUtils")
# 
# library(EcolUtils)
# 
# data(dune)
# data(dune.env)
# adonis.pair(vegdist(dune),dune.env$Management)
```

#COLORS
```{r}

posterCol1=  "#8488AC"
posterCol2 = "#9FB498"
posterCol3 = "#BA6646"
posterCol4 = "#565656"
posterCol5 = "#865856"
posterCol6 = "#B99796"
posterCol7 = "#F8E5DC"

viridisCol1 = "#79307D"
viridisCol2 = "#417C8C"
viridsCol3 = "#E57262"


#colors = c( "#79307D","#79307D","#417C8C","#417C8C","#417C8C","#417C8C","#E57262","#79307D","#79307D","#79307D")

colors = c( "#8488AC","#BA6646","#9FB498")

fortyfive_pal = c("#a0cb6b","#8368cb","#c86c69","#cdd3e5","#dab594","#d692d1",
               "#7495c3","#9fdeca","#e2e8c3","#d8a4af","#71bed2","#bca9dd",
               "#8bb598","#e5cbd4","#6d7ecd","#e4d5ca","#a8dfa5","#a0bada",
               "#cbca6f","#c6926f","#cce7e6","#81a4b2","#ca69c3","#76bd75",
               "#d37dae","#abb684","#6ecda4","#a88bbb","#8cbbb6","#d18191",
               "#c1d8c2","#dac4e6","#e0adce","#a96dca","#e2bdb6","#aacedd",
               "#9992d9","#e0d0ab","#a2abdf","#b88b9e","#b4a382","#dba294",
               "#c9dd9b","#c8b67c","#b9847b")



```



#FUNCTIONS
```{r}

remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")



pairwise.adonis <- function(x,factors, sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  co = combn(unique(factors),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(co[1,elem],co[2,elem]),] ~ factors[factors %in% c(co[1,elem],co[2,elem])] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}



log1p_trans <- function(x,na.rm = FALSE) {
  if (na.rm) {
    x <- x[!is.na(x)]  # Removing NA values if na.rm is TRUE
  }
  return(log10(x+1))
}


#updated adonis2 function

pairwise.adonis2 <- function(x, factors, sim.method = 'bray', p.adjust.m = 'bonferroni') {
  library(vegan)
  
  co <- combn(unique(factors), 2) # all pair combinations
  pairs <- c()
  F.Model <- c()
  R2 <- c()
  p.value <- c()
  
  for (elem in 1:ncol(co)) {
    # Subset the data & factor
    sel <- factors %in% c(co[1, elem], co[2, elem])
    
    # Skip if any group has < 2 samples
    if (min(table(factors[sel])) < 2) {
      next
    }
    
    ad <- adonis2(
      x[sel, ] ~ factors[sel],
      method = sim.method
    )
    
    pairs     <- c(pairs, paste(co[1, elem], 'vs', co[2, elem]))
    F.Model   <- c(F.Model, ad$F[1])
    R2        <- c(R2, ad$R2[1])
    p.value   <- c(p.value, ad$`Pr(>F)`[1])
  }
  
  p.adjusted <- p.adjust(p.value, method = p.adjust.m)
  pairw.res  <- data.frame(pairs, F.Model, R2, p.value, p.adjusted)
  
  return(pairw.res)
}


```


#SET THEMES
```{r}
library(ggplot2)

# Reset global theme and font
theme_set(theme_grey(base_family = "", base_size = 11))

# Remove theme_update() overrides
theme_update(
  axis.text.x        = element_text(),
  axis.text.y        = element_text(),
  axis.ticks.x       = element_line(),
  axis.ticks.y       = element_line(),
  axis.ticks.length.x = unit(0.15, "lines"),
  panel.grid         = element_line(),
  plot.background    = element_rect()
)

# Reset geom defaults that may have inherited "Chivo"
update_geom_defaults("text",  list(family = ""))
update_geom_defaults("label", list(family = ""))

showtext::showtext_auto(FALSE)

```




#Get rel abundances from kraken2 output (out of total reads approach)
```{r}
library(tidyverse)

# ---- Parameters ----
INPUT     <- "Input Data/kraken_ct_0_5_min_hit_3_all_samples.tsv"
FASTQC    <- "Input Data/fastqc_read_counts.tsv"
OUTPUT    <- "Output Spreadsheets/kraken_relative_abundance.tsv"

# ---- Load data ----
kraken <- read_tsv(INPUT, col_types = cols(
  SampleID     = col_character(),
  pct_reads    = col_double(),
  reads_clade  = col_double(),
  reads_direct = col_double(),
  rank         = col_character(),
  taxid        = col_double(),
  name         = col_character()
))

fastqc <- read_tsv(FASTQC, col_types = cols(
  SampleID         = col_character(),
  total_sequences  = col_double()
))



# ---- Filter to species level ----
kraken_rel_ab <- kraken %>%
  filter(rank == "S") %>%
  left_join(fastqc, by = "SampleID") %>%
  mutate(fraction_total_reads = reads_clade / total_sequences) %>% filter(fraction_total_reads > 0.000001) # filter out 0.0001% (0.000001)

# ---- Sanity check ----
# Values won't sum to 1 — that's expected since denominator is ALL reads
# including unclassified and reads assigned above species level
kraken_rel_ab %>%
  group_by(SampleID) %>%
  summarise(total = sum(fraction_total_reads)) %>%
  print(n = Inf)


```



#Look at BLAST results
```{r}


# Read data
blast_report_df <- read.table("Input Data/blast_false_positive_report.tsv",
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE
) %>% 
  
  filter(classification != "BLAST_NOT_RUN") %>% 
  
  dplyr::rename("name" = "expected_organism",
                "SampleID" = "sample",
                "taxonomy_id" = "taxid") %>% 
  
  mutate(taxonomy_id = as.integer(taxonomy_id))




```



#Get Genome Coverage from GeCoCheck
```{r}
library(tidyverse)

READ_LENGTH <- 151  # bp

# Read data
df <- read_tsv("Input Data/coverage_checker_output.tsv",
               col_types = cols(
                 `Bowtie2 mean identity of mapped reads (%)` = col_double(),
                 `Bowtie2 mean MAPQ of mapped reads` = col_double()
               )) %>% 

# Clean column names for easier use

  dplyr::rename(
    Sample                  = Sample,
    taxid                   = taxid,
    Species_name            = `Species name`,
    Genome_length_bp        = `Reference genome length (bp)`,
    Kraken_reads            = `Kraken reads assigned`,
    Bowtie2_reads_mapped    = `Bowtie2 reads mapped`,
    Prop_kraken_mapped      = `Proportion kraken reads mapped with Bowtie2`,
    Bowtie2_mean_identity   = `Bowtie2 mean identity of mapped reads (%)`,
    Bowtie2_mean_MAPQ       = `Bowtie2 mean MAPQ of mapped reads`,
    Bowtie2_genome_fraction = `Bowtie2 genome fraction (%)`
  ) %>% 

# Estimate mean depth per sample:
# depth = (reads_mapped * read_length) / genome_length

  mutate(est_mean_depth = (Bowtie2_reads_mapped * READ_LENGTH) / Genome_length_bp)

# ── Summary per taxon across all samples ──────────────────────────────────────
taxon_summary <- df %>%
  group_by(taxid, Species_name, Genome_length_bp) %>%
  summarise(
    mean_est_depth          = mean(est_mean_depth, na.rm = TRUE),
    median_est_depth        = median(est_mean_depth, na.rm = TRUE),
    mean_genome_fraction    = mean(Bowtie2_genome_fraction, na.rm = TRUE),
    median_genome_fraction  = median(Bowtie2_genome_fraction, na.rm = TRUE),
    mean_bowtie2_reads      = mean(Bowtie2_reads_mapped, na.rm = TRUE),
    median_bowtie2_reads    = median(Bowtie2_reads_mapped, na.rm = TRUE),
    n_samples_in_output     = n(),           # samples where Kraken detected it
    n_samples_total   = 250,
    n_samples_detected      = sum(Bowtie2_reads_mapped > 0),
    .groups = "drop"
  ) %>%
  mutate(
    pct_of_kraken_positive  = round(100 * n_samples_detected / n_samples_in_output, 1),
    pct_of_all_samples      = round(100 * n_samples_detected / 250, 1)  # <-- use this for prevalence
  ) %>%
  arrange(desc(mean_est_depth))

print(taxon_summary, n = Inf)
write_tsv(taxon_summary, "Output Spreadsheets/Ginkgo_rpip_taxon_coverage_summary.tsv")

# ── Optional: filter to known pathogens ───────────────────────────────────────
# Edit this list to match your pathogens of interest 
#CHANGE THIS TO RPIP LIST
pathogen_taxids <- c(
  1764, #Mycobacterium avium
  1773, #Mycobacterium tuberculosis
  666, # Vibrio cholerae
  13373, # Burkholderia mallei
  40324, #Stenotrophomonas maltophilia
  648, # Aeromonas caviae
  5763, # Naegleria fowleri
  644, #Aeromonas hydrophila
  1150, # Citrobacter freundii
  287,    # Pseudomonas aeruginosa
  582,    # Morganella morganii
  648,    # Aeromonas caviae
  651,    # Aeromonas media
  729,    # Haemophilus parainfluenzae
  817,    # Bacteroides fragilis
  106654, # Acinetobacter nosocomialis
  202951, # Acinetobacter bouvetii
  13689,  # Sphingomonas paucimobilis
  255507, # Aliarcobacter cibarius
  28198   # Aliarcobacter cryaerophilus
)

pathogen_summary <- taxon_summary %>%
  filter(taxid %in% pathogen_taxids)

cat("\n── Pathogen summary ──\n")
print(pathogen_summary, n = Inf)
write_tsv(pathogen_summary, "Output Spreadsheets/Ginkgo_rpip_pathogen_coverage_summary.tsv")

# ── Per-sample depth for pathogens (for boxplots etc.) ────────────────────────
pathogen_persample <- df %>%
  filter(taxid %in% pathogen_taxids) %>%
  dplyr::select(Sample, taxid, Species_name, Bowtie2_reads_mapped,
         est_mean_depth, Bowtie2_genome_fraction, Prop_kraken_mapped)

write_tsv(pathogen_persample, "Output Spreadsheets/Ginkgo_rpip_pathogen_persample.tsv")
```



#Load in Metadata
```{r}

metadata_df = read.csv("Input Data/metadata_Gi_AD.csv",header = T) %>% 
  
  #change name of Location
  mutate(Location =  sub("^[0-9]+\\.", "", Location)) %>% 
  
  #change month column
  mutate(Month = sub("^[0-9]+_", "", Month)) %>% 

  #change date column
  mutate(Date = sub("^[0-9]+_", "", Date)) %>% 
  mutate(Date = mdy(paste0(Date, " 2024"))) %>% 
  
  dplyr::select(-LocationMonth,-LocationType) %>% 
  
  dplyr::rename("SampleID"="Sample_id")

check_counts_metadata = metadata_df %>% 
  
  dplyr::select("SampleID") %>% distinct()

#municipal versus non municipal list------------------



```


#Merge data
```{r}

#get kracken_merged_df ----------------
kraken_merged_df <- left_join(kraken_rel_ab, metadata_df, by = 'SampleID') %>% 
  distinct() %>% 
  ungroup() %>% 
  
  dplyr::rename("taxonomy_lvl" = "rank",
                "taxonomy_id" = "taxid") %>% 
  
  filter(Type != "Mu_WW")

# high_abun_species = kraken_merged_df %>% 
#   
#   group_by(name) %>% 
#   summarise(fraction_total_reads = mean(fraction_total_reads,na.rm = T)) %>% 
#   ungroup() %>% 
#   
#   #filter out species at low abundance that counting above a specific relativea abundance as true detects
#   filter(fraction_total_reads > 0.0001) %>% 
#   
#   dplyr::select(name) %>% distinct() %>% pull()
  
  
# kraken_merged_df <- kraken_merged_df %>% 
#   filter(name %in% high_abun_species) 

library(writexl)

# kraken_filtered_for_pathogens= kraken_merged_df %>% filter(name %in% rpip_pathogen_species_list)
# 
# write_xlsx(kraken_filtered_for_pathogens,"Input Data/only_pathogens_kraken_ct_0_50_min_hit_3_filtered_out_less_than_0.001%relab_and_used_total_reads_assigned_qcd_and_trimmed_and_dehosted_as_denominator.xlsx")
# 

# all_samples_and_taxid = kraken_filtered_for_pathogens %>% 
#   dplyr::select(SampleID,taxonomy_id) %>% distinct()
# 
# 
# all_pathogens_and_taxid = kraken_filtered_for_pathogens %>% 
#   dplyr::select(name,taxonomy_id) %>% distinct()
#   

check_counts_merged = kraken_merged_df %>% 
  
  filter(Type %in% c("Endcap", "H_WW","Drain")) %>% 
  
  dplyr::select("SampleID") %>% distinct()

```

#Make Read count table from FastQC
```{r}
library(tidyverse)

# ---- Parameters ----
INPUT     <- "Input Data/kraken_ct_0_5_min_hit_3_all_samples.tsv"
FASTQC    <- "Input Data/fastqc_read_counts.tsv"


# ---- Load data ----
kraken <- read_tsv(INPUT, col_types = cols(
  SampleID     = col_character(),
  pct_reads    = col_double(),
  reads_clade  = col_double(),
  reads_direct = col_double(),
  rank         = col_character(),
  taxid        = col_double(),
  name         = col_character()
))

fastqc <- read_tsv(FASTQC, col_types = cols(
  SampleID         = col_character(),
  total_sequences  = col_double()
))



# ---- Get total read counts annotated by kraken and total sequences from fastqc ----
reads_total_df <- kraken %>%
  filter(rank == "S") %>%
  left_join(fastqc, by = "SampleID") %>% 
  
  distinct(SampleID, name, reads_clade,total_sequences) %>% 
  
  group_by(SampleID,total_sequences) %>%
  summarize(total_reads_clade = sum(reads_clade,na.rm = T)) 

# ---- Merge with metadata to filter out municipal ww samples and output as table ----

merged_table_df <- left_join(reads_total_df, metadata_df, by = 'SampleID') %>% 
  distinct() %>% 
  ungroup() %>% 
  
  filter(Type != "Mu_WW") %>% 
  
  dplyr::select(Date,Location,Type,total_reads_clade, total_sequences)  %>% 
  
  
  mutate(
    timepoint = case_when(
      Date >= as.Date("2023-12-01") & Date <= as.Date("2024-04-30") ~ "Pre",
      Date >= as.Date("2024-08-01") & Date <= as.Date("2024-09-30") ~ "1",
      Date >= as.Date("2024-10-01") & Date <= as.Date("2024-11-30") ~ "2",
      Date >= as.Date("2024-12-01") & Date <= as.Date("2024-12-31") ~ "3",
      TRUE ~ NA_character_
    )
  ) %>%

  
  mutate(Type = case_when(Type == "Endcap" ~ "Sewer Biofilm",
                          Type == "H_WW" ~ "Wastewater",
                          Type == "Drain" ~ "Sink Biofilm")) %>% 
  
    # now replace Site names 
  mutate(Location = case_when(Location == "MARIO" ~ "Site 1",
                              Location == "FIONA" ~ "Site 2",
                              Location == "LUIGI" ~ "Site 3",
                              Location == "SHREK" ~ "Site 4",
                              Location == "OSCAR" ~ "Site 5")) %>% 
  
  dplyr::rename(
    "Timepoint" = "timepoint",
    "Total Reads Annotated by Kraken"= "total_reads_clade",
    "Total Reads After Trimming, Dehosting, and QC" = "total_sequences"
  )

write_xlsx(merged_table_df, "Output Spreadsheets/Total Reads Annotated and Total Reads QCd.xlsx")

#get output with sample IDs
nonmunicipal_sampleids_df <- left_join(reads_total_df, metadata_df, by = 'SampleID') %>% 
  distinct() %>% 
  ungroup() %>% 
  
  filter(Type != "Mu_WW") %>% 
  
  dplyr::select(SampleID,total_reads_clade, total_sequences)  %>% 
  
  
  dplyr::rename(
    "Total Reads Annotated by Kraken"= "total_reads_clade",
    "Total Reads After Trimming, Dehosting, and QC" = "total_sequences"
  )


write_xlsx(nonmunicipal_sampleids_df, "Output Spreadsheets/NonmunicipaSamples.xlsx")

#get municipal  sample IDs
municipal_sampleids_df <- left_join(reads_total_df, metadata_df, by = 'SampleID') %>% 
  distinct() %>% 
  ungroup() %>% 
  
  filter(Type == "Mu_WW") %>% 
  
  dplyr::select(SampleID,total_reads_clade, total_sequences)  %>% 
  
  
  dplyr::rename(
    "Total Reads Annotated by Kraken"= "total_reads_clade",
    "Total Reads After Trimming, Dehosting, and QC" = "total_sequences"
  )


write_xlsx(municipal_sampleids_df, "Output Spreadsheets/MunicipaSamples.xlsx")

```



#Count Pathogens and then Calculate percentage of species that are RPIP vs not
```{r}
#---old list----
old_list <-  c(
  "Achromobacter denitrificans","Achromobacter xylosoxidans","Acinetobacter baumannii",
  "Acinetobacter lwoffii","Acinetobacter nosocomialis","Acinetobacter pittii",
  "Actinomyces graevenitzii","Actinomyces israelii","Actinomyces meyeri",
  "Actinomyces naeslundii","Actinomyces odontolyticus","Aeromonas caviae",
  "Aeromonas hydrophila","Aeromonas sobria","Aeromonas veronii",
  "Aggregatibacter actinomycetemcomitans","Aggregatibacter aphrophilus",
  "Arcanobacterium haemolyticum","Bacillus anthracis","Bacillus cereus",
  "Bacillus thuringiensis","Bacteroides fragilis","Bartonella henselae",
  "Bartonella quintana","Bordetella bronchiseptica","Bordetella hinzii",
  "Bordetella holmesii","Bordetella parapertussis","Bordetella pertussis",
  "Brucella abortus","Brucella canis","Brucella melitensis","Brucella suis",
  "Burkholderia gladioli","Burkholderia mallei","Burkholderia pseudomallei",
  "Burkholderia thailandensis","Campylobacter concisus","Capnocytophaga gingivalis",
  "Capnocytophaga leadbetteri","Capnocytophaga sputigena","Cardiobacterium hominis",
  "Cardiobacterium valvarum","Chlamydia pneumoniae","Chlamydia psittaci",
  "Chlamydia trachomatis","Citrobacter freundii","Citrobacter koseri",
  "Corynebacterium diphtheriae","Corynebacterium jeikeium","Corynebacterium striatum",
  "Corynebacterium ulcerans","Coxiella burnetii","Cronobacter sakazakii",
  "Cutibacterium propionicum","Eikenella corrodens","Elizabethkingia anophelis",
  "Elizabethkingia meningoseptica","Enterobacter cloacae","Enterococcus faecalis",
  "Enterococcus faecium","Escherichia coli","Fusobacterium necrophorum",
  "Fusobacterium nucleatum","Gemella haemolysans","Gemella morbillorum",
  "Haemophilus influenzae","Haemophilus parahaemolyticus","Haemophilus parainfluenzae",
  "Haemophilus pittmaniae","Haemophilus haemolyticus","Klebsiella aerogenes",
  "Klebsiella oxytoca","Klebsiella pneumoniae","Klebsiella quasipneumoniae",
  "Klebsiella variicola","Kingella kingae","Leptotrichia buccalis","Listeria monocytogenes",
  "Moraxella catarrhalis","Moraxella osloensis","Morganella morganii",
  "Mycobacterium tuberculosis complex","Mycobacterium avium","Mycoplasma pneumoniae",
  "Neisseria flavescens","Neisseria lactamica","Neisseria meningitidis",
  "Neisseria mucosa","Pasteurella multocida","Proteus mirabilis","Proteus vulgaris",
  "Providencia stuartii","Pseudomonas aeruginosa","Pseudomonas fluorescens",
  "Pseudomonas stutzeri","Ralstonia pickettii","Raoultella ornithinolytica",
  "Raoultella planticola","Salmonella enterica","Serratia marcescens",
  "Staphylococcus aureus","Stenotrophomonas maltophilia","Streptococcus agalactiae",
  "Streptococcus anginosus","Streptococcus constellatus","Streptococcus dysgalactiae",
  "Streptococcus intermedius","Streptococcus mitis","Streptococcus pneumoniae",
  "Streptococcus pyogenes","Treponema denticola","Ureaplasma urealyticum",
  "Ureaplasma parvum","Veillonella parvula","Yersinia enterocolitica","Yersinia pestis"
)
#-----------Get pathogen lists--------------

#new RPIP list
rpip_pathogen_species_list_illumina_terms_new_version <- c(
  # --- Viruses ---
  "Coxsackievirus A", "Coxsackievirus B", "Cytomegalovirus (CMV)", "EBV", 
  "Enterovirus A71", "Enterovirus D68", "Herpes simplex virus 1 (HSV-1)", 
  "Human adenovirus B", "Human adenovirus C", "Human adenovirus E", 
  "Human bocavirus 1", "Human Coronavirus 229E", "Human Coronavirus HKU1", 
  "Human Coronavirus NL63", "Human Coronavirus OC43", "Human metapneumovirus", 
  "Human parainfluenza virus 1", "Human parainfluenza virus 2", 
  "Human parainfluenza virus 3", "Human parainfluenza virus 4", 
  "Human parechovirus", "Influenza A virus (H1N1)", "Influenza A virus (H3N2)", 
  "Influenza A virus (H5N1)", "Influenza A virus (H7N9)", "Influenza A virus (H9N2)", 
  "Influenza B virus", "Influenza C virus", "Measles Virus", 
  "MERS coronavirus (MERS-CoV)", "Mumps virus", "Parvovirus B19", 
  "Respiratory Syncytial Virus A", "Respiratory Syncytial Virus B", 
  "Rhinovirus A", "Rhinovirus B", "Rhinovirus C", "Rubella virus", 
  "SARS coronavirus", "SARS-CoV-2 (2019-nCoV)", "Varicella-zoster virus (HHV-3)",
  
  # --- Bacteria ---
  "Achromobacter denitrificans", "Bartonella henselae", "Cardiobacterium hominis", "Elizabethkingia meningoseptica", "Haemophilus influenzae", "Mycobacterium gordonae",
  "Achromobacter xylosoxidans", "Bartonella quintana", "Cardiobacterium valvarum", "Enterobacter cloacae complex", "Haemophilus parahaemolyticus", "Mycobacterium kansasii",
  "Acinetobacter baumannii", "Bordetella bronchiseptica", "Chlamydia pneumoniae", "Enterococcus faecalis", "Haemophilus parainfluenzae", "Mycobacterium malmoense",
  "Acinetobacter lwoffii", "Bordetella hinzii", "Chlamydia psittaci", "Enterococcus faecium", "Haemophilus pittmaniae", "Mycobacterium parascrofulaceum",
  "Acinetobacter nosocomialis", "Bordetella holmesii", "Chlamydia trachomatis", "Escherichia coli", "Hafnia alvei", "Mycobacterium scrofulaceum",
  "Acinetobacter pittii", "Bordetella parapertussis", "Chromobacterium violaceum", "Eubacterium brachy", "Klebsiella variicola", "Mycobacterium szulgai",
  "Actinomyces graevenitzii", "Bordetella pertussis", "Citrobacter freundii", "Eubacterium limosum", "Kytococcus sedentarius", "Mycobacterium tuberculosis",
  "Actinomyces israelii", "Bordetella petrii", "Citrobacter koseri", "Eubacterium nodatum", "Leclercia adecarboxylata", "Mycobacterium xenopi",
  "Actinomyces meyeri", "Brucella abortus", "Corynebacterium diphtheriae", "Finegoldia magna", "Legionella anisa", "Mycobacteroides abscessus",
  "Actinomyces naeslundii", "Brucella canis", "Corynebacterium jeikeium", "Francisella tularensis", "Legionella feeleii", "Mycobacteroides chelonae",
  "Actinomyces odontolyticus", "Brucella melitensis", "Corynebacterium propinquum", "Fusobacterium necrophorum", "Legionella longbeachae", "Mycobacteroides immunogenum",
  "Aeromonas caviae", "Brucella suis", "Corynebacterium pseudodiphtheriticum", "Fusobacterium nucleatum", "Legionella maceachernii", "Mycoplasma pneumoniae",
  "Aeromonas hydrophila", "Burkholderia cepacia complex", "Corynebacterium pseudotuberculosis", "Gemella haemolysans", "Legionella pneumophila", "Neisseria flavescens",
  "Aeromonas sobria", "Burkholderia gladioli", "Corynebacterium striatum", "Gemella morbillorum", "Legionella wadsworthii", "Neisseria lactamica",
  "Aeromonas veronii", "Burkholderia glumae", "Corynebacterium ulcerans", "Gordonia araii", "Leptospira interrogans", "Neisseria meningitidis",
  "Aggregatibacter actinomycetemcomitans", "Burkholderia mallei", "Coxiella burnetii", "Gordonia bronchialis", "Leptotrichia buccalis", "Neisseria mucosa",
  "Aggregatibacter aphrophilus", "Burkholderia pseudomallei", "Cronobacter sakazakii", "Haemophilus haemolyticus", "Listeria monocytogenes", "Nocardia abscessus",
  "Arcanobacterium haemolyticum", "Burkholderia thailandensis", "Delftia acidovorans", "Kingella kingae", "Moraxella catarrhalis", "Nocardia arthritidis",
  "Bacillus anthracis", "Campylobacter concisus", "Dialister pneumosintes", "Klebsiella aerogenes", "Moraxella osloensis", "Nocardia beijingensis",
  "Bacillus cereus", "Capnocytophaga gingivalis", "Dolosigranulum pigrum", "Klebsiella oxytoca", "Morganella morganii", "Nocardia brasiliensis",
  "Bacillus thuringiensis", "Capnocytophaga leadbetteri", "Eikenella corrodens", "Klebsiella pneumoniae", "Mycobacterium avium complex", "Nocardia cyriacigeorgica",
  "Bacteroides fragilis", "Capnocytophaga sputigena", "Elizabethkingia anophelis", "Klebsiella quasipneumoniae", "Mycobacterium fortuitum", "Nocardia farcinica",
  "Nocardia nova", "Pediococcus acidilactici", "Pseudomonas aeruginosa", "Serratia marcescens", "Streptococcus intermedius", "Ureaplasma urealyticum",
  "Nocardia otitidiscaviarum", "Peptostreptococcus anaerobius", "Pseudomonas fluorescens", "Shewanella putrefaciens", "Streptococcus mitis", "Veilonella parvula",
  "Nocardia transvalensis", "Prevotella buccae", "Pseudomonas stutzeri", "Slackia exigua", "Streptococcus pneumoniae", "Williamsia muralis",
  "Nocardia veterana", "Prevotella intermedia", "Ralstonia pickettii", "Sphingomonas paucimobilis", "Streptococcus pyogenes", "Yersinia enterocolitica",
  "Ochrobactrum anthropi", "Prevotella melaninogenica", "Raoultella ornithinolytica", "Staphylococcus aureus", "Tatlockia micdadei", "Yersinia pestis",
  "Orientia tsutsugamushi", "Prevotella pleuritidis", "Raoultella planticola", "Stenotrophomonas maltophilia", "Treponema denticola",
  "Pandoraea pulmonicola", "Proteus mirabilis", "Rhodococcus hoagii", "Streptococcus agalactiae", "Tropheryma whipplei",
  "Pantoea agglomerans", "Proteus penneri", "Rickettsia rickettsii", "Streptococcus anginosus", "Tsukamurella pulmonis",
  "Parvimonas micra", "Proteus vulgaris", "Rothia mucilaginosa", "Streptococcus constellatus", "Tsukamurella tyrosinosolvens",
  "Pasteurella multocida", "Providencia stuartii", "Salmonella enterica", "Streptococcus dysgalactiae", "Ureaplasma parvum",

  # --- Fungi ---
"Alternaria alternata", "Blastomyces dermatitidis", "Curvularia lunata", "Lomentospora prolificans", 
  "Pneumocystis jirovecii", "Sarocladium kiliense", "Alternaria infectoria", "Candida auris", 
  "Exophiala dermatitidis", "Microascus cinereus", "Purpureocillium lilacinum", "Scedosporium apiospermum", 
  "Apophysomyces elegans", "Cladophialophora bantiana", "Fusarium oxysporum", "Microascus cirrosus", 
  "Rasamsonia aegroticola", "Schizophyllum commune", "Aspergillus flavus", "Coccidioides immitis", 
  "Fusarium proliferatum", "Microascus paisii", "Rasamsonia argillacea", "Scopulariopsis brevicaulis", 
  "Aspergillus fumigatus", "Coccidioides posadasii", "Fusarium solani", "Mucor circinelloides", 
  "Rhizomucor pusillus", "Sporothrix schenckii", "Aspergillus nidulans", "Cryptococcus gattii", 
  "Fusarium verticillioides", "Mucor indicus", "Rhizopus azygosporus", "Syncephalastrum racemosum", 
  "Aspergillus niger", "Cryptococcus neoformans", "Histoplasma capsulatum", "Mucor racemosus", 
  "Rhizopus microsporus", "Talaromyces marneffei", "Aspergillus terreus", "Cunninghamella bertholletiae", 
  "Lichtheimia corymbifera", "Paecilomyces variotii", "Rhizopus oryzae", "Trichosporon asahii", 
  "Aspergillus versicolor", "Curvularia geniculata", "Lichtheimia ramosa", "Paracoccidioides brasiliensis", 
  "Saksenaea vasiformis"
)

# RPIP Pathogen List: Kraken2/NCBI Nomenclature Version
rpip_pathogen_species_list_illumina_terms_new_version_kraken_nomenclature <- c(
  # --- Viruses ---
  # Kraken typically classifies these by species or specific strain if the DB is custom
  "Coxsackievirus A", "Coxsackievirus B", "Cytomegalovirus", "Epstein-Barr virus", 
  "Enterovirus A71", "Enterovirus D68", "Human alphaherpesvirus 1", # HSV-1
  "Human adenovirus B", "Human adenovirus C", "Human adenovirus E", 
  "Human bocavirus 1", "Human coronavirus 229E", "Human coronavirus HKU1", 
  "Human coronavirus NL63", "Human coronavirus OC43", "Human metapneumovirus", 
  "Human parainfluenza virus 1", "Human parainfluenza virus 2", 
  "Human parainfluenza virus 3", "Human parainfluenza virus 4", 
  "Human parechovirus", "Influenza A virus", "Influenza B virus", "Influenza C virus",
  "Measles morbillivirus", "Middle East respiratory syndrome-related coronavirus", 
  "Mumps orthorubulavirus", "Primate erythroparvovirus 1", # Parvovirus B19
  "Human orthopneumovirus", # RSV
  "Human rhinovirus A", "Human rhinovirus B", "Human rhinovirus C", 
  "Rubella virus", "Severe acute respiratory syndrome-related coronavirus", 
  "Severe acute respiratory syndrome coronavirus 2", "Human alphaherpesvirus 3", # VZV

  # --- Bacteria ---
  "Achromobacter denitrificans", "Achromobacter xylosoxidans", "Acinetobacter baumannii",
  "Acinetobacter lwoffii", "Acinetobacter nosocomialis", "Acinetobacter pittii",
  "Actinomyces graevenitzii", "Actinomyces israelii", "Actinomyces meyeri",
  "Actinomyces naeslundii", "Actinomyces odontolyticus", "Aeromonas caviae",
  "Aeromonas hydrophila", "Aeromonas sobria", "Aeromonas veronii",
  "Aggregatibacter actinomycetemcomitans", "Aggregatibacter aphrophilus",
  "Arcanobacterium haemolyticum", "Bacillus anthracis", "Bacillus cereus",
  "Bacillus thuringiensis", "Bacteroides fragilis", "Bartonella henselae",
  "Bartonella quintana", "Bordetella bronchiseptica", "Bordetella hinzii",
  "Bordetella holmesii", "Bordetella parapertussis", "Bordetella pertussis",
  "Bordetella petrii", "Brucella abortus", "Brucella canis", "Brucella melitensis",
  "Brucella suis", "Burkholderia cepacia complex", "Burkholderia gladioli",
  "Burkholderia glumae", "Burkholderia mallei", "Burkholderia pseudomallei",
  "Burkholderia thailandensis", "Campylobacter concisus", "Capnocytophaga gingivalis",
  "Capnocytophaga leadbetteri", "Capnocytophaga sputigena", "Cardiobacterium hominis",
  "Cardiobacterium valvarum", "Chlamydia pneumoniae", "Chlamydia psittaci",
  "Chlamydia trachomatis", "Chromobacterium violaceum", "Citrobacter freundii complex",
  "Citrobacter koseri", "Corynebacterium diphtheriae", "Corynebacterium jeikeium",
  "Corynebacterium propinquum", "Corynebacterium pseudodiphtheriticum",
  "Corynebacterium pseudotuberculosis", "Corynebacterium striatum",
  "Corynebacterium ulcerans", "Coxiella burnetii", "Cronobacter sakazakii",
  "Delftia acidovorans", "Dialister pneumosintes", "Dolosigranulum pigrum",
  "Eikenella corrodens", "Elizabethkingia anophelis", "Elizabethkingia meningoseptica",
  "Enterobacter cloacae complex", "Enterococcus faecalis", "Enterococcus faecium",
  "Escherichia coli", "Eubacterium brachy", "Eubacterium limosum",
  "Eubacterium nodatum", "Finegoldia magna", "Francisella tularensis",
  "Fusobacterium necrophorum", "Fusobacterium nucleatum", "Gemella haemolysans",
  "Gemella morbillorum", "Gordonia araii", "Gordonia bronchialis",
  "Haemophilus haemolyticus", "Haemophilus influenzae", "Haemophilus parahaemolyticus",
  "Haemophilus parainfluenzae", "Haemophilus pittmaniae", "Hafnia alvei",
  "Kingella kingae", "Klebsiella aerogenes", "Klebsiella oxytoca",
  "Klebsiella pneumoniae", "Klebsiella quasipneumoniae", "Klebsiella variicola",
  "Kytococcus sedentarius", "Leclercia adecarboxylata", "Legionella anisa",
  "Legionella feeleii", "Legionella longbeachae", "Legionella maceachernii",
  "Legionella pneumophila", "Legionella wadsworthii", "Leptospira interrogans",
  "Leptotrichia buccalis", "Listeria monocytogenes", "Moraxella catarrhalis",
  "Moraxella osloensis", "Morganella morganii", "Mycobacterium avium complex",
  "Mycolicibacterium fortuitum", "Mycobacterium gordonae", "Mycobacterium kansasii",
  "Mycobacterium malmoense", "Mycobacterium scrofulaceum", "Mycobacterium szulgai",
  "Mycobacterium tuberculosis complex", "Mycobacterium xenopi",
  "Mycobacteroides abscessus", "Mycobacteroides chelonae", "Mycobacteroides immunogenum",
  "Mycoplasma pneumoniae", "Neisseria flavescens", "Neisseria lactamica",
  "Neisseria meningitidis", "Neisseria mucosa", "Nocardia abscessus",
  "Nocardia arthritidis", "Nocardia beijingensis", "Nocardia brasiliensis",
  "Nocardia cyriacigeorgica", "Nocardia farcinica", "Nocardia nova",
  "Nocardia otitidiscaviarum", "Nocardia transvalensis", "Nocardia veterana",
  "Ochrobactrum anthropi", "Orientia tsutsugamushi", "Pandoraea pulmonicola",
  "Pantoea agglomerans", "Parvimonas micra", "Pasteurella multocida",
  "Pediococcus acidilactici", "Peptostreptococcus anaerobius", "Prevotella buccae",
  "Prevotella intermedia", "Prevotella melaninogenica", "Prevotella pleuritidis",
  "Proteus mirabilis", "Proteus penneri", "Proteus vulgaris",
  "Providencia stuartii", "Pseudomonas aeruginosa", "Pseudomonas fluorescens",
  "Pseudomonas stutzeri", "Ralstonia pickettii", "Raoultella ornithinolytica",
  "Raoultella planticola", "Mycolicibacterium hoagii", "Rickettsia rickettsii",
  "Rothia mucilaginosa", "Salmonella enterica", "Serratia marcescens",
  "Shewanella putrefaciens", "Slackia exigua", "Sphingomonas paucimobilis",
  "Staphylococcus aureus", "Stenotrophomonas maltophilia", "Streptococcus agalactiae",
  "Streptococcus anginosus", "Streptococcus constellatus", "Streptococcus dysgalactiae",
  "Streptococcus intermedius", "Streptococcus mitis", "Streptococcus pneumoniae",
  "Streptococcus pyogenes", "Legionella micdadei", "Treponema denticola",
  "Tropheryma whipplei", "Tsukamurella pulmonis", "Tsukamurella tyrosinosolvens",
  "Ureaplasma parvum", "Ureaplasma urealyticum", "Veillonella parvula",
  "Williamsia muralis", "Yersinia enterocolitica", "Yersinia pestis",

  # --- Fungi ---
  "Alternaria alternata", "Blastomyces dermatitidis", "Curvularia lunata", 
  "Lomentospora prolificans", "Pneumocystis jirovecii", "Sarocladium kiliense", 
  "Alternaria infectoria", "Candida auris", "Exophiala dermatitidis", 
  "Microascus cinereus", "Purpureocillium lilacinum", "Scedosporium apiospermum", 
  "Apophysomyces elegans", "Cladophialophora bantiana", "Fusarium oxysporum", 
  "Microascus cirrosus", "Rasamsonia aegroticola", "Schizophyllum commune", 
  "Aspergillus flavus", "Coccidioides immitis", "Fusarium proliferatum", 
  "Microascus paisii", "Rasamsonia argillacea", "Scopulariopsis brevicaulis", 
  "Aspergillus fumigatus", "Coccidioides posadasii", "Fusarium solani", 
  "Mucor circinelloides", "Rhizomucor pusillus", "Sporothrix schenckii", 
  "Aspergillus nidulans", "Cryptococcus gattii", "Fusarium verticillioides", 
  "Mucor indicus", "Rhizopus azygosporus", "Syncephalastrum racemosum", 
  "Aspergillus niger", "Cryptococcus neoformans", "Histoplasma capsulatum", 
  "Mucor racemosus", "Rhizopus microsporus", "Talaromyces marneffei", 
  "Aspergillus terreus", "Cunninghamella bertholletiae", "Lichtheimia corymbifera", 
  "Paecilomyces variotii", "Rhizopus oryzae", "Trichosporon asahii", 
  "Aspergillus versicolor", "Curvularia geniculata", "Lichtheimia ramosa", 
  "Paracoccidioides brasiliensis", "Saksenaea vasiformis"
)

#old RPIP list (our version)

# RPIP Pathogen Panel Vector
rpip_pathogen_species_list_illumina_terms_old_version <- c(
  "Achromobacter denitrificans", "Achromobacter xylosoxidans", "Acinetobacter baumannii", 
  "Acinetobacter lwoffii", "Acinetobacter nosocomialis", "Acinetobacter pittii", 
  "Actinomyces graevenitzii", "Actinomyces israelii", "Actinomyces meyeri", 
  "Actinomyces naeslundii", "Actinomyces odontolyticus", "Aeromonas caviae", 
  "Aeromonas hydrophila", "Aeromonas sobria", "Aeromonas veronii", 
  "Aggregatibacter actinomycetemcomitans", "Aggregatibacter aphrophilus", 
  "Arcanobacterium haemolyticum", "Bacillus anthracis", "Bacillus cereus", 
  "Bacillus thuringiensis", "Bacteroides fragilis", "Bartonella henselae", 
  "Bartonella quintana", "Bordetella bronchiseptica", "Bordetella hinzii", 
  "Bordetella holmesii", "Bordetella parapertussis", "Bordetella pertussis", 
  "Bordetella petrii", "Brucella abortus", "Brucella canis", "Brucella melitensis", 
  "Brucella suis", "Burkholderia cepacia complex", "Burkholderia gladioli", 
  "Burkholderia glumae", "Burkholderia mallei", "Burkholderia pseudomallei", 
  "Burkholderia thailandensis", "Campylobacter concisus", "Capnocytophaga gingivalis", 
  "Capnocytophaga leadbetteri", "Capnocytophaga sputigena", "Cardiobacterium hominis", 
  "Cardiobacterium valvarum", "Chlamydia pneumoniae", "Chlamydia psittaci", 
  "Chlamydia trachomatis", "Chromobacterium violaceum", "Citrobacter freundii complex", 
  "Citrobacter koseri", "Corynebacterium diphtheriae", "Corynebacterium jeikeium", 
  "Corynebacterium propinquum", "Corynebacterium pseudodiphtheriticum", 
  "Corynebacterium pseudotuberculosis", "Corynebacterium striatum", 
  "Corynebacterium ulcerans", "Coxiella burnetii", "Cronobacter sakazakii", 
  "Delftia acidovorans", "Dialister pneumosintes", "Dolosigranulum pigrum", 
  "Eikenella corrodens", "Elizabethkingia anophelis", "Elizabethkingia meningoseptica", 
  "Enterobacter cloacae complex", "Enterococcus faecalis", "Enterococcus faecium", 
  "Escherichia coli", "Eubacterium brachy", "Eubacterium limosum", 
  "Eubacterium nodatum", "Finegoldia magna (Peptostreptococcus magnus)", 
  "Francisella tularensis", "Fusobacterium necrophorum", "Fusobacterium nucleatum", 
  "Gemella haemolysans", "Gemella morbillorum", "Gordonia araii", "Gordonia bronchialis", 
  "Haemophilus haemolyticus", "Haemophilus influenzae", "Haemophilus parahaemolyticus", 
  "Haemophilus parainfluenzae", "Haemophilus pittmaniae", "Hafnia alvei", 
  "Kingella kingae", "Klebsiella aerogenes (Enterobacter aerogenes)", 
  "Klebsiella oxytoca", "Klebsiella pneumoniae", "Klebsiella quasipneumoniae", 
  "Klebsiella variicola", "Kytococcus sedentarius", "Leclercia adecarboxylata", 
  "Legionella anisa", "Legionella feeleii", "Legionella longbeachae", 
  "Legionella maceachernii", "Legionella pneumophila", "Legionella wadsworthii", 
  "Leptospira interrogans", "Leptotrichia buccalis", "Listeria monocytogenes", 
  "Moraxella catarrhalis", "Moraxella osloensis", "Morganella morganii", 
  "Mycobacterium avium complex", "Mycobacterium fortuitum (Mycolicibacterium fortuitum)", 
  "Mycobacterium gordonae", "Mycobacterium kansasii", "Mycobacterium malmoense", 
  "Mycobacterium scrofulaceum", "Mycobacterium simiae complex", "Mycobacterium szulgai", 
  "Mycobacterium tuberculosis complex", "Mycobacterium xenopi", 
  "Mycobacteroides abscessus (Mycobacterium abscessus)", 
  "Mycobacteroides chelonae (Mycobacterium chelonae)", 
  "Mycobacteroides immunogenum (Mycobacterium immunogenum)", "Mycoplasma pneumoniae", 
  "Neisseria flavescens", "Neisseria lactamica", "Neisseria meningitidis", 
  "Neisseria mucosa", "Nocardia abscessus", "Nocardia arthritidis", 
  "Nocardia beijingensis", "Nocardia brasiliensis", "Nocardia cyriacigeorgica", 
  "Nocardia farcinica", "Nocardia nova", "Nocardia otitidiscaviarum", 
  "Nocardia transvalensis", "Nocardia veterana", "Ochrobactrum anthropi", 
  "Orientia tsutsugamushi", "Pandoraea pulmonicola", "Pantoea agglomerans", 
  "Parvimonas micra", "Pasteurella multocida", "Pediococcus acidilactici", 
  "Peptostreptococcus anaerobius", "Prevotella buccae", "Prevotella intermedia", 
  "Prevotella melaninogenica", "Prevotella pleuritidis", "Proteus mirabilis", 
  "Proteus penneri", "Proteus vulgaris", "Providencia stuartii", 
  "Pseudomonas aeruginosa", "Pseudomonas fluorescens", "Pseudomonas stutzeri", 
  "Ralstonia pickettii", "Raoultella ornithinolytica", "Raoultella planticola", 
  "Rhodococcus hoagii (Rhodococcus equi)", "Rickettsia rickettsii", 
  "Rothia mucilaginosa", "Salmonella enterica", "Serratia marcescens", 
  "Shewanella putrefaciens", "Slackia exigua", "Sphingomonas paucimobilis", 
  "Staphylococcus aureus", "Stenotrophomonas maltophilia", "Streptococcus agalactiae", 
  "Streptococcus anginosus", "Streptococcus constellatus", "Streptococcus dysgalactiae", 
  "Streptococcus intermedius", "Streptococcus mitis", "Streptococcus pneumoniae", 
  "Streptococcus pyogenes", "Tatlockia micdadei (Legionella micdadei)", 
  "Treponema denticola", "Tropheryma whipplei", "Tsukamurella pulmonis", 
  "Tsukamurella tyrosinosolvens", "Ureaplasma parvum", "Ureaplasma urealyticum", 
  "Veillonella parvula", "Williamsia muralis", "Yersinia enterocolitica", 
  "Yersinia pestis", "Alternaria alternata", "Alternaria infectoria", 
  "Apophysomyces elegans", "Aspergillus flavus", "Aspergillus fumigatus", 
  "Aspergillus nidulans", "Aspergillus niger", "Aspergillus terreus", 
  "Aspergillus versicolor", "Blastomyces dermatitidis", 
  "Byssochlamys spectabilis (Paecilomyces variotii)", "Candida auris (Clavispora auris)", 
  "Cladophialophora bantiana", "Coccidioides immitis", "Coccidioides posadasii", 
  "Cryptococcus neoformans/Cryptococcus gattii", "Cunninghamella bertholletiae", 
  "Curvularia geniculata", "Curvularia lunata", "Exophiala dermatitidis", 
  "Fusarium oxysporum", "Fusarium proliferatum", "Fusarium solani", 
  "Fusarium verticillioides", "Histoplasma capsulatum", "Lichtheimia corymbifera", 
  "Lichtheimia ramosa", "Lomentospora prolificans (Scedosporium prolificans)", 
  "Microascus cinereus (Scopulariopsis cinereus)", "Microascus cirrosus (Scopulariopsis paisii)", 
  "Microascus pairsii (Scopulariopsis brumptii)", "Mucor circinelloides", 
  "Mucor indicus", "Mucor racemosus", "Paracoccidioides brasiliensis", 
  "Pneumocystis jirovecii", "Purpureocillium lilacinum (Paecilomyces lilacinus)", 
  "Rasamsonia aegroticola", "Rasamsonia argillacea", "Rhizomucor pusillus", 
  "Rhizopus azygosporus", "Rhizopus microsporus", "Rhizopus oryzae (Rhizopus arrhizus)", 
  "Saksenaea vasiformis", "Sarocladium kiliense", "Scedosporium apiospermum", 
  "Schizophyllum commune", "Scopulariopsis brevicaulis", "Sporothrix schenckii", 
  "Syncephalastrum racemosum", "Talaromyces marneffei", "Trichosporon asahii", 
  "Coxsackievirus A", "Coxsackievirus B", "Cytomegalovirus (CMV)", 
  "Enterovirus A71", "Enterovirus D68", "Epstein-Barr virus (EBV)", 
  "Herpes simplex virus 1 (HSV-1)", "Human adenovirus B", "Human adenovirus C", 
  "Human adenovirus E", "Human bocavirus 1", "Human coronavirus 229E", 
  "Human coronavirus HKU1", "Human coronavirus NL63", "Human coronavirus OC43", 
  "Human herpesvirus 6 (HHV-6)", "Human metapneumovirus", "Human parainfluenza virus 1", 
  "Human parainfluenza virus 2", "Human parainfluenza virus 3", 
  "Human parainfluenza virus 4", "Human parechovirus", "Human rhinovirus A", 
  "Human rhinovirus B", "Human rhinovirus C", "Influenza A virus (H1N1)", 
  "Influenza A virus (H3N2)", "Influenza A virus (H5N1)", "Influenza A virus (H7N9)", 
  "Influenza A virus (H9N2)", "Influenza B virus", "Influenza C virus", 
  "MERS coronavirus (MERS-CoV)", "Measles virus", "Mumps virus", 
  "Parvovirus B19", "Respiratory syncytial virus A", "Respiratory syncytial virus B", 
  "Rubella virus", "SARS coronavirus", "SARS-CoV-2 (2019-nCoV)", 
  "Varicella-zoster virus (VZV)"
)

rpip_pathogen_species_list_kraken_nomenclature_old_version <- c(
  # --- Bacteria ---
  "Achromobacter denitrificans", "Achromobacter xylosoxidans", "Acinetobacter baumannii", 
  "Acinetobacter lwoffii", "Acinetobacter nosocomialis", "Acinetobacter pittii", 
  "Actinomyces graevenitzii", "Actinomyces israelii", "Actinomyces meyeri", 
  "Actinomyces naeslundii", "Actinomyces odontolyticus", "Aeromonas caviae", 
  "Aeromonas hydrophila", "Aeromonas sobria", "Aeromonas veronii", 
  "Aggregatibacter actinomycetemcomitans", "Aggregatibacter aphrophilus", 
  "Arcanobacterium haemolyticum", "Bacillus anthracis", "Bacillus cereus", 
  "Bacillus thuringiensis", "Bacteroides fragilis", "Bartonella henselae", 
  "Bartonella quintana", "Bordetella bronchiseptica", "Bordetella hinzii", 
  "Bordetella holmesii", "Bordetella parapertussis", "Bordetella pertussis", 
  "Bordetella petrii", "Brucella abortus", "Brucella canis", "Brucella melitensis", 
  "Brucella suis", "Burkholderia cepacia complex", "Burkholderia gladioli", 
  "Burkholderia glumae", "Burkholderia mallei", "Burkholderia pseudomallei", 
  "Burkholderia thailandensis", "Campylobacter concisus", "Capnocytophaga gingivalis", 
  "Capnocytophaga leadbetteri", "Capnocytophaga sputigena", "Cardiobacterium hominis", 
  "Cardiobacterium valvarum", "Chlamydia pneumoniae", "Chlamydia psittaci", 
  "Chlamydia trachomatis", "Chromobacterium violaceum", "Citrobacter freundii complex", 
  "Citrobacter koseri", "Corynebacterium diphtheriae", "Corynebacterium jeikeium", 
  "Corynebacterium propinquum", "Corynebacterium pseudodiphtheriticum", 
  "Corynebacterium pseudotuberculosis", "Corynebacterium striatum", 
  "Corynebacterium ulcerans", "Coxiella burnetii", "Cronobacter sakazakii", 
  "Delftia acidovorans", "Dialister pneumosintes", "Dolosigranulum pigrum", 
  "Eikenella corrodens", "Elizabethkingia anophelis", "Elizabethkingia meningoseptica", 
  "Enterobacter cloacae complex", "Enterococcus faecalis", "Enterococcus faecium", 
  "Escherichia coli", "Eubacterium brachy", "Eubacterium limosum", 
  "Eubacterium nodatum", "Finegoldia magna", "Francisella tularensis", 
  "Fusobacterium necrophorum", "Fusobacterium nucleatum", "Gemella haemolysans", 
  "Gemella morbillorum", "Gordonia araii", "Gordonia bronchialis", 
  "Haemophilus haemolyticus", "Haemophilus influenzae", "Haemophilus parahaemolyticus", 
  "Haemophilus parainfluenzae", "Haemophilus pittmaniae", "Hafnia alvei", 
  "Kingella kingae", "Klebsiella aerogenes", "Klebsiella oxytoca", 
  "Klebsiella pneumoniae", "Klebsiella quasipneumoniae", "Klebsiella variicola", 
  "Kytococcus sedentarius", "Leclercia adecarboxylata", "Legionella anisa", 
  "Legionella feeleii", "Legionella longbeachae", "Legionella maceachernii", 
  "Legionella pneumophila", "Legionella wadsworthii", "Leptospira interrogans", 
  "Leptotrichia buccalis", "Listeria monocytogenes", "Moraxella catarrhalis", 
  "Moraxella osloensis", "Morganella morganii", "Mycobacterium avium complex", 
  "Mycolicibacterium fortuitum", "Mycobacterium gordonae", "Mycobacterium kansasii", 
  "Mycobacterium malmoense", "Mycobacterium scrofulaceum", "Mycobacterium simiae complex", 
  "Mycobacterium szulgai", "Mycobacterium tuberculosis complex", "Mycobacterium xenopi", 
  "Mycobacteroides abscessus", "Mycobacteroides chelonae", "Mycobacteroides immunogenum", 
  "Mycoplasma pneumoniae", "Neisseria flavescens", "Neisseria lactamica", 
  "Neisseria meningitidis", "Neisseria mucosa", "Nocardia abscessus", 
  "Nocardia arthritidis", "Nocardia beijingensis", "Nocardia brasiliensis", 
  "Nocardia cyriacigeorgica", "Nocardia farcinica", "Nocardia nova", 
  "Nocardia otitidiscaviarum", "Nocardia transvalensis", "Nocardia veterana", 
  "Mycolicibacterium anthropi", "Orientia tsutsugamushi", "Pandoraea pulmonicola", 
  "Pantoea agglomerans", "Parvimonas micra", "Pasteurella multocida", 
  "Pediococcus acidilactici", "Peptostreptococcus anaerobius", "Prevotella buccae", 
  "Prevotella intermedia", "Prevotella melaninogenica", "Prevotella pleuritidis", 
  "Proteus mirabilis", "Proteus penneri", "Proteus vulgaris", 
  "Providencia stuartii", "Pseudomonas aeruginosa", "Pseudomonas fluorescens", 
  "Pseudomonas stutzeri", "Ralstonia pickettii", "Raoultella ornithinolytica", 
  "Raoultella planticola", "Mycolicibacterium hoagii", "Rickettsia rickettsii", 
  "Rothia mucilaginosa", "Salmonella enterica", "Serratia marcescens", 
  "Shewanella putrefaciens", "Slackia exigua", "Sphingomonas paucimobilis", 
  "Staphylococcus aureus", "Stenotrophomonas maltophilia", "Streptococcus agalactiae", 
  "Streptococcus anginosus", "Streptococcus constellatus", "Streptococcus dysgalactiae", 
  "Streptococcus intermedius", "Streptococcus mitis", "Streptococcus pneumoniae", 
  "Streptococcus pyogenes", "Legionella micdadei", "Treponema denticola", 
  "Tropheryma whipplei", "Tsukamurella pulmonis", "Tsukamurella tyrosinosolvens", 
  "Ureaplasma parvum", "Ureaplasma urealyticum", "Veillonella parvula", 
  "Williamsia muralis", "Yersinia enterocolitica", "Yersinia pestis",

  # --- Fungi ---
  "Alternaria alternata", "Alternaria infectoria", "Apophysomyces elegans", 
  "Aspergillus flavus", "Aspergillus fumigatus", "Aspergillus nidulans", 
  "Aspergillus niger", "Aspergillus terreus", "Aspergillus versicolor", 
  "Blastomyces dermatitidis", "Paecilomyces variotii", "Candida auris", 
  "Cladophialophora bantiana", "Coccidioides immitis", "Coccidioides posadasii", 
  "Cryptococcus neoformans", "Cryptococcus gattii", "Cunninghamella bertholletiae", 
  "Curvularia geniculata", "Curvularia lunata", "Exophiala dermatitidis", 
  "Fusarium oxysporum", "Fusarium proliferatum", "Fusarium solani", 
  "Fusarium verticillioides", "Histoplasma capsulatum", "Lichtheimia corymbifera", 
  "Lichtheimia ramosa", "Lomentospora prolificans", "Microascus cinereus", 
  "Microascus cirrosus", "Microascus pairsii", "Mucor circinelloides", 
  "Mucor indicus", "Mucor racemosus", "Paracoccidioides brasiliensis", 
  "Pneumocystis jirovecii", "Purpureocillium lilacinum", "Rasamsonia aegroticola", 
  "Rasamsonia argillacea", "Rhizomucor pusillus", "Rhizopus azygosporus", 
  "Rhizopus microsporus", "Rhizopus oryzae", "Saksenaea vasiformis", 
  "Sarocladium kiliense", "Scedosporium apiospermum", "Schizophyllum commune", 
  "Scopulariopsis brevicaulis", "Sporothrix schenckii", "Syncephalastrum racemosum", 
  "Talaromyces marneffei", "Trichosporon asahii",

  # --- Viruses ---
  "Coxsackievirus A", "Coxsackievirus B", "Human betaherpesvirus 5", # CMV
  "Enterovirus A71", "Enterovirus D68", "Human gammaherpesvirus 4", # EBV
  "Human alphaherpesvirus 1", # HSV-1
  "Human adenovirus B", "Human adenovirus C", "Human adenovirus E", 
  "Human bocavirus 1", "Human coronavirus 229E", "Human coronavirus HKU1", 
  "Human coronavirus NL63", "Human coronavirus OC43", "Human betaherpesvirus 6", # HHV-6
  "Human metapneumovirus", "Human parainfluenza virus 1", "Human parainfluenza virus 2", 
  "Human parainfluenza virus 3", "Human parainfluenza virus 4", "Human parechovirus", 
  "Human rhinovirus A", "Human rhinovirus B", "Human rhinovirus C", 
  "Influenza A virus", "Influenza B virus", "Influenza C virus", 
  "Middle East respiratory syndrome-related coronavirus", "Measles morbillivirus", 
  "Mumps orthorubulavirus", "Primate erythroparvovirus 1", # B19
  "Human orthopneumovirus", # RSV
  "Rubella virus", "Severe acute respiratory syndrome-related coronavirus", # SARS
  "Severe acute respiratory syndrome coronavirus 2", # SARS-CoV-2
  "Human alphaherpesvirus 3" # VZV
)

#print differences between old and new lists:
# Clean function to remove parentheses and extra spaces for better matching
clean_list <- function(x) {
  x <- gsub(" \\(.*\\)", "", x) # Remove everything in parentheses
  x <- gsub("/.*", "", x)       # Handle the Cryptococcus slash
  trimws(tolower(x))            # Lowercase and trim whitespace
}

# Find differences
only_in_A <- rpip_pathogen_species_list_illumina_terms_new_version[!(clean_list(rpip_pathogen_species_list_illumina_terms_new_version) %in% clean_list(rpip_pathogen_species_list_illumina_terms_old_version))]
only_in_B <- rpip_pathogen_species_list_illumina_terms_old_version[!(clean_list(rpip_pathogen_species_list_illumina_terms_old_version) %in% clean_list(rpip_pathogen_species_list_illumina_terms_new_version))]

print(only_in_A)
print(only_in_B)

#create pathogen species list
rpip_pathogen_species_list= rpip_pathogen_species_list_kraken_nomenclature_old_version


#taxa_detected = kraken_merged_df %>% dplyr::distinct(name) %>% pull(name)

# taxa_detected[1:999]
# taxa_detected[1000:1999]
# taxa_detected[2000:2900]

detected_pathogens =# R-compatible vector of pathogen species from the provided taxa
 c(
  "Klebsiella michiganensis",
  "Klebsiella quasipneumoniae",
  "Klebsiella pneumoniae",
  "Klebsiella variicola",
  "Klebsiella oxytoca",
  "Raoultella ornithinolytica",
  "Raoultella planticola",
  "Citrobacter freundii",
  "Citrobacter farmeri",
  "Enterobacter asburiae",
  "Enterobacter roggenkampii",
  "Enterobacter cloacae",
  "Enterobacter kobei",
  "Enterobacter hormaechei",
  "Escherichia coli",
  "Salmonella enterica",
  "Serratia marcescens",
  "Aeromonas hydrophila",
  "Aeromonas veronii",
  "Aeromonas salmonicida",
  "Aeromonas dhakensis",
  "Aeromonas sobria",
  "Pseudomonas aeruginosa",
  "Acinetobacter baumannii",
  "Acinetobacter nosocomialis",
  "Acinetobacter junii",
  "Stenotrophomonas maltophilia",
  "Haemophilus parainfluenzae",
  "Laribacter hongkongensis",
  "Bacteroides fragilis",
  "Elizabethkingia anophelis",
  "Elizabethkingia bruuniana",
  "Myroides odoratimimus",
  "Prevotella melaninogenica",
  "Clostridioides difficile",
  "Enterococcus faecium",
  "Enterococcus cecorum",
  "Streptococcus pneumoniae",
  "Streptococcus anginosus",
  "Streptococcus mitis",
  "Finegoldia magna",
  "Naegleria fowleri",

  "Acinetobacter pittii",
  "Acinetobacter lwoffii",
  "Acinetobacter seifertii",
  "Morganella morganii",
  "Proteus mirabilis",
  "Brucella anthropi",
  "Aliarcobacter butzleri",
  "Campylobacter coli",
  "Campylobacter rectus",
  "Streptococcus agalactiae",
  "Streptococcus intermedius",
  "Listeria monocytogenes",
  "Burkholderia mallei",
  "Burkholderia pseudomallei",
  "Burkholderia cepacia",
  "Mycolicibacterium fortuitum",
  "Tsukamurella tyrosinosolvens",
  "Riemerella anatipestifer",

  "Fusobacterium vincentii",
  "Ralstonia mannitolilytica",
  "Arcobacter porcinus",
  "Acinetobacter haemolyticus",
  "Aeromonas taiwanensis",
  "Aeromonas jandaei",
  "Dysgonomonas capnocytophagoides",
  "Clostridium perfringens",
  "Enterobacter ludwigii",
  "Leclercia adecarboxylata",
  "Elizabethkingia miricola",
  "Bacillus cereus",
  "Klebsiella aerogenes",
  "Bacteroides difficilis",
  "Gardnerella sp. DNF00536",
  "Bordetella bronchiseptica",
  "Bordetella pertussis",
  "Bordetella parapertussis",
  "Bordetella trematum",
  "Achromobacter xylosoxidans",
  "Burkholderia cenocepacia",
  "Sneathia vaginalis",
  "Rothia dentocariosa",
  "Fusarium verticillioides",
  "Comamonas kerstersii",
  "Sphingobacterium multivorum",
  "Abiotrophia defectiva",
  "Schaalia odontolytica",
  "Ralstonia pickettii",
  "Actinomyces israelii",
  "Lawsonella clevelandensis",
  "Sphingomonas paucimobilis",
  "Fusobacterium gonidiaformans",
  "Empedobacter falsenii",
  "Sphingobacterium spiritivorum",
  "Pandoraea pnomenusa",
  "Chryseobacterium indologenes",
  "Granulicatella adiacens",
  "Gemella haemolysans",
  "Prevotella intermedia",
  "Corynebacterium jeikeium",
  "Proteus vulgaris",
  "Burkholderia contaminans",
  "Vibrio cholerae",
  "Escherichia albertii",

  "Mycobacteroides abscessus",
  "Dialister pneumosintes",
  "Shewanella putrefaciens",
  "Acinetobacter colistiniresistens",
  "Klebsiella africana",
  "Citrobacter braakii",
  "Pseudomonas ceruminis",
  "Pseudomonas asiatica",
  "Fusobacterium mortiferum",
  "Providencia rettgeri",
  "Shigella flexneri",
  "Dermatobacter hominis",
  "Citrobacter koseri",
  "Desulfovibrio fairfieldensis",
  "Citrobacter amalonaticus",
  "Hafnia paralvei",
  "Burkholderia multivorans",
  "Mycobacteroides immunogenum",
  "Kluyvera ascorbata",
  "Haemophilus pittmaniae",
  "Campylobacter ureolyticus",
  "Burkholderia vietnamiensis",
  "Clostridium butyricum",
  "Enterococcus casseliflavus",
  "Selenomonas sputigena",
  "Actinomyces naeslundii",
  "Burkholderia aenigmatica",
  "Streptococcus pyogenes",
  "Prevotella nigrescens",
  "Prevotella dentalis",
  "Porphyromonas asaccharolytica",
  "Fusobacterium animalis",
  "Sneathia sanguinegens",
  "Citrobacter werkmanii",
  "Veillonella dispar",
  "Corynebacterium striatum",
  "Kerstersia gyiorum",
  "Atlantibacter hermannii",
  "Filifactor alocis",
  "Acinetobacter septicus",
  "Streptococcus suis",

  "Rhodococcoides fascians",         # Plant pathogen, emerging human clinical isolate
  "Fusarium musae",                  # Fungal pathogen (keratitis/clinical)
  "Acinetobacter gandensis",         # Opportunistic pathogen
  "Brevundimonas vancanneytii",      # Rare clinical isolate
  "Bacteroides finegoldii",          # Clinical isolate (abscesses/blood)
  "Ralstonia holmesii",              # Opportunistic (septicemia)
  "Prevotella denticola",            # Oral/periodontal pathogen
  "Acinetobacter portensis",         # Emerging clinical species
  "Pseudomonas fluorescens",         # Opportunistic/transfusion-associated
  "Vandammella animalimorsus",       # Emerging (animal bite wounds)
  "Tsukamurella pulmonis",           # Lung/systemic infections
  "Anaerofustis stercorihominis",    # Human clinical isolate
  "Capnocytophaga leadbetteri",      # Oral pathogen
  "Hoylesella buccalis",             # Oral/periodontal pathogen
  "Aerococcus mictus",               # Emerging urinary/clinical pathogen
  "Legionella lytica",               # Atypical pneumonia (Legionellosis)
  "Enterobacter wuhouensis",         # Emerging clinical isolate
  "Achromobacter ruhlandii",         # Opportunistic (Cystic Fibrosis)
  "Metamycoplasma hominis",          # Urogenital pathogen
  "Chlamydia",                       # Genus-level (Major pathogen group)
  "Parachlamydia acanthamoebae",     # Emerging respiratory pathogen
  "Streptobacillus moniliformis",    # Rat-bite fever
  "Paracoccus yeei",                 # Rare opportunistic pathogen
  "Mycobacterium avium",             # Opportunistic (MAC complex)
  "Mycolicibacillus koreensis",      # Emerging mycobacterial pathogen
  "Dietzia timorensis",              # Rare clinical isolate
  "Brevibacterium casei",            # Opportunistic (catheter-related)
  "Cellulomonas hominis",            # Human clinical isolate
  "Cellulosimicrobium cellulans",    # Opportunistic pathogen
  "Fusarium keratoplasticum",        # Major fungal pathogen (keratitis)
  "Prevotella corporis",             # Anaerobic infection isolate
  "Prevotella bivia",                # OB/GYN associated pathogen
  "Desulfobulbus oralis",            # Periodontal disease
  "Gordonia bronchialis",            # Respiratory/sternal infections
  "Clostridium baratii",             # Botulism-associated (Type F)
  "Saccharomyces cerevisiae",        # Opportunistic (Fungemia)
  "Citrobacter gillenii",            # Opportunistic pathogen
  "Corynebacterium diphtheriae",     # Diphtheria (Major pathogen)
  "Cronobacter sakazakii",           # Infant meningitis/enterocolitis
  "Raoultella terrigena",            # Opportunistic pathogen
  "Rothia kristinae",                # Endocarditis/opportunistic
  "Streptococcus constellatus",      # Abscess/endocarditis
  "Winkia neuii",                    # Emerging clinical isolate
  "Staphylococcus aureus",           # Major pathogen
  "Gordonia otitidis",               # Ear/clinical infections
  "Mycobacterium intracellulare",    # Opportunistic (MAC complex)
  "Mycolicibacterium neoaurum",      # Catheter-related bacteremia
  "Corynebacterium xerosis",         # Opportunistic pathogen
  "Murdochiella vaginalis",          # Vaginal/clinical isolate
  "Shigella dysenteriae",            # Dysentery (Major pathogen)
  "Pauljensenia hongkongensis",      # Emerging anaerobic pathogen
  "Roseomonas gilardii",             # Opportunistic/bacteremia
  "Roseomonas mucosa",               # Opportunistic clinical isolate
  "Chryseobacterium gleum",           # Emerging opportunistic pathogen



  "Mycobacterium kansasii",
  "Mycobacterium kubicae",
  "Mycobacterium tuberculosis",
  "Mycolicibacterium mucogenicum",
  "Simkania negevensis",
  "Brucella pseudogrignonensis",
  "Mycobacterium marinum",
  "Brucella pseudintermedia",
  "Sphingomonas sanguinis",
  "Sphingomonas parapaucimobilis",
  "Stenotrophomonas sp. BIO128-Bstrain",
  "Fusarium oxysporum",
  "Mycobacterium paraintracellulare",
  "Chryseobacterium gambrini",
  "Providencia sp. PROV162",
  "Mycoplasma sp. P36-A1",
  "Mycobacterium noviomagense",
  "Intestinibacter bartlettii",
  "Haematobacter massiliensis",
  "Brevundimonas vesicularis",
  "Sphingobacterium cellulitidis",
  "Legionella pneumophila",
  "Corynebacterium durum",
  "Actinomyces oris",
  "Actinomyces massiliensis",
  "Slackia exigua",
  "Fretibacterium fastidiosum",
  "Fusarium falciforme",
  "Nakaseomyces glabratus",
  "Mycobacterium riyadhense",
  "Mycobacterium canetti",
  "Parascardovia denticolens",
  "Fusarium poae",
  "Microbacterium schleiferi",
  "Dermacoccus nishinomiyaensis",
  "Lautropia mirabilis",
  "Brucella pituitosa",
  "Anaerostipes caccae",
  "Umbribacter vaginalis",
  "Microbacterium arborescens",
  "Empedobacter stercoris",
  "Janibacter terrae"
)



#---------Look at differences between total detected species and RPIP pathogens by detection frequency------------------
#get unique list of total detected at species level
detected_species = kraken_merged_df %>% filter(taxonomy_lvl== "S") %>% distinct(name) %>% pull(name)


#
detection_frequency_df = kraken_merged_df %>% 
  filter(name %in% c(detected_species,rpip_pathogen_species_list)) %>% 
  mutate(rel_ab = as.numeric(fraction_total_reads)) %>% 
  filter(rel_ab>0) %>% 
  distinct(SampleID, name,rel_ab) %>% 
  
  # 1. Expand to include zeros for every SampleID/Species combination
  # Note: If some SampleIDs are missing entirely, use: 
  # complete(SampleID = unique(bracken_merged_df$SampleID), name, fill = list(rel_ab = 0))
  #complete(SampleID, name, fill = list(rel_ab = 0)) %>%
  complete(SampleID = unique(kraken_merged_df$SampleID), name, fill = list(rel_ab = 0)) %>% 
  
  #calculate detection frequency across all samples
  group_by(name) %>% 
  mutate(total = n(),
         total_positive = sum(rel_ab>0),
         perc = (total_positive/total) * 100) %>% 
  ungroup() %>% 
  
  mutate(Type = case_when(name %in% rpip_pathogen_species_list ~ "RPIP",
                          .default = "non-RPIP Species"))
  

plot = detection_frequency_df %>% 
  ggplot(aes(x = Type, y = perc)) + 
  geom_boxplot()+
  
  theme_bw()+
  
  labs(x = "",y= "Detection Frequency (%)")

print(plot)

#---------Look at differences in RPIP and non RPIP pathogens by detection frequency-------

detection_frequency_df = kraken_merged_df %>% 
  filter(name %in% c(detected_pathogens,rpip_pathogen_species_list)) %>% 
  mutate(rel_ab = as.numeric(fraction_total_reads)) %>% 
  filter(rel_ab>0) %>% 
  distinct(SampleID, name,rel_ab) %>% 
  
  # 1. Expand to include zeros for every SampleID/Species combination
  # Note: If some SampleIDs are missing entirely, use: 
  # complete(SampleID = unique(bracken_merged_df$SampleID), name, fill = list(rel_ab = 0))
  complete(SampleID, name, fill = list(rel_ab = 0)) %>%
  
  #calculate detection frequency across all samples
  group_by(name) %>% 
  mutate(total = n(),
         total_positive = sum(rel_ab>0),
         perc = (total_positive/total) * 100) %>% 
  ungroup() %>% 
  
  mutate(Type = case_when(name %in% rpip_pathogen_species_list ~ "RPIP",
                          .default = "non-RPIP Pathogens"))
  

plot = detection_frequency_df %>% 
  ggplot(aes(x = Type, y = perc)) + 
  geom_boxplot()+
  
  theme_bw()+
  
  labs(x = "",y= "Detection Frequency (%)")

print(plot)

#---------Look at differences in RPIP and non RPIP pathogens by rel abundance-------

mean_rel_ab_df = kraken_merged_df %>% 
  filter(name %in% c(detected_pathogens,rpip_pathogen_species_list)) %>% 
  mutate(rel_ab = as.numeric(fraction_total_reads)) %>% 
  filter(rel_ab>0) %>% 
  distinct(SampleID, name,rel_ab) %>% 
  
  mutate(Type = case_when(name %in% rpip_pathogen_species_list ~ "RPIP",
                          .default = "non-RPIP Pathogens")) %>% 
  
  
  #calculate detection frequency across all samples
  group_by(name,Type) %>% 
  mutate(mean_rel_ab = mean(rel_ab,na.rm = TRUE)) %>% 
  ungroup() 
  

plot = mean_rel_ab_df %>% 
  ggplot(aes(x = name, y = mean_rel_ab)) + 
  # Set linewidth to make the boxplot outline thin
  # Set outlier.shape = NA so we don't double-count points
  geom_boxplot(linewidth = 0.2, outlier.shape = NA) +
  
  # Add jittered points: 
  # size adjusts the dot size, alpha makes them slightly transparent
  # width controls how much they spread out horizontally
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.4, color = "darkgrey") +
  
  facet_grid(Type ~ ., scales = "free") +
  coord_flip() +
  theme_classic() +
  labs(x = "", y = "Mean Relative Abundance") +
  theme(
    axis.text.y = element_text(size = 5), # Bumped up slightly, 2 is microscopic!
    strip.text.y = element_text(angle = 0) # Makes facet labels easier to read
  )

print(plot)

#---------Look at differences in RPIP and non RPIP pathogens by rel abundance-------

mean_rel_ab_df = kraken_merged_df %>% 
  filter(name %in% c(detected_pathogens,rpip_pathogen_species_list)) %>% 
  mutate(rel_ab = as.numeric(fraction_total_reads)) %>% 
  filter(rel_ab>0) %>% 
  distinct(SampleID, name,rel_ab) %>% 
  
  mutate(Type = case_when(name %in% rpip_pathogen_species_list ~ "RPIP",
                          .default = "non-RPIP Pathogens")) %>% 
  
  
  #calculate detection frequency across all samples
  group_by(name,Type) %>% 
  mutate(mean_rel_ab = mean(rel_ab,na.rm = TRUE)) %>% 
  ungroup() %>% 
  
  filter(mean_rel_ab > 0.0001)
  

plot = mean_rel_ab_df %>% 
  ggplot(aes(x = name, y = rel_ab)) + 
  # Set linewidth to make the boxplot outline thin
  # Set outlier.shape = NA so we don't double-count points
  geom_boxplot(linewidth = 0.2, outlier.shape = NA) +
  
  # Add jittered points: 
  # size adjusts the dot size, alpha makes them slightly transparent
  # width controls how much they spread out horizontally
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.4, color = "darkgrey") +
  
  facet_grid(Type ~ ., scales = "free") +
  coord_flip() +
  theme_classic() +
  labs(x = "", y = "Relative Abundance") +
  theme(
    axis.text.y = element_text(size = 5), # Bumped up slightly, 2 is microscopic!
    strip.text.y = element_text(angle = 0) # Makes facet labels easier to read
  )

print(plot)

#------Sum rel ab of all RPIP pathogens and multiply total reads per sample and that would be total reads for (stacked bar plot)------------

#--------Stacked bar plot of assinged reads to RPIP vs not------------
stacked_bar_df = kraken_merged_df %>% 
  filter(name %in% c(detected_species,rpip_pathogen_species_list)) %>% 
  mutate(rel_ab = as.numeric(fraction_total_reads)) %>% 
  distinct(SampleID, name,rel_ab) %>% 
  
  mutate(Type = case_when(name %in% rpip_pathogen_species_list ~ "RPIP",
                          .default = "non-RPIP Species")) %>% 
  
  
  mutate(total = sum(rel_ab,na.rm =T )) %>% 
  
  group_by(Type) %>% 
  mutate(frac = sum(rel_ab,na.rm = TRUE)) %>% 
  ungroup() %>% 
  
  mutate(perc = (frac/total)) %>% 
  
  distinct(Type,frac,perc) %>% 
  mutate(x_axis = "Value")
  
  
  
  


plot = stacked_bar_df %>% 
  ggplot(aes(x = x_axis, y = perc,fill = Type)) + 

  geom_col(position = "stack")+
  
  theme_classic() +
  labs(x = "", y = "Rel Abundance") +
  theme(
    axis.text.y = element_text(size = 5), # Bumped up slightly, 2 is microscopic!
    strip.text.y = element_text(angle = 0) # Makes facet labels easier to read
  )

print(plot)

# Save final heatmap
ggsave(
  plot = plot,
  filename = "Biofilm Project Figures/RPIP vs non RPIP reads.png",
  width = 3,
  height = 4,
  units = "in"
)

#------------------final pathogen list---------

#create final pathogen species list plus notable others: Clostridium perfrigens, Campy coli 
rpip_pathogen_species_list= c(rpip_pathogen_species_list_kraken_nomenclature_old_version)

rpip_pathogen_species_list_additional= c(rpip_pathogen_species_list_kraken_nomenclature_old_version,"Naegleria fowleri","Campylobacter coli","Clostridium perfrigens","Mycobacterium avium complex", detected_pathogens)

```

#Make Spreadsheet with filtered kraken output for only relevant pathogens for BLASTing
```{r}

spreadsheet_to_blast = kraken_rel_ab %>% 
  
  filter(name %in% rpip_pathogen_species_list)


#write_xlsx(spreadsheet_to_blast,"Output Spreadsheets/kraken_sample_taxid_list_fixed_with_0.0001_rel_ab_threshold.xlsx")

write_xlsx(kraken_merged_df,"Output Spreadsheets/kraken_sample_taxid_list_fixed_with_0.0001%_rel_ab_threshold.xlsx")

taxids_total = spreadsheet_to_blast %>% 
  dplyr::select(taxid) %>% distinct()
```

#REMOVE BLASTED FALSE POSITIVES
```{r}

kraken_merged_df_filtered  = kraken_merged_df %>% full_join(blast_report_df,by = c("SampleID", "taxonomy_id","name")) 
  
  # #REMOVE QUALITATIVE FALSE POSITIVES
  # filter(name != "Raoultella ornithinolytica", name != "Raoultella planticola")

  # drop_na(fraction_total_reads) %>%
  # 
  # filter(classification == "TRUE_POSITIVE"|
  #        classification == "NO_DATA" | is.na(classification) == TRUE)
  # 
  # # filter(classification!= "FALSE_POSITIVE",
  # #          classification != "UNCERTAIN")
  # 

library(writexl)
#write_xlsx(kraken_merged_df_filtered,"Input Data/kraken_ct_0_15_min_hit_2_filtered_out_less_than_0.005%relab_and_used_total_reads_assigned qcd_and_trimmed_and_dehosted_as_denominator.xlsx")
```




#Check Dates
```{r}

##Goal: go into 16s data and figure out which dates for each Sites werre included, then create a new column called tiempoint for metagenomic data and for each sample label it as timepoint 1a/b/c


#-----------------load 16s dataframe ------------
sixteens_df <- read_xlsx(
  "Input Data/16S_Metadata.xlsx"
) %>% 
  dplyr::select(sample_id,sample_type,sample_date,timepoint,corresponding_sewer,sub_timepoint) %>% 
  drop_na(sample_id) %>% drop_na(sample_date) %>% 
  mutate(
    sample_date = as.numeric(sample_date),
    sample_date = as.Date(sample_date, origin = "1899-12-30")) %>% 
  dplyr::select(sample_type,sample_date,timepoint,corresponding_sewer,sub_timepoint) %>%  
  
  mutate(sample_type = case_when(sample_type == "w" ~ "H_WW",
                                 sample_type == "e" ~ "Endcap",
                                 sample_type == "d" ~ "Drain",
                                 .default = sample_type)) %>% 
  filter(sample_type %in% c("H_WW", "Endcap", "Drain")) %>% 
  
  dplyr::rename("Type" = "sample_type",
         "Date" = "sample_date",
         "Location" = "corresponding_sewer")

#specify dates

#------------Don't Filter Dates --------------------------

# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame


# RPIP taxa
rpip_targets <-  rpip_pathogen_species_list


# ---------- STEP 1: Extract and clean ----------
heatmap_df <- kraken_merged_df_filtered %>%
  filter(taxonomy_lvl == "S") %>% 
  mutate(rel_ab = as.numeric(fraction_total_reads)) %>% 
  filter(Location %in% locations) %>% 
  filter(name %in% c(rpip_targets)) %>% 
  filter(Type %in% c("Drain", "Endcap", "H_WW")) %>% 

  
  mutate(
    Type = factor(Type, levels = c("Drain", "Endcap", "H_WW")),
    TypeLabel = case_when(
      Type == "H_WW" ~ "Hospital\nWastewater", 
      TRUE ~ as.character(Type)
    ),
    # Make Date a string in MM/DD/YYYY format
    DateString = format(as.Date(Date), "%m/%d/%Y")
  ) %>% 
  
  left_join(sixteens_df,by = c("Type","Location","Date"))



#-----------------Check Dates with Bar Chart ----------------


check_dates = heatmap_df %>% 
  dplyr::select(SampleID,Location, Type,Date,timepoint,sub_timepoint) %>% distinct()


check_dates_plot <- check_dates %>%
  mutate(Date = as.Date(Date)) %>%
  dplyr::count(Location, timepoint, Date, Type) %>%
  mutate(Date = factor(Date, levels = sort(unique(Date))))

ggplot(check_dates_plot,
       aes(x = Date, y = n, fill = Type)) +
  geom_col(width = 0.8) +
  # facet_grid(
  #   rows = vars(Location),
  #   cols = vars(timepoint),
  #   scales = "free_x",
  #   space  = "free_x"
  # ) +
  facet_grid(Location ~ timepoint, scales = "free_x", space = "free_x") +
  scale_fill_viridis_d(option = "D") +
  #scale_x_date(date_labels = "%m/%d/%Y") +
  labs(
    x = "Sampling date",
    y = "Number of samples",
    fill = "Sample type",
    title = "Samples collected per date by location and timepoint"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA),
    axis.title = element_text(face = "bold")
  )

#summarize # of samples for each type, etc
total_by_type_df = check_dates_plot %>% 
  group_by(Type) %>% 
  summarise(total_by_type = sum(n)) %>% 
  ungroup()

#-----------------Check Dates with Bar Chart and Manually Changed Timepoint for MGX data ----------------
check_dates = heatmap_df %>% 
  dplyr::select(SampleID,Location, Type,Date) %>% distinct()

check_dates_plot <- check_dates %>%
  mutate(Date = as.Date(Date)) %>%
  dplyr::count(Location, Date, Type) %>%
  
  mutate(
    timepoint = case_when(
      Date >= as.Date("2023-12-01") & Date <= as.Date("2024-04-30") ~ "Pre",
      Date >= as.Date("2024-08-01") & Date <= as.Date("2024-09-30") ~ "1",
      Date >= as.Date("2024-10-01") & Date <= as.Date("2024-11-30") ~ "2",
      Date >= as.Date("2024-12-01") & Date <= as.Date("2024-12-31") ~ "3",
      TRUE ~ NA_character_
    )
  ) %>%

  
  mutate(Type = case_when(Type == "Endcap" ~ "Sewer Biofilm",
                          Type == "H_WW" ~ "Wastewater",
                          Type == "Drain" ~ "Sink Biofilm")) %>% 
  
    # now replace Site names 
  mutate(Location = case_when(Location == "MARIO" ~ "Site 1",
                              Location == "FIONA" ~ "Site 2",
                              Location == "LUIGI" ~ "Site 3",
                              Location == "SHREK" ~ "Site 4",
                              Location == "OSCAR" ~ "Site 5")) %>% 
  
  mutate(SeqType = "MGX") %>% 
  
  dplyr::select(Location,Date,Type,timepoint,SeqType,n)

#combine with 16 s df
sixteens_df_for_combining  = read_xlsx(
  "Input Data/16S_Metadata.xlsx"
) %>% 
  dplyr::select(sample_id,sample_type,sample_date,timepoint,corresponding_sewer,sub_timepoint) %>% 
  drop_na(sample_id) %>% drop_na(sample_date) %>% 
  mutate(
    sample_date = as.numeric(sample_date),
    sample_date = as.Date(sample_date, origin = "1899-12-30")) %>% 
  dplyr::select(sample_type,sample_date,timepoint,corresponding_sewer,sub_timepoint) %>%  
  
  mutate(sample_type = case_when(sample_type == "w" ~ "H_WW",
                                 sample_type == "e" ~ "Endcap",
                                 sample_type == "d" ~ "Drain",
                                 .default = sample_type)) %>% 
  
  dplyr::rename("Type" = "sample_type",
         "Date" = "sample_date",
         "Location" = "corresponding_sewer") %>% 
  
  mutate(Date = as.Date(Date)) %>%
  
  filter(Location %in% locations) %>% 
  
  dplyr::count(Location, Date, Type,timepoint)  %>% 

  
  mutate(Type = case_when(Type == "Endcap" ~ "Sewer Biofilm",
                          Type == "H_WW" ~ "Wastewater",
                          Type == "Drain" ~ "Sink Biofilm",
                          Type == "tap" ~ "Tap",
                          .default = Type)) %>% 
  
  # now replace Site names 
  mutate(Location = case_when(Location == "MARIO" ~ "Site 1",
                              Location == "FIONA" ~ "Site 2",
                              Location == "LUIGI" ~ "Site 3",
                              Location == "SHREK" ~ "Site 4",
                              Location == "OSCAR" ~ "Site 5")) %>% 
  
  mutate(SeqType = "16S") %>% 
  
  dplyr::select(Location,Date,Type,timepoint,SeqType,n)

combined_df = bind_rows(check_dates_plot, sixteens_df_for_combining)


plot = ggplot(combined_df,
       aes(x = Date, y = n, fill = Type)) +
  geom_col(width = 0.8) +
  #facet_grid(rows = vars(Location),cols = vars(timepoint),scales = "free_x",space  = "free_x") +
  facet_grid(Location ~ timepoint + SeqType, scales = "free_x", space = "free_x") +
  scale_fill_viridis_d(option = "D") +
  #scale_x_date(date_labels = "%m/%d/%Y") +
  labs(
    x = "Sampling date",
    y = "Number of samples",
    fill = "Sample type",
    title = "Samples collected per date by location and timepoint"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA),
    axis.title = element_text(face = "bold")
  )

print(plot)

# Save final heatmap
ggsave(
  plot = plot,
  filename = "Biofilm Project Figures/Check Dates Shared with Biofilm By Time Point with 16S .png",
  width = 9,
  height = 6,
  units = "in"
)

#summarize # of samples for each type, etc
total_by_type_df = combined_df %>% 
  group_by(Type,SeqType) %>% 
  summarise(total_by_sample_type_andseqtype = sum(n)) %>% 
  ungroup()

total_by_type_df = combined_df %>% 
  group_by(Type,SeqType,timepoint) %>% 
  summarise(total_by_sample_type_andseqtype = sum(n)) %>% 
  ungroup()





#-------------Filter dates down to only samples with paired 16S data-------------
dates = c("2024-09-09", "2024-09-11", "2024-09-12", "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09")

# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame


# RPIP taxa
rpip_targets <-  rpip_pathogen_species_list

# ---------- STEP 1: Extract and clean ----------
# 
# twelve_nine = kraken_merged_df_filtered %>% 
#   
#   #left_join(sixteens_df,by = c("Type","Location","Date"))%>%
#   
#   dplyr::select(SampleID, Date, Type, Location) %>% distinct() 
# 
#   filter(taxonomy_lvl == "S") %>% 
#   mutate(rel_ab = as.numeric(fraction_total_reads)) %>% 
#   filter(Location %in% locations) %>% 
#   filter(name %in% c(rpip_targets)) %>% 
#   filter(Type %in% c("Drain", "Endcap", "H_WW")) %>% 
#   
#   #filter dates
#   filter(
#     (Type == "H_WW" & Date %in% dates) |
#     (Type != "H_WW")
#   ) %>% 
#   
#   mutate(
#     Type = factor(Type, levels = c("Drain", "Endcap", "H_WW")),
#     TypeLabel = case_when(
#       Type == "H_WW" ~ "Hospital\nWastewater", 
#       TRUE ~ as.character(Type)
#     ),
#     # Make Date a string in MM/DD/YYYY format
#     DateString = format(as.Date(Date), "%m/%d/%Y")
#   ) 

heatmap_df <- kraken_merged_df_filtered %>%
  filter(taxonomy_lvl == "S") %>% 
  mutate(rel_ab = as.numeric(fraction_total_reads)) %>% 
  filter(Location %in% locations) %>% 
  filter(name %in% c(rpip_targets)) %>% 
  filter(Type %in% c("Drain", "Endcap", "H_WW")) %>% 
  
  #filter dates
  filter(
    (Type == "H_WW" & Date %in% dates) |
    (Type != "H_WW")
  ) %>% 
  
  mutate(
    Type = factor(Type, levels = c("Drain", "Endcap", "H_WW")),
    TypeLabel = case_when(
      Type == "H_WW" ~ "Hospital\nWastewater", 
      TRUE ~ as.character(Type)
    ),
    # Make Date a string in MM/DD/YYYY format
    DateString = format(as.Date(Date), "%m/%d/%Y")
  ) %>% 
  
  left_join(sixteens_df,by = c("Type","Location","Date")) %>% 
  
  drop_na(timepoint)


#-----------------Check Dates with Bar Chart ----------------
check_dates = heatmap_df %>% 
  dplyr::select(SampleID,Location, Type,Date,timepoint,sub_timepoint) %>% distinct()


check_dates_plot <- check_dates %>%
  mutate(Date = as.Date(Date)) %>%
  dplyr::count(Location, timepoint, Date, Type) %>%
  mutate(Date = factor(Date, levels = sort(unique(Date))))

plot = ggplot(check_dates_plot,
       aes(x = Date, y = n, fill = Type)) +
  geom_col(width = 0.8) +
  #facet_grid(rows = vars(Location),cols = vars(timepoint),scales = "free_x",space  = "free_x") +
  facet_grid(Location ~ timepoint, scales = "free_x", space = "free_x") +
  scale_fill_viridis_d(option = "D") +
  #scale_x_date(date_labels = "%m/%d/%Y") +
  labs(
    x = "Sampling date",
    y = "Number of samples",
    fill = "Sample type",
    title = "Samples collected per date by location and timepoint"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA),
    axis.title = element_text(face = "bold")
  )

# Save final heatmap
ggsave(
  plot = plot,
  filename = "Biofilm Project Figures/Kraken2 CT 0.5 Check Dates Shared with Biofilm By Time Point.png",
  width = 9,
  height = 6,
  units = "in"
)


#------------------------------------Filter for dates closest to biofilm collection-----------------------
biofilm_dates <- heatmap_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
  mutate(
    biofilm_ref_date = Date
  ) %>% 
    ungroup() %>% 
  dplyr::distinct(Location,timepoint,biofilm_ref_date)

filtered_dates <- heatmap_df %>%
  mutate(Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 
  

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  )

#plot

check_dates = filtered_dates %>% 
  dplyr::select(SampleID,Location, Type,Date,timepoint,sub_timepoint) %>% distinct()


check_dates_plot <- check_dates %>%
  mutate(Date = as.Date(Date)) %>%
  dplyr::count(Location, timepoint, Date, Type) %>%
  mutate(Date = factor(Date, levels = sort(unique(Date))))

ggplot(check_dates_plot,
       aes(x = Date, y = n, fill = Type)) +
  geom_col(width = 0.8) +
  #facet_grid(rows = vars(Location),cols = vars(timepoint),scales = "free_x",space  = "free_x") +
  facet_grid(Location ~ timepoint, scales = "free_x", space = "free_x") +
  scale_fill_viridis_d(option = "D") +
  #scale_x_date(date_labels = "%m/%d/%Y") +
  labs(
    x = "Sampling date",
    y = "Number of samples",
    fill = "Sample type",
    title = "Samples collected per date by location and timepoint"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA),
    axis.title = element_text(face = "bold")
  )

  
#---------new dataframe filtering approach-----------------------
heatmap_df <- kraken_merged_df_filtered %>%
  filter(taxonomy_lvl == "S") %>% 
  mutate(rel_ab = as.numeric(fraction_total_reads)) %>% 
  filter(Location %in% locations) %>% 
  filter(name %in% c(rpip_targets)) %>% 
  filter(Type %in% c("Drain", "Endcap", "H_WW")) %>% 
  
  # #filter dates
  filter(
    (Type == "H_WW" & Date %in% dates) |
    (Type != "H_WW")
  ) %>%
  
  mutate(
    Type = factor(Type, levels = c("Drain", "Endcap", "H_WW")),
    TypeLabel = case_when(
      Type == "H_WW" ~ "Hospital\nWastewater", 
      TRUE ~ as.character(Type)
    ),
    # Make Date a string in MM/DD/YYYY format
    DateString = format(as.Date(Date), "%m/%d/%Y")
  ) %>% 
  
  left_join(sixteens_df,by = c("Type","Location","Date")) %>% 
  
  drop_na(timepoint)


biofilm_dates <- heatmap_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
  mutate(
    biofilm_ref_date = Date
  ) %>% 
    ungroup() %>% 
  dplyr::distinct(Location,timepoint,biofilm_ref_date)

filtered_dates <- heatmap_df %>%
  mutate(Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  )


analysis_df <- kraken_merged_df_filtered %>%
  filter(
    taxonomy_lvl == "S",
    Location %in% locations,
    name %in% rpip_targets,
    Type %in% c("Drain", "Endcap", "H_WW"),
     (Type == "H_WW" & Date %in% dates) | (Type != "H_WW")
  ) %>%
  mutate(
    rel_ab = as.numeric(fraction_total_reads),
    Date   = as.Date(Date),
    Type   = factor(Type, levels = c("Drain", "Endcap", "H_WW"))
  ) %>%
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>%
  drop_na(timepoint)


biofilm_dates <- analysis_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
  mutate(
    biofilm_ref_date = Date
  ) %>% 
    ungroup() %>% 
  dplyr::distinct(Location,timepoint,biofilm_ref_date)

filtered_dates <- analysis_df %>%
  mutate(Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  ) %>% 
  
  dplyr::select(-best_date,-date_diff_days,-biofilm_ref_date)


```

#NMDS and PERMANOVA
##Microbiome NMDS only on samples with species-level call-RPIP PATHOGENS ONLY-ONLY SHARED DATES WITH BIOFILM
```{r}
#specify dates
dates = c("2024-09-09", "2024-09-11", "2024-09-12", "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09")


# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")


# RPIP taxa
rpip_targets <-  rpip_pathogen_species_list

analysis_df <- kraken_merged_df_filtered %>%
  filter(
    taxonomy_lvl == "S",
    Location %in% locations,
    name %in% rpip_targets,
    Type %in% c("Drain", "Endcap", "H_WW"),
     (Type == "H_WW" & Date %in% dates) | (Type != "H_WW")
  ) %>%
  mutate(
    rel_ab = as.numeric(fraction_total_reads),
    Date   = as.Date(Date),
    Type   = factor(Type, levels = c("Drain", "Endcap", "H_WW"))
  ) %>%
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>%
  drop_na(timepoint)


biofilm_dates <- analysis_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
  mutate(
    biofilm_ref_date = Date
  ) %>% 
    ungroup() %>% 
  dplyr::distinct(Location,timepoint,biofilm_ref_date)

filtered_dates <- analysis_df %>%
  mutate(Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  ) %>% 
  
  dplyr::select(-best_date,-date_diff_days,-biofilm_ref_date)

metadata = filtered_dates %>%
  
  #filter dates
  dplyr::select(SampleID,Location,Type,Month,Date) %>% distinct() %>% 
  
  
  column_to_rownames(var = "SampleID")

# 1) Extract species-level rows (where clade_name contains "s__" OR length of pipe-split >=7)
genus_df <- filtered_dates %>%
  
  filter(name %in% rpip_targets)  



# make dataframe with a column per sample and values are rpkm
transposed = genus_df %>% dplyr::select(SampleID,name,rel_ab) %>% 
  
  #as numeric
  mutate(rel_ab = as.numeric(rel_ab)) %>% 
  
  group_by(SampleID, name) %>% 
  summarize(rel_ab = sum(rel_ab, na.rm = TRUE), .groups = "drop") %>%
  ungroup() %>% 
  
  pivot_wider(
    names_from = SampleID,   # new column names come from sample_id
    values_from = rel_ab        # values in cells come from rpkm
  ) %>% 
  
  #replace NAs with zeros
  mutate_all( ~replace(., lengths(.)==0, 0)) 


#store data in counts variable and only show fraction reads 
counts = transposed %>% dplyr::select(-name) %>% t() 

# Extract row names as a separate vector
row_names <- rownames(counts)


#remove samples
#counts <- counts[!(rownames(counts) %in% c("S209", "S210","S1","S2","S3")), ]

#make sure it's a numeric dataframe
counts = as.data.frame(counts) %>% mutate(across(everything(), ~ as.numeric(as.character(.))))

#replace NAs with zero
counts[is.na(counts)] = 0
                       
# If too few samples remain, think twice
cat("Samples with species-level data:", nrow(counts), "\n")   

# Optional taxon filtering: remove very rare species
min_prevalence <- 0.05  # species present in >=5% samples
present_in <- colSums(counts > 0)
keep_taxa <- present_in >= (min_prevalence * nrow(counts))
comm_sp_filt <- counts[, keep_taxa, drop = FALSE]

# Hellinger transform (recommended for species abundance)
comm_sp_hell <- decostand(comm_sp_filt, method = "hellinger")

#get the distance matrix
dist1 <- vegdist(counts,method = "bray")
#dist1 = vegdist(comm_sp_hell, method = "bray") # if we want Hellinger transformation and removal of rare species
argMDS <- isoMDS(dist1, trace = F)
scores <- as.data.frame(argMDS)
scores <- scores[,1:2]

# Assuming 'scores' and 'metadata' are your dataframes
common_rows <- intersect(rownames(counts), rownames(metadata))

# Filter 'scores' dataframe to keep only rows present in 'metadata'
filtered_counts <- counts[rownames(counts) %in% common_rows, ]


# Filter 'scores' dataframe to keep only rows present in 'metadata'
filtered_metadata <- metadata[rownames(metadata) %in% common_rows, ]



# meta MDS and create envfit variable
vare.mds <- metaMDS(filtered_counts, trace = FALSE, distance = "bray")
#library(vegan)
ef <- envfit(vare.mds, filtered_metadata, permu = 999,na.rm = TRUE)


#plot results
plot(vare.mds, display = "sites")


#extract envfit arrows data and filter out non significant one 
#en_coord_cont = as.data.frame(scores(ef, "vectors")) * ordiArrowMul(ef) 
#en_coord_cont$pval <- ef[["vectors"]][["pvals"]]
#en_coord_cont <- filter(en_coord_cont,pval<=0.10)

#extract nmds scores and associate to metadata 
nmds_scores <- as.data.frame(scores(vare.mds)$sites)

# Assuming your_data is your dataframe
nmds_scores$SampleID <- rownames(nmds_scores)
rownames(nmds_scores) <- NULL # now get rid of row names

full_metadata = metadata %>% rownames_to_column(var = "SampleID") %>% 
  
  #create new variable combined with site and sample type
  mutate("site_type"=paste0(Location,"_",Type))

nmdsMerged = left_join(nmds_scores,full_metadata, by = "SampleID") 


#specify subareas we want to incude in plot
locations = c("FIONA","SHREK","OSCAR","MARIO","LUIGI")


nmds_plotg <- nmdsMerged %>% 
  
  #filter locations if needed
  #filter(Location %in% locations) %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = site_type)) +
  geom_point() +
  stat_ellipse(linewidth = 0.75) +
  
  #change colors
  #scale_color_manual(values = c(
    #"4" = "lightsalmon", "12" = "grey",
   #                             )) +
  #scale_shape_manual(values = c(1, 2
                                #, 3, 4, 5, 6, 16, 17, 18, 19, 20, 21,22, 29
                                #)) +
  theme_bw() +
  labs(
    title = "NMDS Plot",
    x = "NMDS1",
    y = "NMDS2",
    color = "site_type",
    size = "Collection Date"
  ) + 
  #geom_vline(xintercept = 0,linetype = "dashed")+ geom_hline(yintercept = 0,linetype = "dashed")+
  #geom_segment(aes(x = 0, y = 0, xend =NMDS1, yend = NMDS2), data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  #geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", fontface = "bold", label = row.names(en_coord_cont)) +
  theme(axis.text.y = element_text(angle = 45, hjust = 1)) + theme(axis.text.y = element_blank())

# print plot
nmds_plotg


#prettier NMDS plot with shapes and colors
nmds_plotg <- nmdsMerged %>% 
  
  #filter out municipal ww
  filter(Type != "Mu_WW") %>% 
  
  #rename some variables
  mutate(Type = case_when(Type == "Mu_WW"~ "Municipal WW",
                          Type == "Endcap" ~ "Branch Biofilm",
                          Type == "Drain" ~ "Sink Biofilm",
                          Type == "H_WW" ~ "Wastewater",
                          .default = Type)) %>% 
  
  # now replace Site names 
  mutate(Location = case_when(Location == "MARIO" ~ "Site 1",
                              Location == "FIONA" ~ "Site 2",
                              Location == "LUIGI" ~ "Site 3",
                              Location == "SHREK" ~ "Site 4",
                              Location == "OSCAR" ~ "Site 5")) %>% 
  
  #filter locations if needed
  #filter(Location %in% locations) %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = Type, shape = Location, group = Type)) +
  geom_point() +
  stat_ellipse(linewidth = 0.75) +
  
  #change colors
  #scale_color_manual(values = c(
    #"4" = "lightsalmon", "12" = "grey",
   #                             )) +
  #scale_shape_manual(values = c(1, 2
                                #, 3, 4, 5, 6, 16, 17, 18, 19, 20, 21,22, 29
                                #)) +
  theme_bw() +
  labs(
    title = "NMDS Plot",
    x = "NMDS1",
    y = "NMDS2",
    color = "site_type",
    size = "Collection Date"
  ) + 
  
  #manualy change colors
  scale_color_manual(values = c("Sink Biofilm" = "#79307D", "Branch Biofilm" = "#417C8C", 
                                "Municipal WW" = fortyfive_pal[[3]], 
                                "Wastewater" = "#E57262"
                                )) +
  #scale_shape_manual(values = c("Drain" = 1, "Endcap" = 2, "Municipal WW" = 3,  "Hospital WW" = 4 )) +
  #geom_vline(xintercept = 0,linetype = "dashed")+ geom_hline(yintercept = 0,linetype = "dashed")+
  #geom_segment(aes(x = 0, y = 0, xend =NMDS1, yend = NMDS2), data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  #geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", fontface = "bold", label = row.names(en_coord_cont)) +
  labs(title = "Kraken2 NMDS Analysis",
       color = "",
       shape = "") +
  theme(axis.text.y = element_text(angle = 45, hjust = 1)) + 
  theme(axis.text.y = element_blank(),
        axis.title = element_text(size= 20,face = "bold",color = "black"),
        legend.text = element_text(size = 15,face = "bold", color = "black"))

# print plot
nmds_plotg

# Save the combined plot as a PNG file
ggsave(plot = nmds_plotg, path = "Biofilm Project Figures", "Kraken2 CT 0.5- Microbiome Bracken derived NMDS (species level) by Site and Location-ONLY SHARED DATES WITH BIOFILM-RPIP only.png", dpi = 300, height = 5, width = 7.5, units = "in")

ggsave(plot = nmds_plotg, path = "Biofilm Project Figures", "Kraken2 CT 0.5- Microbiome Bracken derived NMDS (species level) by Site and Location-ONLY SHARED DATES WITH BIOFILM-RPIP only.svg", height = 6, width = 7)


```




## Run permanova on NMDS Microbiome-ONLY SHARED DATES WITH BIOFILM
```{r}

#specify subareas we want to incude
locations = c("FIONA","SHREK","OSCAR","MARIO","LUIGI")

# Find the row indices where "Influent" or "Septic" appear in the "subarea" column
rows <- which(filtered_metadata[, 1] %in% locations)

# Subset the matrix based on the selected rows (so filtering for specific sites)
subset_metadata <- filtered_metadata[rows, ] %>% 
  
  #create new variable combined with site and sample type
  mutate("site_type"=paste0(Location,"_",Type))

# Get row names that are present in both data frames
common_rows <- intersect(rownames(subset_metadata), rownames(filtered_counts))

subset_adonis = filtered_counts[common_rows,]

  

adonis2(subset_adonis ~ Location + Type , subset_metadata, permutations = 9999, by = "margin")


#run pairwise adonis using function
pairwise_adonis_microbiome = pairwise.adonis2(subset_adonis, subset_metadata$Location, sim.method = 'bray', p.adjust.m = 'bonferroni')

#adonis.pair(vegdist(subset_adonis),subset_metadata$subarea)


#run pairwise adonis using function
pairwise_adonis_microbiome = pairwise.adonis2(subset_adonis, subset_metadata$Type, sim.method = 'bray', p.adjust.m = 'bonferroni')

#adonis.pair(vegdist(subset_adonis),subset_metadata$subarea)


subset_adonisMGX= subset_adonis

subset_metadataMGX= subset_metadata


# Run PERMANOVA
res <- adonis2(subset_adonis ~ Location * Type, subset_metadata,
               permutations = 9999, by = "terms")

# Convert to data frame
res_df <- as.data.frame(res) %>% 
  # Add a column for term names
  mutate(Term = rownames(res),
         # Replace very small p-values with <0.0001 for readability
         `p-value` = ifelse(`Pr(>F)` < 0.0001, "<0.0001", `Pr(>F)`),
         # Round numeric columns nicely
         R2 = round(R2, 4),
         F = round(F, 3),
         SumOfSqs = round(SumOfSqs, 3)) %>%
  # Select and reorder columns
  dplyr::select(Term, Df, SumOfSqs, R2, F, `p-value`)

# View the data frame
res_df


library(dplyr)
library(ggplot2)
library(scales)

df_plot <- res_df %>%
  filter(Term != "Total") %>%        # remove Total row
  mutate(Percent = R2)     %>%           # R2 already proportions
  mutate(Term = case_when(Term == "Residual" ~ "Unknown",
                          .default = Term))

df_plot$Term <- factor(df_plot$Term,
                       levels = c("Location", "Type", "Location:Type", "Unknown"))


plot = ggplot(df_plot, aes(x = "Explained Variation", y = Percent, fill = Term)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = NULL,
       y = "Explained Variation (%)") +
  theme_classic()

print(plot)

# Save with ggsave
ggsave("Biofilm Project Figures/Kraken CT 0.5 All microbiome all dates Explained Permanova Variation-oNLY DATES SHARED WITH BIOFILM.svg", plot = plot, width = 4, height = 6)

# Save with ggsave
ggsave("Biofilm Project Figures/Kraken CT 0.5 All microbiome all dates Explained Permanova Variation-oNLY DATES SHARED WITH BIOFILM.png", plot = plot, width = 4, height = 6,dpi = 300)




#
library(vegan)

dist_mat <- vegdist(subset_adonis, method = "bray")

dispersion <- betadisper(dist_mat, subset_metadata$Type)
anova(dispersion)
permutest(dispersion)

library(ggplot2)

disp_df <- data.frame(
  Distance = dispersion$distances,
  Type = subset_metadata$Type
)

ggplot(disp_df, aes(x = Type, y = Distance)) +
  geom_boxplot()

```

#Diversity
##Diversity of Microbiome-RPIP Pathogens only- ONLY SHARED DATES WITH BIOFILM
```{r}
#print number of taxa

n_taxa_df =kraken_rel_ab %>%
  group_by(SampleID) %>%
  summarise(n_taxa = sum(fraction_total_reads > 0)) %>%
  arrange(n_taxa)


#specify dates
dates = c("2024-09-09", "2024-09-11", "2024-09-12", "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09")



# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")


# RPIP taxa
rpip_targets <-  rpip_pathogen_species_list

analysis_df <- kraken_merged_df_filtered %>%
  filter(
    taxonomy_lvl == "S",
    Location %in% locations,
    name %in% rpip_pathogen_species_list,
    Type %in% c("Drain", "Endcap", "H_WW"),
     (Type == "H_WW" & Date %in% dates) | (Type != "H_WW")
  ) %>%
  mutate(
    rel_ab = as.numeric(fraction_total_reads),
    Date   = as.Date(Date),
    Type   = factor(Type, levels = c("Drain", "Endcap", "H_WW"))
  ) %>%
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>%
  drop_na(timepoint)


biofilm_dates <- analysis_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
  mutate(
    biofilm_ref_date = Date
  ) %>% 
    ungroup() %>% 
  dplyr::distinct(Location,timepoint,biofilm_ref_date)

filtered_df <- analysis_df %>%
  mutate(Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  ) %>% 
  
  dplyr::select(-best_date,-date_diff_days,-biofilm_ref_date)


# 1. sample IDs that pass all your upstream filters
sample_ids <- filtered_df %>% distinct(SampleID) %>% pull(SampleID)

# 2. build wide from the un-RPIP-filtered species table, restricted to those samples
wide <- kraken_merged_df_filtered %>%
  filter(
    taxonomy_lvl == "S",
    SampleID %in% sample_ids
    # NOTE: no `name %in% rpip_pathogen_species_list` here
  ) %>%
  dplyr::select(SampleID, name, reads_clade) %>%
  mutate(reads_clade = as.integer(reads_clade)) %>%
  group_by(SampleID, name) %>%
  summarise(reads_clade = sum(reads_clade, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = name, values_from = reads_clade, values_fill = 0) %>%
  column_to_rownames("SampleID")

# 3. rarefy DOWN (pick a deliberate cutoff if min is tiny)
depths  <- rowSums(wide)
raremin <- min(depths[depths > 0])     # consider summary(depths) and a hard floor instead
wide_mat <- as.matrix(wide[depths >= raremin, ])
storage.mode(wide_mat) <- "integer"

Srare <- vegan::rrarefy(wide_mat, raremin)   # only this call — no second overwrite

shannon  <- vegan::diversity(Srare, index = "shannon")
simpson  <- vegan::diversity(Srare, index = "simpson")
richness <- vegan::specnumber(Srare)
evenness <- ifelse(richness > 1, shannon / log(richness), NA_real_)

diversity <- tibble(
  SampleID = rownames(Srare),
  Shannon  = shannon,
  Simpson  = simpson,
  Richness = richness,
  Evenness = evenness
) %>%
  left_join(metadata_df, by = "SampleID") %>%
  distinct() %>%
  mutate(site_type = paste0(Location, "_", Type))


# Calculate the median values for each group
medians <- aggregate(Shannon ~ site_type, diversity, median)

# Reorder the levels of the 'Group' factor based on the median values
diversity$site_type <- factor(diversity$site_type, levels = medians$site_type[order(medians$Shannon)])


# site_type_diversity_plot = ggbetweenstats(
#   data = diversity,
#   x    = site_type, 
#   y    = Shannon, 
#   outlier.tagging = T,
#   type = "nonparametric", 
#   p.adjust.method = "bonferroni",
#   conf.level = 0.95,
#   pairwise.display = "none",
#   mean.plotting = TRUE,
#   mean.ci = FALSE,
#   #point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.4, size = 3, stroke = 0),
#   #median.point.args = list(size = 0.2, color = "pink"),
#   mean.label.args = list(size = 0.5)
# ) + 
#   ylab("Shannon Diversity for Microbiome")+
#   xlab("Site and Sample Type")+
#   labs(color = "")+
#   guides(color = "none")+
#   theme_classic()+
#   theme(axis.text.x= element_text(angle = 45,size = 7,vjust = 0.9,hjust = 0.9,color= "black"))
# 
# print(site_type_diversity_plot)


# # Save the combined plot as a PNG file
# ggsave(plot = site_type_diversity_plot, path = "Biofilm Project Figures", "Kraken2 CT 0.5 Based Microbiome Diversity by Sample Type and Location.png", dpi = 300, height = 6, width = 7, units = "in")


#new ggstats plot wiht diversity by sample type

tidied_diversity = diversity %>%
  
  filter(Type != "Mu_WW") %>% 
  
  mutate(Type = case_when(Type =="H_WW"~ "Hospital WW",Type == "Mu_WW" ~ "Municipal_WW",
                          Type =="Endcap" ~ "Sewer Biofilm",
                          Type == "Drain" ~ "Sink Biofilm",
                          .default = Type))

# Calculate the median values for each group
medians <- aggregate(Shannon ~ Type, tidied_diversity, median)

# Reorder the levels of the 'Group' factor based on the median values
tidied_diversity$Type <- factor(tidied_diversity$Type, levels = medians$Type[order(medians$Shannon)])


Location_diversity_plot =ggbetweenstats(
  data = tidied_diversity,
  x    = Type, 
  y    = Shannon, 
  outlier.tagging = T,
  type = "nonparametric", 
  p.adjust.method = "bonferroni",
  conf.level = 0.95,
  #pairwise.display = "none",
  mean.plotting = TRUE,
  mean.ci = FALSE,
  #point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.4, size = 3, stroke = 0),
  #median.point.args = list(size = 0.2, color = "pink"),
  mean.label.args = list(size = 0.5)
) + 
  ylab("Shannon Diversity for Microbiome")+
  xlab("")+
  labs(color = "")+
  guides(color = "none")+
  theme_classic()+
  theme(axis.text.x= element_text(size = 9,face = "bold",color= "black"))


print(Location_diversity_plot)

# Save the combined plot as a PNG file
ggsave(plot = Location_diversity_plot,path = "Biofilm Project Figures", "Kraken2 CT 0.5: Kraken-Based Microbiome Diversity by Sample Type-ONLY SHARED DATES WITH BIOFILM-RPIP only.png", dpi = 300, height = 6, width = 7, units = "in")



table(rowSums(wide))                          # distribution of "richness"
diversity %>% filter(Shannon == 0) %>%        # check which samples
  left_join(tibble(SampleID = rownames(wide),
                   n_species = rowSums(wide)),
            by = "SampleID")
```


#Venn Diagrams
## Venn diagram for sites shared-ONLY SHARED DATES WITH BIOFILM-RPIP only
```{r}

#specify dates
dates = c("2024-09-09", "2024-09-11", "2024-09-12", "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09")

# ---- Prepare genus-level matrix ----
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")

rpip_targets <- rpip_pathogen_species_list

analysis_df <- kraken_merged_df_filtered %>%
  filter(
    taxonomy_lvl == "S",
    Location %in% locations,
    name %in% rpip_targets,
    Type %in% c("Drain", "Endcap", "H_WW"),
     (Type == "H_WW" & Date %in% dates) | (Type != "H_WW")
  ) %>%
  mutate(
    rel_ab = as.numeric(fraction_total_reads),
    Date   = as.Date(Date),
    Type   = factor(Type, levels = c("Drain", "Endcap", "H_WW"))
  ) %>%
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>%
  drop_na(timepoint)


biofilm_dates <- analysis_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
   mutate(
    biofilm_ref_date = Date
  ) %>% 
    ungroup() %>% 
  dplyr::distinct(Location,timepoint,biofilm_ref_date)

filtered_df <- analysis_df %>%
  mutate(Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  ) %>% 
  
  dplyr::select(-best_date,-date_diff_days,-biofilm_ref_date)


library(dplyr)
library(ggvenn)

# Colors
earthy_green  <- "#9FB498"
earthy_purple <- "#8488AC"
earthy_orange <- "#E57262"

viridisCol1 = "#79307D"
viridisCol2 = "#417C8C"
viridsCol3 = "#E57262"

# Subareas to include
locations <- c("FIONA","SHREK","OSCAR","MARIO","LUIGI")



#overlap for total pathogens across all Sites------------------------


# Species-level rows only
venn_df <- filtered_df %>%

  filter(name %in% rpip_targets) 



  df <- venn_df %>%
    ungroup() %>%
    distinct(Type, SampleID, name, rel_ab) %>%
    filter(rel_ab != 0) %>% 
    mutate(presence = 1) %>% dplyr::select(-rel_ab) %>% distinct()
  
  df <- filtered_df %>%
  filter(name %in% rpip_targets) %>%
  distinct(Type, SampleID, name, rel_ab) %>%
  filter(rel_ab != 0) %>%
  dplyr::select(Type, name) %>%
  distinct()
  
  sets_total <- list(
  Sink_Biofilm   = df %>% filter(Type == "Drain") %>% pull(name),
  Sewer_Biofilm = df %>% filter(Type == "Endcap") %>% pull(name),
  Hospital_WW   = df %>% filter(Type == "H_WW") %>% pull(name)
)

all_species <- unique(unlist(sets_total))

presence_mat <- sapply(sets_total, function(s) all_species %in% s)

region_df <- tibble(
  region = apply(presence_mat, 1, function(x)
    paste(which(x), collapse = "_"))
) %>%
  dplyr::count(region, name = "n")

  
print(region_df)


library(eulerr)

fit_total <- euler(c(
  "Sink_Biofilm"   = region_df$n[region_df$region == "1"],
  "Sewer_Biofilm" = region_df$n[region_df$region == "2"],
  "Hospital_WW"   = region_df$n[region_df$region == "3"],

  "Sink_Biofilm&Sewer_Biofilm" =
    region_df$n[region_df$region == "1_2"],

  "Sink_Biofilm&Hospital_WW" =
    region_df$n[region_df$region == "1_3"],

  "Sewer_Biofilm&Hospital_WW" =
    region_df$n[region_df$region == "2_3"],

  "Sink_Biofilm&Sewer_Biofilm&Hospital_WW" =
    region_df$n[region_df$region == "1_2_3"]
))

plot(
  fit_total,
  fills = c(viridisCol1, viridisCol2, viridsCol3),
  edges = TRUE,
  labels = list(font = 2, cex = 0.9),
  quantities = list(type = "counts", font = 2, cex = 0.9)
)

# Capture whatever is currently on the device
g <- grid.grab()

# Save with ggsave
ggsave("Biofilm Project Figures/Kraken CT 0.5 Based Euler diagram Total -Counts -ONLY SHARED DATES WITH BIOFILM-RPIP only.svg", plot = g, width = 6, height = 6)

# Save with ggsave
ggsave("Biofilm Project Figures/Kraken CT 0.5 Based Euler diagram Total -Counts -ONLY SHARED DATES WITH BIOFILM-RPIP only.png", plot = g, width = 6, height = 6,dpi = 300)

#percentages
plot(
  fit_total,
  fills = c(viridisCol1, viridisCol2, viridsCol3),
  edges = TRUE,
  labels = list(font = 2, cex = 0.9),
  quantities = list(type = "percent", font = 2, cex = 0.9)
)

# Capture whatever is currently on the device
g <- grid.grab()

# Save with ggsave
ggsave("Biofilm Project Figures/Kraken CT 0.5 Based Euler diagram Total -Percents -ONLY SHARED DATES WITH BIOFILM-RPIP only.svg", plot = g, width = 6, height = 6)

# Save with ggsave
ggsave("Biofilm Project Figures/Kraken CT 0.5 Based Euler diagram Total -Percents -ONLY SHARED DATES WITH BIOFILM-RPIP only-percents.png", plot = g, width = 6, height = 6,dpi = 300)

  
  # Build Venn input
  ggvenn_df <- list(
    Sink_Biofilm   = df %>% filter(Type == "Drain") %>% pull(name),
    Sewer_Biofilm = df %>% filter(Type == "Endcap") %>% pull(name),
    Hospital_WW  = df %>% filter(Type == "H_WW") %>% pull(name)
  )

  # Plot with a title for clarity
  print(
    ggvenn(
      ggvenn_df,
      fill_color = c(viridisCol1, viridisCol2, viridsCol3),
      stroke_size = 0.5,
      set_name_size = 4
    ) +
    ggtitle("Species overlap")
  )

```



#Core Sewer Biofilm Microbiome
##Reset Themes
```{r}
library(ggplot2)

# Reset global theme and font
theme_set(theme_grey(base_family = "", base_size = 11))

# Remove theme_update() overrides
theme_update(
  axis.text.x        = element_text(),
  axis.text.y        = element_text(),
  axis.ticks.x       = element_line(),
  axis.ticks.y       = element_line(),
  axis.ticks.length.x = unit(0.15, "lines"),
  panel.grid         = element_line(),
  plot.background    = element_rect()
)

# Reset geom defaults that may have inherited "Chivo"
update_geom_defaults("text",  list(family = ""))
update_geom_defaults("label", list(family = ""))

showtext::showtext_auto(FALSE)

```

## Core microbiome of sewer biofilm plots-ONLY DATES SHARED WITH BIOFILM
```{r}
library("devtools")
#install_github("microbiome/microbiome")
library(microbiome)
library(RColorBrewer)
library(reshape)


# data(peerj32)
# 
# # Rename the data
# pseq <- peerj32$phyloseq
# 
# # Calculate compositional version of the data
# # (relative abundances)
# pseq.rel <- microbiome::transform(pseq, "compositional")
# 
# head(prevalence(pseq.rel, detection = 1/100, sort = TRUE))
# 
# p.seq.prev = prevalence(pseq.rel, detection = 1/100, sort = TRUE)

# -----------------------Full Heatmap Code by Month-----------------------

#specify dates
dates = c("2024-09-09", "2024-09-11", "2024-09-12", "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09")


# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")


# RPIP taxa
rpip_targets <-  rpip_pathogen_species_list

# ---------------------------STEP 1: Extract and clean-------------------------
analysis_df <- kraken_merged_df_filtered %>%
  filter(
    taxonomy_lvl == "S",
    Location %in% locations,
    name %in% rpip_targets,
    Type %in% c("Drain", "Endcap", "H_WW"),
     (Type == "H_WW" & Date %in% dates) | (Type != "H_WW")
  ) %>%
  mutate(
    rel_ab = as.numeric(fraction_total_reads),
    #rel_ab = log1p(rel_ab),
    Date   = as.Date(Date),
    Type   = factor(Type, levels = c("Drain", "Endcap", "H_WW"))
  ) %>%
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>%
  drop_na(timepoint)


biofilm_dates <- analysis_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
  mutate(
      biofilm_ref_date = Date
    ) %>% 
    ungroup() %>% 
  dplyr::distinct(Location,timepoint,biofilm_ref_date)


filtered_df <- analysis_df %>%
  mutate(Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  ) %>% 
  
  dplyr::select(-best_date,-date_diff_days,-biofilm_ref_date)


heatmap_df <- filtered_df %>%
  filter(Type %in% c("Endcap")) %>% 
  
  mutate(
    presence = if_else(rel_ab > 0, 1, 0),
    Type = case_when(
      Type == "H_WW"   ~ "Hospital Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain"  ~ "Sink Biofilm",
      .default = Type
    ),
    Month = format(as.Date(Date), "%m-%Y")  # X-axis by month
  ) %>% 
  
  dplyr::select(SampleID,name,Location,Type,Date, rel_ab,presence) %>% 
  
   # now replace Site names 
  mutate(Location = case_when(Location == "MARIO" ~ "Site 1",
                              Location == "FIONA" ~ "Site 2",
                              Location == "LUIGI" ~ "Site 3",
                              Location == "SHREK" ~ "Site 4",
                              Location == "OSCAR" ~ "Site 5",
                              .default = Location))





library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# -----------------------STEP 2: Compute prevalence across relative abundance thresholds-----------------------
thresholds <- seq(0.0005, 0.01, by = 0.0005)  # 0% to 5% relative abundance

#thresholds <- seq(0.2, 5, by = 0.25)  # 0% to 2% relative abundance

prevalence_df <- expand.grid(
  name = unique(heatmap_df$name),
  threshold = thresholds
) %>%
  as_tibble() %>%
  left_join(
    heatmap_df %>% mutate(rel_ab_pct = rel_ab * 100),  # convert to %
    by = "name"
  ) %>% 
  group_by(name, threshold) %>%
  summarise(
  prevalence = mean(rel_ab_pct > threshold, na.rm = TRUE) * 100
  ,.groups = "drop"
  )

prevalence_df_filtered <- prevalence_df %>%
  group_by(name) %>%
  filter(prevalence[threshold == 0.0005] > 0.1) %>%
  ungroup()


# ------------------------STEP 3: Order species by prevalence at highest threshold-------------------------
species_order <- prevalence_df_filtered %>%
  group_by(name) %>% 
  summarise(sum_prevalence = sum(prevalence,na.rm = TRUE)) %>%
  ungroup() %>% 
  arrange(sum_prevalence)  %>%  # lowest prevalence at top, highest at bottom
  pull(name)

# species_order <- prevalence_df %>%
#   filter(threshold == mean(threshold,na.rm = TRUE)) %>%
#   arrange(prevalence) %>%  # lowest prevalence at top, highest at bottom
#   pull(name)

prevalence_df_filtered$name <- factor(prevalence_df_filtered$name, levels = species_order)

# -------------------------STEP 4: Core microbiome heatmap-----------------------
# Gray palette: dark = high prevalence, light = low prevalence
# gray_palette <- gray(seq(0.8, 0.2, length.out = 5))
# 
# fortyfive_pal = c("#a0cb6b","#8368cb","#c86c69","#cdd3e5","#dab594","#d692d1",
#                "#7495c3","#9fdeca","#e2e8c3","#d8a4af","#71bed2","#bca9dd",
#                "#8bb598","#e5cbd4","#6d7ecd","#e4d5ca","#a8dfa5","#a0bada",
#                "#cbca6f","#c6926f","#cce7e6","#81a4b2","#ca69c3","#76bd75",
#                "#d37dae","#abb684","#6ecda4","#a88bbb","#8cbbb6","#d18191",
#                "#c1d8c2","#dac4e6","#e0adce","#a96dca","#e2bdb6","#aacedd",
#                "#9992d9","#e0d0ab","#a2abdf","#b88b9e","#b4a382","#dba294",
#                "#c9dd9b","#c8b67c","#b9847b")
# 
# 
# med_col  <- "#a0cb6b"
# high_col <- "#8368cb"
# 
# med_col  <- "#e2bdb6"   # soft peach
# high_col <- "#6d7ecd"   # periwinkle blue

med_col  <- "#9fdeca"   # light teal
high_col <- "#084081"   # keep your dark anchor OR use "#6d7ecd"
high_col = "#6d7ecd"
low_col = "#f7fcf0"

plot = ggplot(prevalence_df_filtered, aes(x = threshold, y = name, fill = prevalence)) +
geom_tile(color = "white", linewidth = 0.01)+
  #scale_fill_gradientn(colors = gray_palette, name = "Prevalence (%)") +
  
scale_fill_gradientn(
  colors = c("#f7fcf0", med_col, high_col),  # green-blue ramp (very safe)
  #colors = c("#FDE725", "#21908C", "#440154"),
  values = scales::rescale(c(0, 50, 100)),
  limits = c(0, 100),
  breaks = seq(0, 100, 10),
  labels = function(x) paste0(x, "%"),
  name = "Prevalence"
)+
# scale_fill_gradientn(
#   #colors = rev(RColorBrewer::brewer.pal(11, "Spectral")),
#   #colors = viridis(100, option = "magma"),
#   # colors = viridis(100, option = "cividis"),  # especially strong for color blindness
#   #colors = viridis(100, option = "plasma"),
#   #colours = rev(brewer.pal(5, "Spectral")),
#   limits = c(0, 100),
#   breaks = seq(0, 100, 10),
#   labels = function(x) paste0(x, "%"),
#   name = "Prevalence"
# )+
    scale_x_continuous(
    labels = function(x) paste0(x, "%")
  ) +
  # scale_fill_distiller(
  #   palette = "Spectral",
  #   direction = -1,   # reverse if you want high = red
  #   name = "Prevalence (%)"
  # ) +
  #scale_fill_viridis_c(option = "viridis", direction = -1, name = "Prevalence (%)") +
  labs(
    x = "Detection Threshold",
    y = "",
    title = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 6.5),
    axis.ticks = element_line(color = "white"),
    axis.text.y = element_text(size = 7.5, face = "italic"),
    axis.title = element_text(size = 11),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 10)
  )

print(plot)

# Save with ggsave
ggsave("Biofilm Project Figures/Kraken CT 0.5 Based Sewer Biofilm Core Microbiome Heat Map -Counts -ONLY SHARED DATES WITH BIOFILM-RPIP only.svg", plot = plot, width = 5.5, height = 3)

# Save with ggsave
ggsave("Biofilm Project Figures/Kraken CT 0.5 Based Sewer Biofilm Core Microbiome Heat Map -Counts -ONLY SHARED DATES WITH BIOFILM-RPIP only.png", plot = plot, width = 5.5, height = 3,dpi = 300)

```




#16S Data
##get data
```{r}
fig_dir <- 'Biofilm Project Figures/'

ps <- readRDS("Input Data/phyloseq_object_updated.rds")


#ps <- filter_taxa(ps, function(x) sum(x) > 0, TRUE) #keep only seqs that had reads
# Keep only taxa with at least one read
keep_taxa <- taxa_sums(ps) > 0

# Prune taxa that have 0 reads
ps_pruned <- prune_taxa(keep_taxa, ps)

ps_nofactor <- ps_pruned

sample_type_levels <- c('d', 'e', 'tap', 'w')
corresponding_sewer_levels <- c( 'FIONA', 'LUIGI', 'MARIO', 'OSCAR', 'SHREK')
timepoint_levels <- c('1', '2', '3', '4')
sub_timepoint_levels <- c('a', 'b')

sample_data(ps_pruned)$sample_type <- factor(sample_data(ps_pruned)$sample_type, order= FALSE, levels=sample_type_levels)
sample_data(ps_pruned)$corresponding_sewer <- factor(sample_data(ps_pruned)$corresponding_sewer, order= FALSE, levels=corresponding_sewer_levels)
sample_data(ps_pruned)$timepoint <- factor(sample_data(ps_pruned)$timepoint, order= TRUE, levels=timepoint_levels)
sample_data(ps_pruned)$sub_timepoint <- factor(sample_data(ps_pruned)$sub_timepoint, order= FALSE, levels=sub_timepoint_levels)

# normalize to percents
ps_perc <- transform_sample_counts(ps_pruned, function(x) 100*x/sum(x))
ps_perc_nofactor <- transform_sample_counts(ps_nofactor, function(x) 100*x/sum(x))

#summarize counts
seq_summary <- sample_data(ps_pruned) %>% group_by(corresponding_sewer, timepoint) %>% dplyr::summarise(n=n())

```

## Rarefaction
```{r}


otu_tab <- as.data.frame(otu_table(ps_pruned))

# named color vector mapped to ORIGINAL codes
sewer_colors <- c(
  "tap" = "darkorchid",
  "e"   = "coral1",
  "w"   = "goldenrod2",
  "d"   = "lightblue"
)

# relabel sample types
sample_types <- as.character(sample_data(ps_pruned)$sample_type)

sample_types_named <- dplyr::recode(sample_types,
  "tap" = "Tap",
  "e"   = "Sewer Biofilm",
  "w"   = "Wastewater",
  "d"   = "Sink Biofilm"
)

# colors in same order as samples
sample_colors <- sewer_colors[sample_types]

# legend labels + colors
legend_labels <- c("Tap", "Sewer Biofilm", "Wastewater", "Sink Biofilm")
legend_colors <- sewer_colors[c("tap","e","w","d")]

# save (ggsave equivalent for base plot)
jpeg(filename = file.path(fig_dir, "ASV_by_sample_type2.jpeg"),
     width = 10, height = 8, units = "in", res = 150)

rarecurve(otu_tab,
          step = 100,
          col = sample_colors,
          xlab = "Reads",
          ylab = "Amplicon Sequence Variants (ASVs)",
          lwd = 2,
          ylim = c(0, 800),
          label = FALSE)

legend("topright",
       legend = legend_labels,
       col = legend_colors,
       lty = 1,
       lwd = 2)

dev.off()
```

##Core Microbiome Analysis using core microbiome package
```{r}
#---------Libraries----------------
library("devtools")
#install_github("microbiome/microbiome")
library(microbiome)

#------------Test Phyloseq Object-----------
# data(peerj32)
# 
# # Rename the data
# pseq_test <- peerj32$phyloseq
# 
# # Calculate compositional version of the data
# # (relative abundances)
# pseq.rel_test <- microbiome::transform(pseq_test, "compositional")
# 
# otu_table(pseq.rel_test)      # abundance matrix
# tax_table(pseq.rel_test)      # taxonomy
# sample_data(pseq.rel_test)   # metadata

#------Load Phyloseq Object----------



#subset dataset to only sewer biofilm
ps_sewer_biofilm <- subset_samples(
  ps_pruned,
  sample_type %in% c("e")
)

ps_sewer_biofilm <- prune_taxa(
  taxa_sums(ps_sewer_biofilm) > 0,
  ps_sewer_biofilm
)


pseq = ps_sewer_biofilm

ps_genus <- pseq %>%
  microbiome::aggregate_taxa("Genus")

# Calculate compositional version of the data
# (relative abundances)

pseq.rel <- microbiome::transform(ps_genus, "compositional")

#---------Get Population Frequencies-----------

#Relative population frequencies; at 1% compositional abundance threshold:
head(prevalence(pseq.rel, detection = 1/100, sort = TRUE))

#head(prevalence(pseq.rel_test, detection = 1/100, sort = TRUE))

#Absolute population frequencies (sample count):
head(prevalence(pseq.rel, detection = 1/100, sort = TRUE, count = TRUE))

#If you only need the names of the core taxa, do as follows. This returns the taxa that exceed the given prevalence and detection thresholds
core.taxa.standard <- core_members(pseq.rel, detection = 0, prevalence = 50/100)

#A full phyloseq object of the core microbiota is obtained as follows:
pseq.core <- core(pseq.rel, detection = 0, prevalence = .5)

#We can also collapse the rare taxa into an “Other” category
pseq.core2 <- aggregate_rare(pseq.rel, "Genus", detection = 0, prevalence = .5)

#Retrieving the core taxa names from the phyloseq object:
core.taxa <- taxa(pseq.core)

#Total core abundance in each sample (sum of abundances of the core members):
#core.abundance <- sample_sums(core(pseq.rel, detection = .01, prevalence = .95))

#----------Core Line Plots-----------
# With compositional (relative) abundances
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
 #ggplot(d) + geom_point(aes(x, y)) + scale_x_continuous(trans="log10", limits=c(NA,1))


plot_core(pseq.rel, 
          prevalences = prevalences, 
          detections = det, 
          plot.type = "lineplot") + 
  xlab("Relative Abundance (%)")

#----------Core Heat Maps-----------------
# Core with compositionals:
library(RColorBrewer)
library(reshape)

prevalences <- seq(.05, 1, .05)

detections <- round(10^seq(log10(0.01), log10(.2), length = 9), 3);detections

# Also define gray color palette
gray <- gray(seq(0,1,length=5))


#Added pseq.rel, I thin... must be checked if it was in the the rednred version,; where it is initialized
#pseq.rel<- microbiome::transform(pseq, 'compositional')
#min-prevalence gets the 100th highest prevalence
p <- plot_core(pseq.rel,
               plot.type = "heatmap", 
               colours = gray,
               prevalences = prevalences, 
               detections = detections, 
               min.prevalence = prevalence(pseq.rel, sort = TRUE)[100]) +
  labs(x = "Detection Threshold\n(Relative Abundance (%))") +
    
  #Adjusts axis text size and legend bar height
  theme(axis.text.y= element_text(size=8, face="italic"),
        axis.text.x.bottom=element_text(size=8),
        axis.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))

print(p)

ggsave(path = "Biofilm Project Figures",file= "Amplicon based Core Microbiome Heat map of Sewer Biofilm-Top Taxa.png", plot=p, width=5, height=4.5,dpi = 400)
  ggsave(path = "Biofilm Project Figures",file="Amplicon based Core Microbiome Heat map of Sewer Biofilm-Top Taxa.png", plot=p, width=5, height=4.5)

# -------------Core with absolute counts and horizontal view--------------
  
  # Extract genus column from tax_table as a vector
genera <- as.character(tax_table(pseq.rel)[, "Genus"])

# Create logical vector: TRUE for known genera
known_genera <- !is.na(genera) & genera != "Unknown" & genera != "unclassified" & genera != ""

# Prune taxa
pseq.rel.filtered <- prune_taxa(known_genera, pseq.rel)

pseq.rel.filtered <- microbiome::transform(pseq.rel.filtered, "compositional")


  
detections <- seq(from = 50, to = round(max(abundances(pseq)), -1), by = 100); detections


      
detections <- seq(from = 0, to = round(max(abundances(pseq)), -1), by = 10); detections

detections <- seq(from = 50, to = round(max(abundances(pseq))/10, -1), by = 100); detections

#Define prevalence thresholds
prevalences <- seq(0.05, 1, 0.05)  # 5% to 100%

p <- plot_core(pseq.rel.filtered, 
               plot.type = "heatmap",
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = 0.2,
               horizontal = FALSE) +  # Genus on y-axis
  #scale_x_discrete(labels = detection_labels) +
  theme(axis.text.x= element_text(size=2, face="italic", hjust=1),
        axis.text.y= element_text(size=3),
        axis.title = element_text(size=6),
        legend.text = element_text(size=5),
        legend.title = element_text(size=7))

print(p)

ggsave(path = "Biofilm Project Figures",file="Amplicon based Core Microbiome Heat map of Sewer Biofilm-old plot--All Taxa.png", plot=p, width=5, height=7,dpi = 400)
ggsave(path = "Biofilm Project Figures",file="Amplicon based Core Microbiome Heat map of Sewer Biofilm-old plot--All Taxa.png", plot=p, width=5, height=6)

x_levels <- levels(p$data$Detection);x_levels

x_labels <- x_levels[round(seq(1, length(x_levels), length.out = 6))]

p <- p +
  scale_x_discrete(
    breaks = x_labels,
    labels = function(x) signif(as.numeric(x), 2)  # rounds/shortens numbers
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # rotate for readability


print(p)

max_det <- max(as.numeric(levels(p$data$Detection)))

p <- plot_core(pseq.rel.filtered, 
               plot.type = "heatmap",
               prevalences = prevalences,
               #detections = detections,
               detections = 20,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = 0.5,
               horizontal = FALSE) 

x_levels <- levels(p$data$Detection);x_levels

x_labels <- x_levels[round(seq(1, length(x_levels), length.out = 6))]

max_det <- max(as.numeric(levels(p$data$Detection)))

#filter dataframe for greater than 70% and rel abundance greater than 0.5%-
thresholds <-c(seq(0.0025, 0.05, by = 0.0025)) # 0% to 5% relative abundance

med_col  <- "#a0cb6b"
high_col <- "#8368cb"

med_col  <- "#e2bdb6"   # soft peach
high_col <- "#6d7ecd"   # periwinkle blue

med_col  <- "#9fdeca"   # light teal
high_col <- "#084081"   # keep your dark anchor OR use "#6d7ecd"
high_col = "#6d7ecd"

p <- plot_core(
  pseq.rel.filtered, 
  plot.type = "heatmap",
  prevalences = prevalences,
  detections = thresholds,
  colours = colorRampPalette(
     #colors = c("#FDE725", "#21908C", "#440154")
    c("#f7fcf0", med_col, high_col)
  )(100),
  min.prevalence = 0.75,
  horizontal = FALSE
) +  
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 11),
    legend.text = element_text(size = 8)
  )

print(p)

ggsave(path = "Biofilm Project Figures",file="Amplicon based Core Microbiome Heat map of Sewer Biofilm-All Taxa.png", plot=p, width=5, height=3,dpi = 400)
ggsave(path = "Biofilm Project Figures",file="Amplicon based Core Microbiome Heat map of Sewer Biofilm-All Taxa.svg", plot=p, width=5, height=3)

#----------Get list of genera at prevalence greater than 70% and rel abundance greater than 0.5%-----------
prev <- prevalence(pseq.rel.filtered, detection = 0.005)  # 0.5%


prev_df <- data.frame(
  Genus = names(prev),
  Prevalence = prev
)

pathogens_high_prev <- prev_df %>%
  dplyr::filter(Prevalence > 0.7)

pathogens = pathogens_high_prev %>% dplyr::select(Genus) %>% distinct(Genus) %>% pull(Genus)


```



##Only Sewer Biofilm Stacked bars using 16S data COLLAPSED ACROSS TIME top 44 genera by mean relative abundance colored) -- x axis is unique sample type and date combination as row numbers -WITH UPDATED DATE CLASSIFICATION BY TIMEPOINT
```{r}

#specify dates
dates = c("2024-09-09", "2024-09-11", "2024-09-12", "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09")


# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame


library(dplyr)
library(ggplot2)
library(forcats)
library(viridis)
library(ggh4x)

# ---------- STEP 1: Extract and clean ----------

library(phyloseq)
library(microbiome)
library(dplyr)
library(tidyr)


# 
# sample_data(ps)$sample_type <- as.character(sample_data(ps)$sample_type)
# ps_by_type <- merge_samples(ps, "sample_type")
# 
# sample_data(ps_by_type)$sample_type <- sample_names(ps_by_type)
# ps_by_type_genus <- tax_glom(ps_by_type, taxrank = "Genus")
# ps_by_type_genus_perc <- transform_sample_counts(ps_by_type_genus, function(x) 100 * x / sum(x)) #converts counts to percentages per sample (relative abundance).
# 
# ps_sample_filt1.2 <- filter_taxa(ps_by_type_genus_perc, filterfun(kOverA(1, 2)), TRUE) #keep taxa with ≥2% abundance in at least 1 sample.

#df_long <- psmelt(ps_sample_filt1.2)  # creates data.frame with OTU abundances + metadata

library(reshape2)


library(phyloseq)

library(phyloseq)
library(microbiome)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

# --- 0) Make sure these exist in sample_data(ps) ---
# corresponding_sewer (your Site codes like MARIO/FIONA/etc)
# sample_type (w/e/d/tap)

sample_data(ps)$corresponding_sewer <- as.character(sample_data(ps)$corresponding_sewer)
sample_data(ps)$sample_type         <- as.character(sample_data(ps)$sample_type)

# Optional: map Site codes to labels *before* merging
sample_data(ps)$Site <- dplyr::recode(sample_data(ps)$corresponding_sewer,
                                     "MARIO" = "Site 1",
                                     "FIONA" = "Site 2",
                                     "LUIGI" = "Site 3",
                                     "SHREK" = "Site 4",
                                     "OSCAR" = "Site 5",
                                     .default = sample_data(ps)$corresponding_sewer)

# --- 1) Agglomerate to Genus (so you're plotting genus stacks) ---
ps_genus <- tax_glom(ps, taxrank = "Genus")

# --- 2) Create grouping var that ignores time (Site × sample_type) ---
# Use "__" so we can split safely later
sample_data(ps_genus)$Site_Type <- with(sample_data(ps_genus),
                                       paste(Site, sample_type, sep = "__"))

# --- 3) Collapse over time: sums counts within each Site_Type ---
ps_collapsed <- merge_samples(ps_genus, "Site_Type")

# --- 4) Convert to relative abundance so each bar sums to 1 ---
ps_collapsed_rel <- transform_sample_counts(ps_collapsed, function(x) x / sum(x))

# --- 5) Melt + split Site/Type back out ---
library(tidyr)
library(dplyr)

df_update <- psmelt(ps_collapsed_rel) %>% 
  mutate(Site_Type = as.character(Sample)) %>%   # <- use merged sample name
  separate(Site_Type, into = c("Site", "Type_code"), sep = "__", remove = FALSE) %>%
  mutate(
    Type = case_when(
      Type_code == "w"   ~ "Wastewater",
      Type_code == "e"   ~ "Sewer Biofilm",
      Type_code == "d"   ~ "Sink Biofilm",
      Type_code == "tap" ~ "Tap",
      TRUE ~ Type_code
    )
  ) %>%
  dplyr::select(Site, Type, Type_code, Genus, Abundance)  # keep only what matters
  




# Columns: Sample, OTU, Abundance, taxonomy ranks, sample_data columns


# # ---------- New Way: Assign letters time points ----------
# drain_samples <- df_long %>%
#   filter(Type == "d") %>%
#   distinct(Location, DateString, Type,timepoint,sub_timepoint) %>% 
#   group_by(Location, timepoint,sub_timepoint) %>%
#   mutate(SampleLetter = seq_len(n()),
#          NewDateString = paste0(timepoint,"-",sub_timepoint, "-", SampleLetter)) %>%
#   ungroup() %>%
#   dplyr::select(Location, DateString, Type, NewDateString)

# # ---------- Assign letters for sewer biofilm replicates ----------
# sewer_samples  <- df_long %>%
#   filter(Type == "e") %>%
#   distinct(Location, DateString, Type,timepoint,sub_timepoint) %>% 
#   group_by(Location, timepoint,sub_timepoint) %>%
#   mutate(SampleLetter = seq_len(n()),
#          NewDateString = paste0(timepoint, "-",sub_timepoint,"-", SampleLetter)
#          #NewDateString = timepoint
# ) %>%
#   ungroup() %>%
#   dplyr::select(Location, DateString, Type, NewDateString)

# # ---------- Assign letters for ww replicates ----------
# ww_samples  <- df_long %>%
#   filter(Type == "w") %>%
#   distinct(Location, DateString, Type,timepoint,sub_timepoint) %>% 
#   group_by(Location, timepoint,sub_timepoint) %>%
#   mutate(SampleLetter = seq_len(n()),
#          NewDateString = paste0(timepoint,"-",sub_timepoint, "-", SampleLetter)
#           #NewDateString = timepoint
#         ) %>%
#   ungroup() %>%
#   dplyr::select(Location, DateString, Type, NewDateString)
# 
# # ---------- Assign letters for tap replicates ----------
# tap_samples  <- df_long %>%
#   mutate(sub_timepoint = as.character(sub_timepoint)) %>%   # ← CRITICAL
#   mutate(sub_timepoint = ifelse(is.na(sub_timepoint), "x", sub_timepoint)) %>%
#   filter(Type == "tap") %>%
#   distinct(Location, DateString, Type,timepoint,sub_timepoint) %>% 
#   group_by(Location, timepoint,sub_timepoint) %>%
#   mutate(SampleLetter = seq_len(n()),
#          NewDateString = paste0(timepoint,"-",sub_timepoint, "-", SampleLetter)
#           #NewDateString = timepoint
#          ) %>%
#   ungroup() %>%
#   dplyr::select(Location, DateString, Type, NewDateString)
# 
# # Combine drain and sewer samples
# replicate_samples <- bind_rows(
#   drain_samples,
#   ww_samples,
#   tap_samples,
#   sewer_samples %>% dplyr::select(Location, DateString, Type, NewDateString)
# )
# 
# # ---------- Join back to main dataframe ----------
# df_update <- df_long %>%
#   left_join(replicate_samples, by = c("Location", "DateString", "Type")) %>%
#   mutate(
#     NewDateString = ifelse(is.na(NewDateString) == TRUE, DateString, NewDateString)
#   )


# ---------- Collapse into top species ----------

n_genus = 44

# get top globally by total rel_ab
top_genera <- df_update %>%
  group_by(Genus) %>%
  summarise(total_ab = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_ab)) %>% 
  slice_head(n = n_genus) %>% 
  pull(Genus)

# assign species or "Other"
plot_df <- df_update %>%
  dplyr::select(Genus, Site, Type, Abundance) %>%
  mutate(
    class = factor(ifelse(Genus %in% top_genera, Genus, "Other")),
    Type  = factor(Type, levels = c("Wastewater", "Sewer Biofilm","Sink Biofilm","Tap"))
  ) %>% 
  
  mutate(Type = dplyr::recode(Type,
                       "Sink Biofilm" = "Sink\nBiofilm",
                       "Sewer Biofilm" = "Sewer\nBiofilm"))

# ---------- Colors for Plot ----------

# Assign colors: first 15 get a palette, "Other" gets grey
ref = "Other"
myColors <- fortyfive_pal[1:length(levels(plot_df$class))]
names(myColors) <- levels(plot_df$class)
myColors[names(myColors)==ref] <- "grey"


##actual plotting-----------
library(dplyr)
library(ggplot2)
library(ggh4x)
library(ggtext)   # for element_markdown


# italicize legend labels
genera_labels <- setNames(
  ifelse(levels(plot_df$class) == "Other",
         "Other",
         paste0("*", levels(plot_df$class), "*")),   # markdown italics
  levels(plot_df$class)
)


# Plot using Event (factor) on x → no gaps within each facet


library(scales)  # for percent_format()


plot = ggplot(plot_df, aes(x = Site, y = Abundance, fill = class)) +
  geom_col(position = "stack") +
  facet_nested(. ~ Type, scales = "free_x", space = "free_x", switch = "y") +
  scale_fill_manual(values = myColors,
                    labels = genera_labels) +
  # scale_y_continuous(
  #   breaks = c(0, 0.25, 0.50, 0.75),   # tick marks at 25, 50, 75
  #   labels = c("0","25", "50", "75"),   # plain numbers, no %
  #   limits = c(0, 1)                 # optional: keep to 100%
  # ) +
  labs(x = NULL, y = "Relative Abundance", fill = "", title = "") +
theme_minimal(base_size = 11) +
  guides(fill  = guide_legend(ncol = 2))+
theme(
  panel.background = element_rect(fill = "white", colour = "black"),
  strip.background = element_rect(fill = "white", colour = "black"),
  strip.text = element_text(
    face = "bold", 
    hjust = 0.5,          # center horizontally
    vjust = 0.5,          # center vertically
    lineheight = 0.9      # tighter line spacing if wrapped
  ),
  #axis.text.x        = element_blank(),
  axis.text.x = element_text(size = 5),
  axis.text.y = element_text(size = 6),
  panel.spacing      = unit(0.03, "lines"),
  legend.position    = "right",
  legend.key.size    = unit(0.3, "cm"),
  legend.text        = ggtext::element_markdown(size = 5.5),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  #panel.grid.major.y = element_line(color = "grey80"),
  panel.grid.minor.y = element_blank(),
  panel.ontop        = FALSE   # <- ensures geoms (bars) are drawn over grid lines
)


print(plot)

  ggsave(path = "Biofilm Project Figures",file="16Sstacked bar plot of genera detected across all sample types (collapsed over time).png", plot=plot, width=12, height=6,dpi = 400)
  ggsave(path = "Biofilm Project Figures",file="16S stacked bar plot of genera detected across all sample types (collapsed over time).svg", plot=plot, width=12, height=6)
  
#---------------Isolate to only Genera in Sewer Biofilm------------------

#Collapse into top species

n_genus = 44
  
df_sewer_biofilm = df_update %>% filter(Type_code == "e") 

# get top globally by total rel_ab
top_genera <- df_sewer_biofilm %>%
  group_by(Genus) %>%
  summarise(total_ab = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_ab)) %>%
  slice_head(n = n_genus) %>%
  pull(Genus)

# assign species or "Other"
plot_df <- df_sewer_biofilm %>%
  mutate(
    class = factor(ifelse(Genus %in% top_genera, Genus, "Other")),
    Type  = factor(Type, levels = c("Wastewater", "Sewer Biofilm","Sink Biofilm","Tap"))
  ) %>% 
  
  mutate(Type = dplyr::recode(Type,
                       "Sink Biofilm" = "Sink\nBiofilm",
                       "Sewer Biofilm" = "Sewer\nBiofilm"))

#Colors for Plot

# Assign colors: first 15 get a palette, "Other" gets grey
ref = "Other"
myColors <- fortyfive_pal[1:length(levels(plot_df$class))]
names(myColors) <- levels(plot_df$class)
myColors[names(myColors)==ref] <- "grey"


##actual plotting
library(dplyr)
library(ggplot2)
library(ggh4x)
library(ggtext)   # for element_markdown


# italicize legend labels
genera_labels <- setNames(
  ifelse(levels(plot_df$class) == "Other",
         "Other",
         paste0("*", levels(plot_df$class), "*")),   # markdown italics
  levels(plot_df$class)
)


# Plot using Event (factor) on x → no gaps within each facet

library(scales)  # for percent_format()


plot = ggplot(plot_df, aes(x = Site, y = Abundance, fill = class)) +
  geom_col(position = "stack") +
  facet_nested(.~ Type, scales = "free_x", space = "free_x", switch = "y") +
  scale_fill_manual(values = myColors,
                    labels = genera_labels) +
  # scale_y_continuous(
  #   breaks = c(0, 0.25, 0.50, 0.75),   # tick marks at 25, 50, 75
  #   labels = c("0","25", "50", "75"),   # plain numbers, no %
  #   limits = c(0, 1)                 # optional: keep to 100%
  # ) +
  labs(x = NULL, y = "Relative Abundance", fill = "", title = "") +
theme_minimal(base_size = 12) +
  guides(fill = guide_legend(ncol = 2))+
theme(
  panel.background = element_rect(fill = "white", colour = "black"),
  strip.background = element_rect(fill = "white", colour = "black"),
  strip.text = element_text(
    face = "bold", 
    hjust = 0.5,          # center horizontally
    vjust = 0.5,          # center vertically
    lineheight = 0.9      # tighter line spacing if wrapped
  ),
  #axis.text.x        = element_blank(),
  axis.text.x = element_text(size = 5),
  axis.title.x = element_text(size = 5),
  axis.text.y = element_text(size = 10),
  panel.spacing      = unit(0.03, "lines"),
  legend.position    = "right",
  legend.key.size    = unit(0.3, "cm"),
  legend.text        = ggtext::element_markdown(size = 9),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  #panel.grid.major.y = element_line(color = "grey80"),
  panel.grid.minor.y = element_blank(),
  panel.ontop        = FALSE   # <- ensures geoms (bars) are drawn over grid lines
)


print(plot)


ggsave(path = "Biofilm Project Figures",file="16Sstacked bar plot of genera detected-only sewer biofilm-collapsed over time.png", plot=plot, width=6, height=3.5,dpi = 400)
  ggsave(path = "Biofilm Project Figures",file="16S stacked bar plot of genera detected-only sewer biofilm-collapsed over time.svg", plot=plot, width=6, height=3.5)
  

  
  #-------------Create summary table of percent contribution by pathogen and sample type------------------
genera_summary <- plot_df %>%
  group_by(class, Type) %>%
  summarise(
    MinPercent = min(Abundance, na.rm = TRUE) * 100,
    MaxPercent = max(Abundance, na.rm = TRUE) * 100,
    MeanPercent = mean(Abundance, na.rm = TRUE) * 100,
    MedianPercent = median(Abundance, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  dplyr::rename(Genus = class, SampleType = Type) %>%
  arrange(SampleType, desc(MeanPercent))


```

##Only Sewer Biofilm Stacked bars using 16S data (facet nested by Location and Type) top 44 Families by mean relative abundance colored) -- x axis is unique sample type and date combination as row numbers -WITH UPDATED DATE CLASSIFICATION BY TIMEPOINT
```{r}

#specify dates
dates = c("2024-09-09", "2024-09-11", "2024-09-12", "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09")


# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame


library(dplyr)
library(ggplot2)
library(forcats)
library(viridis)
library(ggh4x)

# ---------- STEP 1: Extract and clean ----------

library(phyloseq)
library(microbiome)
library(dplyr)
library(tidyr)

sample_data(ps)$sample_type <- as.character(sample_data(ps)$sample_type)
ps_by_type <- merge_samples(ps, "sample_type")

sample_data(ps_by_type)$sample_type <- sample_names(ps_by_type)
ps_by_type_genus <- tax_glom(ps_by_type, taxrank = "Genus")
ps_by_type_genus_perc <- transform_sample_counts(ps_by_type_genus, function(x) 100 * x / sum(x)) #converts counts to percentages per sample (relative abundance).

ps_sample_filt1.2 <- filter_taxa(ps_by_type_genus_perc, filterfun(kOverA(1, 2)), TRUE) #keep taxa with ≥2% abundance in at least 1 sample.

#df_long <- psmelt(ps_sample_filt1.2)  # creates data.frame with OTU abundances + metadata

library(reshape2)

df_long <- psmelt(ps) %>%  # creates data.frame with OTU abundances + metadata

  dplyr::select(OTU,Sample,Abundance,sample_id,sample_type,sample_date,timepoint,corresponding_sewer,sub_timepoint,qubit_DNA_ng.uL,Kingdom,Phylum,Class, Order, Family, Genus, Species) %>% distinct %>% 
  dplyr::rename("SampleID" = "sample_id",
                "Type" = "sample_type",
                "Location" = "corresponding_sewer",
                "Date" = "sample_date") %>% 
  mutate(Date = as.Date(as.character(Date), format = "%m/%d/%y"),
        # Make Date a string in MM/DD/YYYY format
        Abundance = as.numeric(Abundance),
        DateString = format(as.Date(Date), "%m/%d/%Y"),
        Month = format(as.Date(Date), "%m-%Y")) %>% 
  
    group_by(SampleID) %>%              # group by sample
  mutate(total_ab = sum(Abundance)) %>%  # convert counts to relative abundance (0-1)
  ungroup() %>% 
  
  group_by(SampleID,Type, OTU) %>%              # group by sample
  mutate(rel_ab = Abundance / total_ab) %>%  # convert counts to relative abundance (0-1)
  ungroup() %>% 
  
  dplyr::select(-Abundance)

# Columns: Sample, OTU, Abundance, taxonomy ranks, sample_data columns


# ---------- New Way: Assign letters time points ----------
drain_samples <- df_long %>%
  filter(Type == "d") %>%
  distinct(Location, DateString, Type,SampleID,timepoint,sub_timepoint) %>% 
  group_by(Location, timepoint,sub_timepoint) %>%
  mutate(SampleLetter = seq_len(n()),
         NewDateString = paste0(timepoint,"-",sub_timepoint, "-", SampleLetter)) %>%
  ungroup() %>%
  dplyr::select(SampleID,Location, DateString, Type, NewDateString)

# ---------- Assign letters for sewer biofilm replicates ----------
sewer_samples  <- df_long %>%
  filter(Type == "e") %>%
  distinct(Location, DateString, Type,SampleID,timepoint,sub_timepoint) %>% 
  group_by(Location, timepoint,sub_timepoint) %>%
  mutate(SampleLetter = seq_len(n()),
         NewDateString = paste0(timepoint, "-",sub_timepoint,"-", SampleLetter)
         #NewDateString = timepoint
) %>%
  ungroup() %>%
  dplyr::select(SampleID,Location, DateString, Type, NewDateString)

# ---------- Assign letters for ww replicates ----------
ww_samples  <- df_long %>%
  filter(Type == "w") %>%
  distinct(Location, DateString, Type,SampleID,timepoint,sub_timepoint) %>% 
  group_by(Location, timepoint,sub_timepoint) %>%
  mutate(SampleLetter = seq_len(n()),
         NewDateString = paste0(timepoint,"-",sub_timepoint, "-", SampleLetter)
          #NewDateString = timepoint
        ) %>%
  ungroup() %>%
  dplyr::select(SampleID,Location, DateString, Type, NewDateString)

# ---------- Assign letters for tap replicates ----------
tap_samples  <- df_long %>%
  mutate(sub_timepoint = as.character(sub_timepoint)) %>%   # ← CRITICAL
  mutate(sub_timepoint = ifelse(is.na(sub_timepoint), "x", sub_timepoint)) %>%
  filter(Type == "tap") %>%
  distinct(Location, DateString, Type,SampleID,timepoint,sub_timepoint) %>% 
  group_by(Location, timepoint,sub_timepoint) %>%
  mutate(SampleLetter = seq_len(n()),
         NewDateString = paste0(timepoint,"-",sub_timepoint, "-", SampleLetter)
          #NewDateString = timepoint
         ) %>%
  ungroup() %>%
  dplyr::select(SampleID,Location, DateString, Type, NewDateString)

# Combine drain and sewer samples
replicate_samples <- bind_rows(
  drain_samples,
  ww_samples,
  tap_samples,
  sewer_samples %>% dplyr::select(SampleID, Location, DateString, Type, NewDateString)
)

# ---------- Join back to main dataframe ----------
df_update <- df_long %>%
  left_join(replicate_samples, by = c("SampleID", "Location", "DateString", "Type")) %>%
  mutate(
    NewDateString = ifelse(is.na(NewDateString) == TRUE, DateString, NewDateString)
  )


# ---------- Collapse into top species ----------

n_family = 44

# get top 25 globally by total rel_ab
top_genera <- df_update %>%
  group_by(Family) %>%
  summarise(total_ab = mean(rel_ab, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_ab)) %>%
  slice_head(n = n_family) %>%
  pull(Family)

# assign species or "Other"
plot_df <- df_update %>%
  dplyr::select(Family, SampleID, Date, Location, Type, rel_ab,timepoint,NewDateString) %>%
  mutate(
    Type = case_when(
      Type == "w"   ~ "Wastewater",
      Type == "e" ~ "Sewer Biofilm",
      Type == "d"  ~ "Sink Biofilm",
      Type == "tap" ~ "Tap", 
      .default = Type
    ),
    class = factor(ifelse(Family %in% top_genera, Family, "Other")),
    Type  = factor(Type, levels = c("Wastewater", "Sewer Biofilm","Sink Biofilm","Tap"))
  ) %>% 
  
  # now replace Site names 
  mutate(Location = case_when(Location == "MARIO" ~ "Site 1",
                              Location == "FIONA" ~ "Site 2",
                              Location == "LUIGI" ~ "Site 3",
                              Location == "SHREK" ~ "Site 4",
                              Location == "OSCAR" ~ "Site 5",
                              .default = Location)) %>% 
  
  mutate(Type = dplyr::recode(Type,
                       "Sink Biofilm" = "Sink\nBiofilm",
                       "Sewer Biofilm" = "Sewer\nBiofilm"))

# ---------- Colors for Plot ----------

# Assign colors: first 15 get a palette, "Other" gets grey
ref = "Other"
myColors <- fortyfive_pal[1:length(levels(plot_df$class))]
names(myColors) <- levels(plot_df$class)
myColors[names(myColors)==ref] <- "grey"


##actual plotting-----------
library(dplyr)
library(ggplot2)
library(ggh4x)
library(ggtext)   # for element_markdown


# italicize legend labels
genera_labels <- setNames(
  ifelse(levels(plot_df$class) == "Other",
         "Other",
         paste0("*", levels(plot_df$class), "*")),   # markdown italics
  levels(plot_df$class)
)


# Plot using Event (factor) on x → no gaps within each facet


library(scales)  # for percent_format()




plot = ggplot(plot_df, aes(x = NewDateString, y = rel_ab, fill = class)) +
  geom_col(position = "stack") +
  facet_nested(Location ~ Type, scales = "free_x", space = "free_x", switch = "y") +
  scale_fill_manual(values = myColors,
                    labels = genera_labels) +
  # scale_y_continuous(
  #   breaks = c(0, 0.25, 0.50, 0.75),   # tick marks at 25, 50, 75
  #   labels = c("0","25", "50", "75"),   # plain numbers, no %
  #   limits = c(0, 1)                 # optional: keep to 100%
  # ) +
  labs(x = NULL, y = "Relative Abundance", fill = "", title = "") +
theme_minimal(base_size = 11) +
  guides(fill  = guide_legend(ncol = 1))+
theme(
  panel.background = element_rect(fill = "white", colour = "black"),
  strip.background = element_rect(fill = "white", colour = "black"),
  strip.text = element_text(
    face = "bold", 
    hjust = 0.5,          # center horizontally
    vjust = 0.5,          # center vertically
    lineheight = 0.9      # tighter line spacing if wrapped
  ),
  #axis.text.x        = element_blank(),
  axis.text.x = element_text(size = 5),
  axis.text.y = element_text(size = 6),
  panel.spacing      = unit(0.03, "lines"),
  legend.position    = "right",
  legend.key.size    = unit(0.3, "cm"),
  legend.text        = ggtext::element_markdown(size = 5.5),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  #panel.grid.major.y = element_line(color = "grey80"),
  panel.grid.minor.y = element_blank(),
  panel.ontop        = FALSE   # <- ensures geoms (bars) are drawn over grid lines
)


print(plot)

  ggsave(path = "Biofilm Project Figures",file="16Sstacked bar plot of families detected.png", plot=plot, width=12, height=6,dpi = 400)
  ggsave(path = "Biofilm Project Figures",file="16S stacked bar plot of families detected.svg", plot=plot, width=12, height=6)
  
#---------------Isolate to only Genera in Sewer Biofilm------------------

#Collapse into top species

n_family = 44
  
df_sewer_biofilm = df_update %>% filter(Type == "e")

# get top globally by total rel_ab
top_genera <- df_sewer_biofilm %>%
  group_by(Family) %>%
  summarise(total_ab = mean(rel_ab, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_ab)) %>%
  slice_head(n = n_family) %>%
  pull(Family)

# assign species or "Other"
plot_df <- df_sewer_biofilm %>%
  dplyr::select(Family, SampleID, Date, Location, Type, rel_ab,timepoint,NewDateString) %>%
  mutate(
    Type = case_when(
      Type == "w"   ~ "Wastewater",
      Type == "e" ~ "Sewer Biofilm",
      Type == "d"  ~ "Sink Biofilm",
      Type == "tap" ~ "Tap", 
      .default = Type
    ),
    class = factor(ifelse(Family %in% top_genera, Family, "Other")),
    Type  = factor(Type, levels = c("Wastewater", "Sewer Biofilm","Sink Biofilm","Tap"))
  ) %>% 
  
  # now replace Site names 
  mutate(Location = case_when(Location == "MARIO" ~ "Site 1",
                              Location == "FIONA" ~ "Site 2",
                              Location == "LUIGI" ~ "Site 3",
                              Location == "SHREK" ~ "Site 4",
                              Location == "OSCAR" ~ "Site 5",
                              .default = Location)) %>% 
  
  mutate(Type = dplyr::recode(Type,
                       "Sink Biofilm" = "Sink\nBiofilm",
                       "Sewer Biofilm" = "Sewer\nBiofilm"))

#Colors for Plot

# Assign colors: first 15 get a palette, "Other" gets grey
ref = "Other"
myColors <- fortyfive_pal[1:length(levels(plot_df$class))]
names(myColors) <- levels(plot_df$class)
myColors[names(myColors)==ref] <- "grey"


##actual plotting
library(dplyr)
library(ggplot2)
library(ggh4x)
library(ggtext)   # for element_markdown


# italicize legend labels
family_labels <- setNames(
  ifelse(levels(plot_df$class) == "Other",
         "Other",
         paste0("*", levels(plot_df$class), "*")),   # markdown italics
  levels(plot_df$class)
)


# Plot using Event (factor) on x → no gaps within each facet

library(scales)  # for percent_format()

#rename the timepoints to make a shorter label
plot_df_updated  = plot_df %>% 
  mutate(NewDateString = case_when(NewDateString == "1-a-1" ~ "1",
                                   NewDateString == "2-a-1" ~ "2",
                                   NewDateString == "3-a-1" ~ "3",
                                   .default= NewDateString))


plot = ggplot(plot_df_updated, aes(x = NewDateString, y = rel_ab, fill = class)) +
  geom_col(position = "stack") +
  facet_nested(Location ~ Type, scales = "free_x", space = "free_x", switch = "y") +
  scale_fill_manual(values = myColors,
                    labels = family_labels) +
  # scale_y_continuous(
  #   breaks = c(0, 0.25, 0.50, 0.75),   # tick marks at 25, 50, 75
  #   labels = c("0","25", "50", "75"),   # plain numbers, no %
  #   limits = c(0, 1)                 # optional: keep to 100%
  # ) +
  labs(x = "Timepoint", y = "Relative Abundance", fill = "", title = "") +
theme_minimal(base_size = 11) +
  guides(fill = guide_legend(ncol = 1))+
theme(
  panel.background = element_rect(fill = "white", colour = "black"),
  strip.background = element_rect(fill = "white", colour = "black"),
  strip.text = element_text(
    face = "bold", 
    hjust = 0.5,          # center horizontally
    vjust = 0.5,          # center vertically
    lineheight = 0.9      # tighter line spacing if wrapped
  ),
  #axis.text.x        = element_blank(),
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 6),
  panel.spacing      = unit(0.03, "lines"),
  legend.position    = "right",
  legend.key.size    = unit(0.3, "cm"),
  legend.text        = ggtext::element_markdown(size = 5),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  #panel.grid.major.y = element_line(color = "grey80"),
  panel.grid.minor.y = element_blank(),
  panel.ontop        = FALSE   # <- ensures geoms (bars) are drawn over grid lines
)


print(plot)

ggsave(path = "Biofilm Project Figures",file="16Sstacked bar plot of families detected-only sewer biofilm.png", plot=plot, width=3, height=6,dpi = 400)
  ggsave(path = "Biofilm Project Figures",file="16S stacked bar plot of families detected-only sewer biofilm.svg", plot=plot, width=3, height=6)
  

  
  #-------------Create summary table of percent contribution by pathogen and sample type------------------
family_summary <- plot_df %>%
  group_by(class, Type) %>%
  summarise(
    MinPercent = min(rel_ab, na.rm = TRUE) * 100,
    MaxPercent = max(rel_ab, na.rm = TRUE) * 100,
    MeanPercent = mean(rel_ab, na.rm = TRUE) * 100,
    MedianPercent = median(rel_ab, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  dplyr::rename(Family = class, SampleType = Type) %>%
  arrange(SampleType, desc(MeanPercent))


```







##NMDS and PermANOVA
###Microbiome NMDS only on samples with OTU -level call and PERMANOVA and Explained Variance
```{r}

library(reshape2)

df_long <- psmelt(ps_pruned) %>%  # creates data.frame with OTU abundances + metadata

  dplyr::select(OTU,Sample,Abundance,sample_id,sample_type,sample_date,timepoint,corresponding_sewer,sub_timepoint,qubit_DNA_ng.uL,Kingdom,Phylum,Class, Order, Family, Genus, Species) %>% distinct %>% 
  dplyr::rename("SampleID" = "sample_id",
                "Type" = "sample_type",
                "Location" = "corresponding_sewer",
                "Date" = "sample_date") %>% 
  mutate(Date = as.Date(as.character(Date), format = "%m/%d/%y"),
        # Make Date a string in MM/DD/YYYY format
        Abundance = as.numeric(Abundance),
        DateString = format(as.Date(Date), "%m/%d/%Y"),
        Month = format(as.Date(Date), "%m-%Y")) %>% 
  group_by(SampleID,Type) %>%              # group by sample
  mutate(rel_ab = Abundance / sum(Abundance)) %>%  # convert counts to relative abundance (0-1)
  ungroup() %>% 
  
  dplyr::select(-Abundance) 


metadata = df_long %>%
  dplyr::select(SampleID,Date,Location,Type) %>% distinct() %>% 
  column_to_rownames(var = "SampleID")

# 1) Extract species-level rows (where clade_name contains "s__" OR length of pipe-split >=7)
df <- df_long 

# make dataframe with a column per sample and values are rpkm
transposed = df %>% dplyr::select(SampleID,OTU,rel_ab) %>% 
  
  #as numeric
  mutate(rel_ab = as.numeric(rel_ab)) %>% 
  
  group_by(SampleID, OTU) %>% 
  summarize(rel_ab = sum(rel_ab, na.rm = TRUE), .groups = "drop") %>%
  ungroup() %>% 
  
  pivot_wider(
    names_from = SampleID,   # new column names come from sample_id
    values_from = rel_ab        # values in cells come from rpkm
  ) 
  
  #replace NAs with zeros
  #mutate_all( ~replace(., lengths(.)==0, 0)) 


#store data in counts variable and only show fraction reads 
counts = transposed %>% dplyr::select(-OTU) %>% t() 

# Extract row names as a separate vector
row_names <- rownames(counts)


#remove samples
#counts <- counts[!(rownames(counts) %in% c("S209", "S210","S1","S2","S3")), ]

#make sure it's a numeric dataframe
counts = as.data.frame(counts) %>% mutate(across(everything(), ~ as.numeric(as.character(.))))

#replace NAs with zero
counts[is.na(counts)] = 0
                       
# If too few samples remain, think twice
cat("Samples with species-level data:", nrow(counts), "\n")   

# Optional taxon filtering: remove very rare species
min_prevalence <- 0.005  # species present in >=5% samples
present_in <- colSums(counts > 0)
keep_taxa <- present_in >= (min_prevalence * nrow(counts))
comm_sp_filt <- counts[, keep_taxa, drop = FALSE]

# Hellinger transform (recommended for species abundance)
comm_sp_hell <- decostand(comm_sp_filt, method = "hellinger")

#get the distance matrix
dist1 <- vegdist(counts,method = "bray")
#dist1 = vegdist(comm_sp_hell, method = "bray") # if we want Hellinger transformation and removal of rare species
argMDS <- isoMDS(dist1, trace = F)
scores <- as.data.frame(argMDS)
scores <- scores[,1:2]

# Assuming 'scores' and 'metadata' are your dataframes
common_rows <- intersect(rownames(counts), rownames(metadata))

# Filter 'scores' dataframe to keep only rows present in 'metadata'
filtered_counts <- counts[rownames(counts) %in% common_rows, ]


# Filter 'scores' dataframe to keep only rows present in 'metadata'
filtered_metadata <- metadata[rownames(metadata) %in% common_rows, ]



# meta MDS and create envfit variable
vare.mds <- metaMDS(filtered_counts, trace = FALSE, distance = "bray")
#library(vegan)
ef <- envfit(vare.mds, filtered_metadata, permu = 999,na.rm = TRUE)


#plot results
plot(vare.mds, display = "sites")


#extract envfit arrows data and filter out non significant one 
#en_coord_cont = as.data.frame(scores(ef, "vectors")) * ordiArrowMul(ef) 
#en_coord_cont$pval <- ef[["vectors"]][["pvals"]]
#en_coord_cont <- filter(en_coord_cont,pval<=0.10)

#extract nmds scores and associate to metadata 
nmds_scores <- as.data.frame(scores(vare.mds)$sites)

# Assuming your_data is your dataframe
nmds_scores$SampleID <- rownames(nmds_scores)
rownames(nmds_scores) <- NULL # now get rid of row names

full_metadata = metadata %>% rownames_to_column(var = "SampleID") %>% 
  
  #create new variable combined with site and sample type
  mutate("site_type"=paste0(Location,"_",Type))

nmdsMerged = left_join(nmds_scores,full_metadata, by = "SampleID") 


#specify subareas we want to incude in plot
locations = c("FIONA","SHREK","OSCAR","MARIO","LUIGI")


nmds_plotg <- nmdsMerged %>% 
  
  #filter locations if needed
  #filter(Location %in% locations) %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = site_type)) +
  geom_point() +
  stat_ellipse(linewidth = 0.75) +
  
  #change colors
  #scale_color_manual(values = c(
    #"4" = "lightsalmon", "12" = "grey",
   #                             )) +
  #scale_shape_manual(values = c(1, 2
                                #, 3, 4, 5, 6, 16, 17, 18, 19, 20, 21,22, 29
                                #)) +
  theme_bw() +
  labs(
    title = "NMDS Plot",
    x = "NMDS1",
    y = "NMDS2",
    color = "site_type",
    size = "Collection Date"
  ) + 
  #geom_vline(xintercept = 0,linetype = "dashed")+ geom_hline(yintercept = 0,linetype = "dashed")+
  #geom_segment(aes(x = 0, y = 0, xend =NMDS1, yend = NMDS2), data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  #geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", fontface = "bold", label = row.names(en_coord_cont)) +
  theme(axis.text.y = element_text(angle = 45, hjust = 1)) + theme(axis.text.y = element_blank())

# print plot
nmds_plotg

#prettier NMDS plot with shapes and colors
nmds_plotg <- nmdsMerged %>% 
  
  #filter out municipal ww
  filter(Type != "Mu_WW") %>% 
  
  #rename some variables
  mutate(Type = case_when(Type == "Mu_WW"~ "Municipal WW",
                          Type == "w" ~ "Wastewater",
                          Type == "tap" ~ "Tap Water",
                          .default = Type)) %>% 
  
  #filter locations if needed
  #filter(Location %in% locations) %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = Type, shape = Location, group = Type)) +
  geom_point() +
  stat_ellipse(linewidth = 0.75) +
  
  #change colors
  #scale_color_manual(values = c(
    #"4" = "lightsalmon", "12" = "grey",
   #                             )) +
  #scale_shape_manual(values = c(1, 2
                                #, 3, 4, 5, 6, 16, 17, 18, 19, 20, 21,22, 29
                                #)) +
  theme_bw() +
  labs(
    title = "NMDS Plot",
    x = "NMDS1",
    y = "NMDS2",
    color = "site_type",
    size = "Collection Date"
  ) + 
  
  #manualy change colors
  scale_color_manual(values = c("d" = fortyfive_pal[[1]], "e" = fortyfive_pal[[2]], 
                                "Municipal WW" = fortyfive_pal[[3]], 
                                "Wastewater" = fortyfive_pal[[4]],"Tap Water" = "grey76"
                                )) +
  #scale_shape_manual(values = c("Drain" = 1, "Endcap" = 2, "Municipal WW" = 3,  "Hospital WW" = 4 )) +
  #geom_vline(xintercept = 0,linetype = "dashed")+ geom_hline(yintercept = 0,linetype = "dashed")+
  #geom_segment(aes(x = 0, y = 0, xend =NMDS1, yend = NMDS2), data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  #geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", fontface = "bold", label = row.names(en_coord_cont)) +
  labs(title = "",
       color = "",
       shape = "") +
  theme(axis.text.y = element_text(angle = 45, hjust = 1)) + theme(axis.text.y = element_blank())

# print plot
nmds_plotg



#prettier NMDS plot with shapes and colors
nmds_plotg <- nmdsMerged %>% 
  
  #filter out municipal ww
  filter(Type != "Mu_WW") %>% 
  
  #rename some variables
  mutate(Type = case_when(Type == "Mu_WW"~ "Municipal WW",
                          Type == "e" ~ "Branch Biofilm",
                          Type == "d" ~ "Sink Biofilm",
                          Type == "w" ~ "Wastewater",
                          Type == "tap" ~ "Tap Water",
                          .default = Type)) %>% 
  
  #filter locations if needed
  #filter(Location %in% locations) %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = Type, shape = Location, group = Type)) +
  geom_point() +
  stat_ellipse(linewidth = 0.75) +
  
  #change colors
  #scale_color_manual(values = c(
    #"4" = "lightsalmon", "12" = "grey",
   #                             )) +
  #scale_shape_manual(values = c(1, 2
                                #, 3, 4, 5, 6, 16, 17, 18, 19, 20, 21,22, 29
                                #)) +
  theme_bw() +
  labs(
    title = "NMDS Plot",
    x = "NMDS1",
    y = "NMDS2",
    color = "site_type",
    size = "Collection Date"
  ) + 
  
  #manualy change colors
  scale_color_manual(values = c("Sink Biofilm" = "#79307D", "Branch Biofilm" = "#417C8C", 
                                "Municipal WW" = fortyfive_pal[[3]], 
                                "Wastewater" = "#E57262","Tap Water" = "grey76"
                                )) +
  #scale_shape_manual(values = c("Drain" = 1, "Endcap" = 2, "Municipal WW" = 3,  "Hospital WW" = 4 )) +
  #geom_vline(xintercept = 0,linetype = "dashed")+ geom_hline(yintercept = 0,linetype = "dashed")+
  #geom_segment(aes(x = 0, y = 0, xend =NMDS1, yend = NMDS2), data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  #geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", fontface = "bold", label = row.names(en_coord_cont)) +
  labs(title = "Kraken2 NMDS Analysis",
       color = "",
       shape = "") +
  theme(axis.text.y = element_text(angle = 45, hjust = 1)) + 
  theme(axis.text.y = element_blank(),
        axis.title = element_text(size= 20,face = "bold",color = "black"),
        legend.text = element_text(size = 15,face = "bold", color = "black"))

# print plot
nmds_plotg

# Save the combined plot as a PNG file
ggsave(plot = nmds_plotg, path = "Biofilm Project Figures", "16S- Microbiome Bracken derived NMDS (species level) by Site and Location.png", dpi = 300, height = 5, width = 7.5, units = "in")

ggsave(plot = nmds_plotg, path = "Biofilm Project Figures", "16S- Microbiome Bracken derived NMDS (species level) by Site and Location.svg", height = 6, width = 7)


#----------Run PERMANOVA---------------------


#specify subareas we want to incude
#locations = c("FIONA","SHREK","OSCAR","MARIO","LUIGI")

# Subset the matrix based on the selected rows (so filtering for specific sites)
subset_metadata <- filtered_metadata 

  # 
  # #create new variable combined with site and sample type
  # mutate("site_type"=paste0(Location,"_",Type))


# Get row names that are present in both data frames
common_rows <- intersect(rownames(subset_metadata), rownames(filtered_counts))


subset_adonis = filtered_counts[common_rows,]

subset_metadata <- subset_metadata[common_rows, ]
subset_adonis   <- filtered_counts[common_rows, ]
stopifnot(identical(rownames(subset_adonis), rownames(subset_metadata)))


#by = "terms"
# Reports sequential sums of squares: each term’s contribution is calculated in order.
# Gives you all components (Location, Type, Location:Type) plus residual, so the totals add to 100% of variation.

# Caveat: order matters, but usually the main effects are reported first, then interaction, which is intuitive


# Run PERMANOVA
res <- adonis2(subset_adonis ~ Location * Type, subset_metadata,
               permutations = 9999, by = "terms")

# Convert to data frame
res_df <- as.data.frame(res) %>% 
  # Add a column for term names
  mutate(Term = rownames(res),
         # Replace very small p-values with <0.0001 for readability
         `p-value` = ifelse(`Pr(>F)` < 0.0001, "<0.0001", `Pr(>F)`),
         # Round numeric columns nicely
         R2 = round(R2, 4),
         F = round(F, 3),
         SumOfSqs = round(SumOfSqs, 3)) %>%
  # Select and reorder columns
  dplyr::select(Term, Df, SumOfSqs, R2, F, `p-value`)

# View the data frame
res_df


library(dplyr)
library(ggplot2)
library(scales)

df_plot <- res_df %>%
  filter(Term != "Total") %>%        # remove Total row
  mutate(Percent = R2)     %>%           # R2 already proportions
  mutate(Term = case_when(Term == "Residual" ~ "Unknown",
                          .default = Term))

df_plot$Term <- factor(df_plot$Term,
                       levels = c("Location", "Type", "Location:Type", "Unknown"))


plot = ggplot(df_plot, aes(x = "Explained Variation", y = Percent, fill = Term)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = NULL,
       y = "Explained Variation (%)") +
  theme_classic()


# Save with ggsave
ggsave("Biofilm Project Figures/16S Explained Permanova Variation.svg", plot = plot, width = 4, height = 6)

# Save with ggsave
ggsave("Biofilm Project Figures/16S Explained Permanova Variation.png", plot = plot, width = 4, height = 6,dpi = 300)


#run pairwise adonis using function
pairwise_adonis_microbiome = pairwise.adonis2(subset_adonis, subset_metadata$Location, sim.method = 'bray', p.adjust.m = 'bonferroni'); pairwise_adonis_microbiome

#adonis.pair(vegdist(subset_adonis),subset_metadata$subarea)


#run pairwise adonis using function
pairwise_adonis_microbiome = pairwise.adonis2(subset_adonis, subset_metadata$Type, sim.method = 'bray', p.adjust.m = 'bonferroni'); pairwise_adonis_microbiome

subset_adonis16S= subset_adonis

subset_metadata16S= subset_metadata

#adonis.pair(vegdist(subset_adonis),subset_metadata$subarea)



```



##16S Venn Diagrams
###Venn diagram for 16S
```{r}

# ---- Prepare genus-level matrix ----
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")


library(reshape2)




ps_genus <- ps_pruned %>%
  microbiome::aggregate_taxa("Genus")

# Calculate compositional version of the data
# (relative abundances)

pseq.rel <- microbiome::transform(ps_genus, "compositional")

df_long <- psmelt(pseq.rel) %>%  # creates data.frame with OTU abundances + metadata

  dplyr::select(OTU,Sample,Abundance,sample_id,sample_type,sample_date,timepoint,corresponding_sewer,sub_timepoint,qubit_DNA_ng.uL,Kingdom,Phylum,Class, Order, Family, Genus) %>% distinct %>% 
  dplyr::rename("SampleID" = "sample_id",
                "Type" = "sample_type",
                "Location" = "corresponding_sewer",
                "Date" = "sample_date") %>% 
  mutate(Date = as.Date(as.character(Date), format = "%m/%d/%y"),
        # Make Date a string in MM/DD/YYYY format
        Abundance = as.numeric(Abundance),
        DateString = format(as.Date(Date), "%m/%d/%Y"),
        Month = format(as.Date(Date), "%m-%Y")) %>% 

  mutate(rel_ab = Abundance ) %>%  # convert counts to relative abundance (0-1)

  
  #filter for rows that have Genus
  drop_na(Genus)



library(dplyr)
library(ggvenn)

# Colors
earthy_green  <- "#9FB498"
earthy_purple <- "#8488AC"
earthy_orange <- "#E57262"

viridisCol1 = "#79307D"
viridisCol2 = "#417C8C"
viridsCol3 = "#E57262"

# Subareas to include
locations <- c("FIONA","SHREK","OSCAR","MARIO","LUIGI")



#overlap for total pathogens across all Sites------------------------


# Species-level rows only
venn_df <- df_long 


  df <- df_long %>%
  distinct(Type, SampleID, Genus, rel_ab) %>%
  filter(rel_ab != 0) %>%
  dplyr::select(Type, Genus) %>%
  distinct()
  
  sets_total <- list(
  Sink_Biofilm   = df %>% filter(Type == "d") %>% pull(Genus),
  Sewer_Biofilm = df %>% filter(Type == "e") %>% pull(Genus),
  Tap = df %>% filter(Type == "tap") %>% pull(Genus),
  Hospital_WW   = df %>% filter(Type == "w") %>% pull(Genus)
)

all_species <- unique(unlist(sets_total))

presence_mat <- sapply(sets_total, function(s) all_species %in% s)

region_df <- tibble(
  region = apply(presence_mat, 1, function(x)
    paste(which(x), collapse = "_"))
) %>%
  dplyr::count(region, name = "n")

  
print(region_df)


library(eulerr)

fit_total <- euler(c(
  "Sink_Biofilm"   = region_df$n[region_df$region == "1"],
  "Sewer_Biofilm" = region_df$n[region_df$region == "2"],
  "Hospital_WW"   = region_df$n[region_df$region == "3"],

  "Sink_Biofilm&Sewer_Biofilm" =
    region_df$n[region_df$region == "1_2"],

  "Sink_Biofilm&Hospital_WW" =
    region_df$n[region_df$region == "1_3"],

  "Sewer_Biofilm&Hospital_WW" =
    region_df$n[region_df$region == "2_3"],

  "Sink_Biofilm&Sewer_Biofilm&Hospital_WW" =
    region_df$n[region_df$region == "1_2_3"]
))

plot(
  fit_total,
  fills = c(viridisCol1, viridisCol2, viridsCol3),
  edges = TRUE,
  labels = list(font = 2, cex = 0.9),
  quantities = list(type = "counts", font = 2, cex = 0.9)
)

# Capture whatever is currently on the device
g <- grid.grab()

# Save with ggsave
ggsave("Biofilm Project Figures/16S Based Euler diagram Total -Counts -ONLY SHARED DATES WITH BIOFILM-RPIP only.svg", plot = g, width = 6, height = 6)

# Save with ggsave
ggsave("Biofilm Project Figures/16S Based Euler diagram Total -Counts -ONLY SHARED DATES WITH BIOFILM-RPIP only.png", plot = g, width = 6, height = 6,dpi = 300)

#percentages
plot(
  fit_total,
  fills = c(viridisCol1, viridisCol2, viridsCol3),
  edges = TRUE,
  labels = list(font = 2, cex = 0.9),
  quantities = list(type = "percent", font = 2, cex = 0.9)
)

# Capture whatever is currently on the device
g <- grid.grab()

# Save with ggsave
ggsave("Biofilm Project Figures/16S Based Euler diagram Total -Percents -ONLY SHARED DATES WITH BIOFILM-RPIP only.svg", plot = g, width = 6, height = 6)

# Save with ggsave
ggsave("Biofilm Project Figures/16S Based Euler diagram Total -Percents -ONLY SHARED DATES WITH BIOFILM-RPIP only-percents.png", plot = g, width = 6, height = 6,dpi = 300)

  
  # Build Venn input
  ggvenn_df <- list(
  Sink_Biofilm   = df %>% filter(Type == "d") %>% pull(Genus),
  Sewer_Biofilm = df %>% filter(Type == "e") %>% pull(Genus),
  Tap = df %>% filter(Type == "tap") %>% pull(Genus),
  Wastewater   = df %>% filter(Type == "w") %>% pull(Genus)
  )

  # Plot with a title for clarity
  print(
    ggvenn(
      ggvenn_df,
      fill_color = c(viridisCol1, viridisCol2, viridsCol3),
      stroke_size = 0.5,
      set_name_size = 4
    ) +
    ggtitle("Genera overlap")
  )

```



#Dissimilarity Indices by Sample Type
###BrayCurtis distance by sample type
```{r}
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(ggstatsplot)

#specify dates
dates = c("2024-09-09", "2024-09-11", "2024-09-12", "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09")


# ---- Prepare genus-level matrix ----
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")



library(reshape2)

df_long <- psmelt(ps) %>%  # creates data.frame with OTU abundances + metadata

  dplyr::select(OTU,Sample,Abundance,sample_id,sample_type,sample_date,timepoint,corresponding_sewer,sub_timepoint,qubit_DNA_ng.uL,Kingdom,Phylum,Class, Order, Family, Genus, Species) %>% distinct %>% 
  dplyr::rename("SampleID" = "sample_id",
                "Type" = "sample_type",
                "Location" = "corresponding_sewer",
                "Date" = "sample_date") %>% 
  mutate(Date = as.Date(as.character(Date), format = "%m/%d/%y"),
        # Make Date a string in MM/DD/YYYY format
        Abundance = as.numeric(Abundance),
        DateString = format(as.Date(Date), "%m/%d/%Y"),
        Month = format(as.Date(Date), "%m-%Y")) %>% 
  group_by(SampleID,Type) %>%              # group by sample
  mutate(rel_ab = Abundance / sum(Abundance)) %>%  # convert counts to relative abundance (0-1)
  ungroup() %>% 
  
  dplyr::select(-Abundance)


# Create sample x species matrix
pathogen_mat <- df_long %>%
  dplyr::select(SampleID, OTU, rel_ab) %>% distinct() %>% 
  pivot_wider(names_from = SampleID, values_from = rel_ab, values_fill = 0) %>% 
  column_to_rownames("OTU") %>%
  t() %>%
  as.data.frame()

# ---- Compute Bray–Curtis per Location × Type ----
bc_long <- df_long %>%
  filter(Type %in% c("d", "e", "w","tap"),
         Location %in% locations) %>%
  group_by(Location, Type) %>%
  group_modify(~{
    df_sub <- .x
    ids <- intersect(df_sub$SampleID, rownames(pathogen_mat))
    
    # Need at least 2 samples
    if(length(ids) < 2){
      return(tibble(BrayCurtis = numeric(0)))  
    }
    
    # Compute Bray-Curtis
    d <- vegdist(pathogen_mat[ids, , drop = FALSE], method = "bray")
    d_table <- as.data.frame(as.table(as.matrix(d)))
    colnames(d_table) <- c("Sample1", "Sample2", "BrayCurtis")
    
    # Keep only non-self comparisons
    d_table <- d_table %>% filter(Sample1 != Sample2)
    
    # Return as a tibble (data frame)
    tibble(BrayCurtis = d_table$BrayCurtis)
  }) %>%
  ungroup() %>%
  # Add back grouping variables safely
  mutate(Type = case_when(
           Type == "w" ~ "Wastewater",
           Type == "e" ~ "Sewer Biofilm",
           Type == "d" ~ "Sink Biofilm",
           TRUE ~ as.character(Type)
         ))

# ---- Reorder factor by median ----
medians <- bc_long %>% group_by(Type) %>% summarise(median_bc = median(BrayCurtis))
bc_long$Type <- factor(bc_long$Type, levels = medians$Type[order(medians$median_bc)])

# ---- Plot faceted by Location ----
Location_braycurtis_plot <- ggplot(bc_long, aes(x = Type, y = BrayCurtis)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5) +
  ylab("Bray–Curtis Dissimilarity") +
  xlab("") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 9, face = "bold", color = "black")) +
  facet_wrap(~Location)

print(Location_braycurtis_plot)

# ---- Plot with facet by Location ----
Location_braycurtis_plot <- ggbetweenstats(
  data = bc_long,
  x = Type,
  y = BrayCurtis,
  outlier.tagging = TRUE,
  type = "nonparametric",
  mean.plotting = TRUE,
  mean.ci = FALSE
) + 
  ylab("Bray–Curtis Dissimilarity") +
  xlab("") +
  theme_classic() +
  theme(axis.text.x= element_text(size = 9, face = "bold", color= "black")) 

print(Location_braycurtis_plot)



# Save the combined plot as a PNG file
ggsave(plot = Location_braycurtis_plot,path = "Biofilm Project Figures", "16S: Bray-curtis dissimilarity-ONLY SHARED DATES WITH BIOFILM-RPIP only.png", dpi = 300, height = 6, width = 7, units = "in")



```


#Heat Maps
##RPIP only Relative Abundance Heat Maps 
###set themes
```{r}
library(ggplot2)

# Reset global theme and font
theme_set(theme_grey(base_family = "", base_size = 11))

# Remove theme_update() overrides
theme_update(
  axis.text.x        = element_text(),
  axis.text.y        = element_text(),
  axis.ticks.x       = element_line(),
  axis.ticks.y       = element_line(),
  axis.ticks.length.x = unit(0.15, "lines"),
  panel.grid         = element_line(),
  plot.background    = element_rect()
)

# Reset geom defaults that may have inherited "Chivo"
update_geom_defaults("text",  list(family = ""))
update_geom_defaults("label", list(family = ""))

showtext::showtext_auto(FALSE)

```
### (RPIP only) Top Pathogens per Site median Rel Abundance Heat Map (faceted with locations stacked) (x axis is dates) - ONLY SHARED DATES WITH BIOFILM BUT SAME ORDER AS IF ALL SAMPLES WERE INCLUDED-ordered by top pathogens by median concentration across all Sites in sewer biofilm- count NAs as zero using complete
```{r}
#specify dates
dates = c("2024-09-09", "2024-09-11", "2024-09-12", "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09")

# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame


# RPIP taxa
rpip_targets <-  rpip_pathogen_species_list


# ---------- STEP 1: Extract and clean ----------
heatmap_df <- kraken_merged_df_filtered %>%
  filter(
    taxonomy_lvl == "S",
    Location %in% locations,
    name %in% rpip_targets,
    Type %in% c("Drain", "Endcap", "H_WW")) %>% 
  mutate(
    rel_ab = as.numeric(fraction_total_reads),
    Date   = as.Date(Date),
    Type   = factor(Type, levels = c("Drain", "Endcap", "H_WW"))
  ) %>% 
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>%
  distinct() %>%
  
  mutate(
    Type = factor(Type, levels = c("Drain", "Endcap", "H_WW")),
    TypeLabel = case_when(
      Type == "H_WW" ~ "Wastewater", 
      TRUE ~ as.character(Type)
    ),
    # Make Date a string in MM/DD/YYYY format
    DateString = format(as.Date(Date), "%m/%d/%Y")
  )



# # ---------- STEP 2: Assign letters for drain replicates ----------
# drain_samples <- heatmap_df %>%
#   filter(Type == "Drain") %>%
#   distinct(Location, DateString, Type,SampleID) %>% 
#   group_by(Location, DateString) %>%
#   mutate(SampleLetter = letters[seq_len(n())],
#          NewDateString = paste0(DateString, "-", SampleLetter)) %>%
#   ungroup() %>%
#   dplyr::select(SampleID,Location, DateString, Type, NewDateString)
# 
# # ---------- STEP 2: Assign letters for sewer biofilm replicates ----------
# sewer_samples <- heatmap_df %>%
#   filter(Type == "Endcap") %>%
#   distinct(Location, DateString, Type, SampleID) %>% 
#   mutate(
#     # Convert to Date first
#     Date_parsed = as.Date(DateString, format = "%m/%d/%Y"),
#     # Use case_when consistently
#     NewDateString = case_when(
#       format(Date_parsed, "%m") == "08" ~ "08/2024",
#       TRUE ~ DateString
#     )
#   ) %>% 
#   
#   # group_by(Location) %>%
#   # mutate(SampleLetter = letters[seq_len(n())],
#   #        NewDateString = paste0(DateString, "-", SampleLetter)) %>%
#   # ungroup() %>%
#   dplyr::select(SampleID,Location, DateString, Type, NewDateString)
# 
# # Combine drain and sewer samples
# replicate_samples <- bind_rows(
#   drain_samples,
#   sewer_samples %>% dplyr::select(SampleID, Location, DateString, Type, NewDateString)
# )
# 
# # ---------- Join back to main dataframe ----------
# heatmap_df_update <- heatmap_df %>%
#   left_join(replicate_samples, by = c("SampleID", "Location", "DateString", "Type")) %>% 
#   mutate(
#     NewDateString = ifelse(is.na(NewDateString) == TRUE, DateString, NewDateString)
#   )

# ---------- New Way: Assign letters time points ----------
drain_samples <- heatmap_df %>%
  filter(Type == "Drain") %>%
  distinct(Location, DateString, Type,SampleID,timepoint) %>% 
  group_by(Location, DateString) %>%
  mutate(SampleLetter = letters[seq_len(n())],
         NewDateString = paste0(timepoint, "-", SampleLetter)) %>%
  ungroup() %>%
  dplyr::select(SampleID,Location, DateString, Type, NewDateString)

# ---------- Assign letters for sewer biofilm replicates ----------
sewer_samples  <- heatmap_df %>%
  filter(Type == "Endcap") %>%
  distinct(Location, DateString, Type,SampleID,timepoint) %>% 
  group_by(Location, DateString) %>%
  mutate(SampleLetter = letters[seq_len(n())],
         #NewDateString = paste0(timepoint, "-", SampleLetter)) %>%
         NewDateString = timepoint) %>%
  ungroup() %>%
  dplyr::select(SampleID,Location, DateString, Type, NewDateString)

# ---------- Assign letters for ww replicates ----------
ww_samples  <- heatmap_df %>%
  filter(Type == "H_WW") %>%
  distinct(Location, DateString, Type,SampleID,timepoint) %>% 
  group_by(Location, DateString) %>%
  mutate(SampleLetter = letters[seq_len(n())],
         #NewDateString = paste0(timepoint, "-", SampleLetter)) %>%
          NewDateString = timepoint) %>%
  ungroup() %>%
  dplyr::select(SampleID,Location, DateString, Type, NewDateString)

# Combine drain and sewer samples
replicate_samples <- bind_rows(
  drain_samples,
  ww_samples,
  sewer_samples %>% dplyr::select(SampleID, Location, DateString, Type, NewDateString)
)

# ---------- Join back to main dataframe ----------
heatmap_df_update <- heatmap_df %>%
  left_join(replicate_samples, by = c("SampleID", "Location", "DateString", "Type")) %>%
  mutate(
    NewDateString = ifelse(is.na(NewDateString) == TRUE, DateString, NewDateString)
  )


# ---------- Filter for top species ----------
plot_df <- heatmap_df_update %>%
    mutate(

    Species_italic = paste0("*", name, "*")  # markdown italics
  ) %>% 
  
  mutate(TypeLabel = case_when(TypeLabel == "Endcap" ~ "Sewer\nBiofilm",
                               TypeLabel == "Drain" ~ "Sink\nBiofilm",
                               .default = TypeLabel)) %>% 
  
  mutate(
    TypeLabel = factor(
      TypeLabel,
      levels = c("Wastewater", "Sewer\nBiofilm","Sink\nBiofilm")
    )
  ) %>% 
  
  mutate(
    FillValue = log10(rel_ab)   # keep continuous relative abundance
  )

#Filter top per Site--------Pick top species by total relative abundance ----------
n_species = 10


top_species_df <- plot_df %>%
  filter(TypeLabel == "Sewer\nBiofilm") %>% 
  
  # complete dataframe so NAs are zero
  # This ensures every SampleID has a row for every Species_italic
  tidyr::complete(SampleID, Species_italic, fill = list(rel_ab = 0)) %>%
  
  group_by(Species_italic) %>%
  summarise(median_abundance = median(rel_ab, na.rm = TRUE), .groups = "drop") %>%
  slice_max(order_by = median_abundance, n = n_species) %>%
  arrange(desc(median_abundance)) %>%
  ungroup()

top_species = top_species_df %>% dplyr::select(Species_italic) %>% distinct() %>% pull()

#reorder species to have top per Site
plot_df_top <- plot_df %>%
  semi_join(top_species_df, by = c( "Species_italic")) %>% 
  # mutate(
  #   Species_italic = reorder_within(
  #     Species_italic,
  #     rel_ab,
  #     Location,
  #     fun = median
  #   )
  # ) %>% 
  left_join(
    top_species_df %>% dplyr::select(Species_italic, median_abundance),
    by = "Species_italic"
  ) %>%

  mutate(
    Species_italic = forcats::fct_reorder(
      Species_italic,
      median_abundance,
      .desc = FALSE
    )
  ) %>%
  
  #filter dates
  filter(
    (TypeLabel == "Wastewater" & Date %in% dates) |
    (TypeLabel != "Wastewater")
  ) %>% 
  
  # now replace Site names 
  mutate(Location = case_when(Location == "MARIO" ~ "Site 1",
                              Location == "FIONA" ~ "Site 2",
                              Location == "LUIGI" ~ "Site 3",
                              Location == "SHREK" ~ "Site 4",
                              Location == "OSCAR" ~ "Site 5")) %>% 
  
  drop_na(timepoint) %>%
  mutate(
    Species_italic = forcats::fct_drop(Species_italic)
  )  %>%
  group_by(TypeLabel) %>%
  tidyr::complete(
    Species_italic,
    NewDateString,
    Location,
    fill = list(FillValue = NA)
  ) %>%
  ungroup()

# ---------- STEP 4: Plot as abundance heatmap ----------
fortyfive_pal = c("#a0cb6b","#8368cb","#c86c69","#cdd3e5","#dab594","#d692d1",
               "#7495c3","#9fdeca","#e2e8c3","#d8a4af","#71bed2","#bca9dd",
               "#8bb598","#e5cbd4","#6d7ecd","#e4d5ca","#a8dfa5","#a0bada",
               "#cbca6f","#c6926f","#cce7e6","#81a4b2","#ca69c3","#76bd75",
               "#d37dae","#abb684","#6ecda4","#a88bbb","#8cbbb6","#d18191",
               "#c1d8c2","#dac4e6","#e0adce","#a96dca","#e2bdb6","#aacedd",
               "#9992d9","#e0d0ab","#a2abdf","#b88b9e","#b4a382","#dba294",
               "#c9dd9b","#c8b67c","#b9847b")


med_col  <- "#a0cb6b"
high_col <- "#8368cb"

med_col  <- "#e2bdb6"   # soft peach
high_col <- "#6d7ecd"   # periwinkle blue

med_col  <- "#9fdeca"   # light teal
high_col <- "#084081"   # keep your dark anchor OR use "#6d7ecd"
high_col = "#6d7ecd"
low_col = "#f7fcf0"

heat_map <- ggplot(plot_df_top, aes(x = NewDateString, y = Species_italic, fill = FillValue)) +
  geom_raster() +
  # scale_fill_gradient(low = "white", high = "red", 
  #                     na.value = "white",   # 🔑 white for NA
  #                     name = "Log10(Relative Abundance)") +
  
  scale_fill_gradientn(
    colors = c(low_col, med_col, high_col),
    values = scales::rescale(c(
      min(plot_df_top$FillValue, na.rm = TRUE),
      median(plot_df_top$FillValue, na.rm = TRUE),
      max(plot_df_top$FillValue, na.rm = TRUE)
    )),
    na.value = "white",
    name = "Log10(Relative Abundance)"
  ) +
  
  facet_grid(Location ~ TypeLabel, 
             scales = "free", 
             space = "free", 
             switch = "y",
             #scales = "free_x",
            #space  = "free_x",
             labeller = labeller(TypeLabel = label_wrap_gen(width = 8))
            ) +
  scale_x_discrete(
    labels = plot_df %>%
      distinct(SampleID, Date, TypeLabel) %>%
      arrange(TypeLabel, Date) %>%
      mutate(Label = paste0(Date)) %>%
      dplyr::select(SampleID, Label) %>%
      tibble::deframe()
  ) +
  #scale_y_reordered() +   # 🔑 THIS removes ___SiteNAME
  #labs(x = "", y = "", title = paste0("Top ", n_species, " Species by Summed Relative Abundance Across Dates and Sample Types")) +
  labs(x = "Timepoint", y = "", title = "") +
  theme(
    axis.text.x = element_markdown(angle = 30, size = 7, vjust = 0.9, hjust = 0.9),
    panel.grid = element_blank(),
    axis.text.y = element_markdown(color = "black"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.position = "bottom",
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(size = 10, face = "bold", color = "black")
  )


# Save final heatmap
ggsave(
  plot = heat_map,
  filename = paste0("Biofilm Project Figures/Kraken2 CT 0.5 Top ", n_species, " Species by Median Relative Abundance Across Dates and Sites only Sewer Biofilm -ONLY SHARED DATES WITH BIOFILM.png"),
  width = 5,
  height = 8,
  units = "in"
)

#install.packages("svglite")
library("svglite")

# Save final heatmap
ggsave(
  plot = heat_map,
  filename = "Biofilm Project Figures/Kraken CT 0.5 -RPIP only Kraken Species by Median Relative Abundance Across Dates and Sites only Sewer Biofilm-ONLY SHARED DATES WITH BIOFILM.svg",
  width = 5,
  height = 8
)



plot(heat_map)

```





#CVs and Box Plots
## (RPIP only) Top Pathogens by Median Rel Abundance in Sewer Biofilm Rel Abundance and CVs Ordered by Top Pathogens Overall in Sewer Biofilm (faceted with locations stacked) 
```{r}

# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame


# RPIP taxa
rpip_targets <-  rpip_pathogen_species_list


# ---------- STEP 1: Extract and clean ----------
heatmap_df <- kraken_merged_df_filtered %>%
  filter(taxonomy_lvl == "S") %>% 
  mutate(rel_ab = as.numeric(fraction_total_reads)) %>% 
  filter(Location %in% locations) %>% 
  filter(name %in% c(rpip_targets)) %>% 
  filter(Type %in% c("Drain", "Endcap", "H_WW")) %>% 
  
  
  mutate(
    Type = factor(Type, levels = c("H_WW", "Endcap","Drain")),
    TypeLabel = case_when(
      Type == "H_WW" ~ "Wastewater", 
      TRUE ~ as.character(Type)
    ),
    # Make Date a string in MM/DD/YYYY format
    DateString = format(as.Date(Date), "%m/%d/%Y")
  ) 




# ---------- STEP 2: Assign letters for drain replicates ----------
drain_samples <- heatmap_df %>%
  filter(Type == "Drain") %>%
  distinct(Location, DateString, Type,SampleID) %>% 
  group_by(Location, DateString) %>%
  mutate(SampleLetter = letters[seq_len(n())],
         NewDateString = paste0(DateString, "-", SampleLetter)) %>%
  ungroup() %>%
  dplyr::select(SampleID,Location, DateString, Type, NewDateString)

# ---------- STEP 2: Assign letters for sewer biofilm replicates ----------
sewer_samples <- heatmap_df %>%
  filter(Type == "Endcap") %>%
  distinct(Location, DateString, Type, SampleID) %>% 
  mutate(
    # Convert to Date first
    Date_parsed = as.Date(DateString, format = "%m/%d/%Y"),
    # Use case_when consistently
    NewDateString = case_when(
      format(Date_parsed, "%m") == "08" ~ "08/2024",
      TRUE ~ DateString
    )
  ) %>% 
  
  # group_by(Location) %>%
  # mutate(SampleLetter = letters[seq_len(n())],
  #        NewDateString = paste0(DateString, "-", SampleLetter)) %>%
  # ungroup() %>%
  dplyr::select(SampleID,Location, DateString, Type, NewDateString)

# Combine drain and sewer samples
replicate_samples <- bind_rows(
  drain_samples,
  sewer_samples %>% dplyr::select(SampleID, Location, DateString, Type, NewDateString)
)

# ---------- Join back to main dataframe ----------
heatmap_df_update <- heatmap_df %>%
  left_join(replicate_samples, by = c("SampleID", "Location", "DateString", "Type")) %>% 
  mutate(
    NewDateString = ifelse(is.na(NewDateString) == TRUE, DateString, NewDateString)
  )

# # ---------- Keep only species detected in >1 sample type per Location ----------
# shared_species_per_location <- heatmap_df_update %>%
#   group_by(Location, name) %>%
#   summarise(n_types = n_distinct(Type), .groups = "drop") %>%
#   filter(n_types > 1)
# 
# # Filter main data for those shared species
# heatmap_df_shared <- heatmap_df_update %>%
#   semi_join(shared_species_per_location, by = c("Location", "name"))

# ---------- Keep only species detected in sewer biofilm and ww sample per Location ----------
shared_species_per_location <- heatmap_df_update %>%
  group_by(Location, name) %>%
  summarise(n_types = n_distinct(Type), .groups = "drop") %>%
  filter(n_types > 1)

# Filter main data for those shared species
heatmap_df_shared <- heatmap_df_update %>%
  semi_join(shared_species_per_location, by = c("Location", "name"))

library(tidytext)

# ---------- STEP 3: Filter for top species ----------
plot_df <- heatmap_df_shared %>%
    mutate(

    Species_italic = paste0("*", name, "*")  # markdown italics
  ) %>% 
  
  mutate(TypeLabel = case_when(TypeLabel == "Endcap" ~ "Sewer\nBiofilm",
                               TypeLabel == "Drain" ~ "Sink\nBiofilm",
                               .default = TypeLabel)) %>% 
  
  mutate(
    TypeLabel = factor(
      TypeLabel,
      levels = c("Sink\nBiofilm", "Sewer\nBiofilm", "Wastewater")
    )
  ) %>% 
  
  mutate(
    FillValue = log10(rel_ab)   # keep continuous relative abundance
  )

#---------top summed relative abundance in sewer biofilm per Site ----------
n_species = 10


top_species_df <- plot_df %>%
  filter(TypeLabel == "Sewer\nBiofilm") %>% 
  
  # complete dataframe so NAs are zero
  # This ensures every SampleID has a row for every Species_italic
  tidyr::complete(SampleID, Species_italic, fill = list(rel_ab = 0)) %>%
  
  group_by(Species_italic) %>%
  summarise(median_abundance = median(rel_ab, na.rm = TRUE), .groups = "drop") %>%
  slice_max(order_by = median_abundance, n = n_species) %>%
  arrange(desc(median_abundance)) %>%
  ungroup()

top_species = top_species_df %>% dplyr::select(Species_italic) %>% distinct() %>% pull()

#reorder species to have top per Site
plot_df_top <- plot_df %>%
  semi_join(top_species_df, by = c( "Species_italic")) %>% 
  # mutate(
  #   Species_italic = reorder_within(
  #     Species_italic,
  #     rel_ab,
  #     Location,
  #     fun = median
  #   )
  # ) %>% 
  left_join(
    top_species_df %>% dplyr::select(Species_italic, median_abundance),
    by = "Species_italic"
  ) %>%

  mutate(
    Species_italic = forcats::fct_reorder(
      Species_italic,
      median_abundance,
      .desc = TRUE
    )
  ) %>%
  
  
  # now replace Site names 
  mutate(Location = case_when(Location == "MARIO" ~ "Site 1",
                              Location == "FIONA" ~ "Site 2",
                              Location == "LUIGI" ~ "Site 3",
                              Location == "SHREK" ~ "Site 4",
                              Location == "OSCAR" ~ "Site 5",
                              .default = Location)) 



CVs_df = plot_df_top %>% 
  
  
  #group by name,Type,Location
  group_by(name,Type,Location) %>% 
  mutate(s_dev = sd(rel_ab,na.rm = TRUE)) %>% 
  mutate(avg = mean(rel_ab,na.rm = TRUE)) %>% 
  mutate(cv_ratio = s_dev/avg) %>% 
  mutate(cv = cv_ratio * 100) %>% 
  ungroup() %>% 
  
  mutate(TypeLabel = case_when(TypeLabel == "Sewer\nBiofilm" ~ "Se",
                               TypeLabel == "Sink\nBiofilm" ~ "Si",
                               TypeLabel == "Wastewater" ~ "WW",
                               .default = TypeLabel)) %>% 
  
  mutate(
    TypeLabel = factor(
      TypeLabel,
      levels = c("WW", "Si", "Se")
    )
  ) %>%
  
  # cv category
  mutate(
    CV_category = case_when(
      cv > 100 ~ "Unstable (>100%)",
      cv > 60  & cv <= 100 ~ "Moderately Stable (60–100%)",
      cv > 30  & cv <= 60  ~ "Stable (30–60%)",
      cv > 15  & cv <= 30  ~ "Very Stable (15–30%)",
      cv <= 15 ~ "Highly Stable (<15%)",
      TRUE ~ NA_character_
    )
  ) %>% 
  
    mutate(
    CV_category = factor(
      CV_category,
      levels = c(
        "Unstable (>100%)",
        "Moderately Stable (60–100%)",
        "Stable (30–60%)",
        "Very Stable (15–30%)",
        "Highly Stable (<15%)"
      )
    )
  ) %>% 
  
  #filter to only ww
  filter(TypeLabel=="WW" | TypeLabel == "Se")%>% 
  
  # now replace Site names 
  mutate(Location = case_when(Location == "MARIO" ~ "Site 1",
                              Location == "FIONA" ~ "Site 2",
                              Location == "LUIGI" ~ "Site 3",
                              Location == "SHREK" ~ "Site 4",
                              Location == "OSCAR" ~ "Site 5",
                              .default = Location))


# ---------- STEP 4: Plot as abundance heatmap ----------
heat_map <- ggplot(CVs_df, aes(x = Type, y = Species_italic, fill = CV_category)) +
  geom_raster()+
  #scale_fill_gradient(low = "white", high = "red", name = "CV") +

  scale_fill_manual(
    values = c(
      "Unstable (>100%)"             = "#c86c69",  # muted coral
      "Moderately Stable (60–100%)"     = "#dba294",  # soft warm rose-peach
      "Stable (30–60%)"              = "#cdd3e5",  # pale blue-gray/lavender neutral
      "Very Stable (15–30%)"              = "#9fdeca",  # light teal
      "Highly Stable (<15%)"        = "#6d7ecd"   # muted blue
    ),
    na.value = "white",   # ← THIS is the key line
    name = "Temporal stability\n(CV of concentration)"
  ) +

  facet_grid(Location ~ TypeLabel, scales = "free", space = "free", switch = "y",
             labeller = labeller(TypeLabel = label_wrap_gen(width = 8))) +
  scale_x_discrete(
    labels = plot_df %>%
      distinct(SampleID, Date, TypeLabel) %>%
      arrange(TypeLabel, Date) %>%
      mutate(Label = paste0(Date)) %>%
      dplyr::select(SampleID, Label) %>%
      tibble::deframe()
  ) +
  scale_y_discrete(limits = rev(levels(CVs_df$Species_italic))) +  # 🔑 reverse
  #scale_y_reordered() +   # 👈 THIS removes "___FIONA" etc.
  labs(x = "CV", y = ""
       #, title = paste0("Top ", n_species, " Species by Summed Relative Abundance Across Dates and Sample Types")
       ) +
  theme(
    #axis.text.x = element_markdown(angle = 30, size = 7, vjust = 0.9, hjust = 0.9),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_markdown(color = "black"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.position = "bottom", 
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA),
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(size = 10, face = "bold", color = "black")
  )

print(heat_map)


# Save final heatmap
ggsave(
  plot = heat_map,
  filename = paste0("Biofilm Project Figures/Kraken2 CT 0.5 CV Heat Map Top ", n_species, " Species by Median Relative Abundance Across Dates in Sewer Biofil.png"),
  width = 2.9,
  height = 8,
  units = "in"
)

#install.packages("svglite")
library("svglite")

# Save final heatmap
ggsave(
  plot = heat_map,
  filename = "Biofilm Project Figures/Kraken CT 0.5 CV Heat Map -RPIP only Kraken Top 10 Species by Median Relative Abundance Across Dates in Sewer Biofilm.svg",
  width = 2.9,
  height = 8
)

plot(heat_map)

#show relative abundance box plot
plot_df_box <- plot_df_top %>%
  mutate(
    log_rel_ab = log10(rel_ab + 1e-6)  # avoid -Inf
  )%>% 
  
  # now replace Site names 
  mutate(Location = case_when(Location == "MARIO" ~ "Site 1",
                              Location == "FIONA" ~ "Site 2",
                              Location == "LUIGI" ~ "Site 3",
                              Location == "SHREK" ~ "Site 4",
                              Location == "OSCAR" ~ "Site 5",
                              .default = Location))

species_levels <- top_species_df %>%
  arrange(desc(median_abundance)) %>%
  pull(Species_italic)

y_scaffold <- expand.grid(
  Species_italic = species_levels,   # already globally ordered
  Location = unique(plot_df_box$Location)
) %>%
  mutate(
    Species_italic = factor(Species_italic, levels = species_levels)
  )

plot_df_box_final <- plot_df_box %>%
  mutate(
    Species_italic = factor(Species_italic, levels = species_levels)
  ) %>% 
  mutate(
    TypeLabel = factor(
      TypeLabel,
      levels = c("Wastewater", "Sewer\nBiofilm", "Sink\nBiofilm")
    )
  ) 


# box_plot <- ggplot(
#   plot_df_box_final,
#   aes(x = log_rel_ab, y = Species_italic)
# ) +
#   
# 
#   
#   geom_boxplot(
#     outlier.size = 0.4,
#     width = 0.6,
#     color = "black"
#   ) +
#   
#   
#   
#   facet_grid(
#     Location ~ TypeLabel,
#     scales = "free_y",
#     space = "free",
#     switch = "y"
#   ) +
#   scale_y_reordered() +
#   labs(
#     x = "Log10(Relative Abundance)",
#     y = ""
#   ) +
#   theme(
#     panel.grid = element_blank(),
#     axis.text.y = element_markdown(size = 8),
#     axis.text.x = element_text(size = 8),
#     panel.border = element_rect(color = "black", fill = NA),
#     strip.background = element_rect(fill = "white", color = "white"),
#     strip.text = element_text(size = 10, face = "bold"),
#     legend.position = "none"
#   )

box_plot <- ggplot(
  plot_df_box_final,
  aes(x = log_rel_ab, y = Species_italic)
) +
  
  # 👻 invisible scaffold (THIS is the magic)
  geom_blank(
  data = y_scaffold,
  aes(y = Species_italic),
  inherit.aes = FALSE
  ) +
  
  # Add dots for individual samples
  geom_jitter(
    aes(color = TypeLabel),
    height = 0.15,
    size = 1,
    alpha = 0.8
  ) +
  # Boxplot outline only
  geom_boxplot(
    outlier.shape = NA,  # hide default outliers because we’ll add points
    width = 0.6,
    color = "black",
    fill = NA
  ) +
  
  facet_grid(
  Location ~ TypeLabel,
  scales = "free_y",
  space = "free",
  switch = "y",
  drop = FALSE   # 🔑 THIS IS THE MAGIC
)+
  scale_y_discrete(limits = rev(levels(plot_df_box_final$Species_italic))) +  # 🔑 reverse
  #scale_y_reordered() +
  labs(
    x = "Log10(Relative Abundance)",
    y = "",
    color = "Sample Type"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_markdown(size = 8),
    axis.text.x = element_text(size = 8),
    panel.background = element_rect(fill = "white", color = NA),  # ← fix
    panel.border = element_rect(color = "black", fill = NA),
    plot.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "bottom"
  ) +
  scale_color_manual(values = c( 
    "Sink\nBiofilm" =  "#9fdeca",
    "Sewer\nBiofilm" = "#c86c69",
    "Wastewater" = "#6d7ecd"
  ))


plot(box_plot)


# Save final heatmap
ggsave(
  plot = box_plot,
  filename = paste0("Biofilm Project Figures/Kraken2 CT 0.5 Box Plot- Top ", n_species, " Species by Median Relative Abundance Across Dates in Sewer.png"),
  width = 5,
  height = 8,
  units = "in"
)

#install.packages("svglite")
library("svglite")

# Save final heatmap
ggsave(
  plot = box_plot,
  filename = "Biofilm Project Figures/Kraken CT 0.5 Box Plot -RPIP only Kraken Top 10 Species by Median Relative Abundance Across Dates in Sewer.svg",
  width = 5,
  height = 8
)


# Create a reference table for CVs per pathogen × Type × Location
CV_table <- CVs_df %>%
  group_by(name, Type, Location) %>%
  summarise(
    mean_rel_ab = mean(rel_ab, na.rm = TRUE),
    sd_rel_ab = sd(rel_ab, na.rm = TRUE),
    cv_percent = mean(cv, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(Location, Type, desc(cv_percent))

# View the first rows
head(CV_table)


```

#UpSet Plots
##Faceted Genus Level Upset-ONLY DATES SHARED WITH BIOFILM
```{r}
#specify dates
dates = c("2024-09-09", "2024-09-11", "2024-09-12", "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09")

# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame


# RPIP taxa
rpip_targets <- rpip_pathogen_species_list


# ---------- Extract and clean ----------
analysis_df <- kraken_merged_df_filtered %>%
  filter(
    taxonomy_lvl == "S",
    Location %in% locations,
    name %in% rpip_targets,
    Type %in% c("Drain", "Endcap", "H_WW")
  ) %>%
  mutate(
    rel_ab = as.numeric(fraction_total_reads),
    Date   = as.Date(Date),
    Type   = factor(Type, levels = c("Drain", "Endcap", "H_WW"))
  ) %>%
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>%
  drop_na(timepoint)


biofilm_dates <- analysis_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
   mutate(
    biofilm_ref_date = Date
  ) %>% 
    ungroup() %>% 
  dplyr::distinct(Location,timepoint,biofilm_ref_date)

filtered_df <- analysis_df %>%
  mutate(Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  ) %>% 
  
  dplyr::select(-best_date,-date_diff_days,-biofilm_ref_date)

df <- filtered_df %>%
  
  
  mutate(
    Type = factor(Type, levels = c("Drain", "Endcap", "H_WW")),
    TypeLabel = case_when(
      Type == "H_WW" ~ "Wastewater", 
      TRUE ~ as.character(Type)
    ),
    Month = format(as.Date(Date), "%m-%Y")
  )


# ---------- Filter for present species ----------
plot_df <- df %>%
  #inner_join(top_species_per_location %>% dplyr::select(Location, name),
             #by = c("Location", "name")) %>%
  # mutate(
  #   # Order species by abundance within each location for cleaner y-axis
  #   name = reorder_within(name, rel_ab, Location),
  #   FillValue = rel_ab
  # ) %>% 
  dplyr::select(name,SampleID, Date, Location, Type,rel_ab) %>% 
  
  
  #make presence/absence
  mutate(rel_ab = ifelse(rel_ab>0,1,0)) %>%
  mutate(
    Type  = dplyr::recode(Type,
                   "H_WW"   = "Wastewater",
                   "Endcap" = "Sewer Biofilm",
                   "Drain"  = "Sink Biofilm")
  )%>% 
  
  # now replace Site names 
  mutate(Location = case_when(Location == "MARIO" ~ "Site 1",
                              Location == "FIONA" ~ "Site 2",
                              Location == "LUIGI" ~ "Site 3",
                              Location == "SHREK" ~ "Site 4",
                              Location == "OSCAR" ~ "Site 5"))
  
wide_df <- plot_df %>%
  group_by(Location, name, Type) %>%
  summarise(detect = as.integer(sum(rel_ab) > 0), .groups = "drop") %>%
  pivot_wider(
    names_from = Type,
    values_from = detect,
    values_fill = 0
  )



# install.packages("UpSetR")
# library(UpSetR)



library(ggupset)
library(forcats)

n_class = 45


wide_long <- wide_df %>%
  pivot_longer(cols = c("Sink Biofilm","Sewer Biofilm","Wastewater"),
               names_to = "Type",
               values_to = "detect") %>%
  filter(detect == 1) %>%
  group_by(Location, name) %>%
  summarise(sets = list(unique(Type)), .groups = "drop") %>%
  ungroup() %>% 
  mutate(class = sub(" .*", "", name)) %>%
  # keep top N globally
  mutate(class = fct_lump_n(class, n = n_class, other_level = "Other"))

class_freq <- wide_long %>%
  mutate(class = as.character(class)) %>%   # 🔑 FIX
  dplyr::count(class, sort = TRUE) %>% 
  mutate(rank = dense_rank(dplyr::desc(n))) %>%
  arrange(rank)



keep_classes <- class_freq %>%
  slice_head(n = n_class) %>%
  pull(class)

wide_long_cut <- wide_long %>%
  mutate(class = ifelse(as.character(class) %in% keep_classes,
                        as.character(class),
                        "Other")) %>%
  mutate(class = factor(class))   # re-factor afterSite




# Get levels (after lumping)
class_levels <- levels(wide_long_cut$class)

# Assign colors: first 15 get a palette, "Other" gets grey
palette <- c(fortyfive_pal, "grey70")
names(palette) <- c(setdiff(class_levels, "Other"), "Other")



# build italicized labels
genus_labels <- setNames(
  ifelse(class_levels == "Other",
         "Other",
         paste0("*", class_levels, "*")),   # markdown italics
  class_levels
)

# Plot

plot = ggplot(wide_long, aes(x = sets, fill = class)) +
  geom_bar() +
  scale_x_upset() +
  facet_wrap(~Location, ncol = 5) +
  labs(
    fill = "Genus",
    x = "",
    y = "Number of Pathogen Species\n\ Shared Between Groups"
  ) +
  scale_fill_manual(values = palette, labels = genus_labels) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_rect(colour = "black"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_markdown(size = 7),  # <- enable markdown
    axis.title.y = element_text(margin = margin(r = 5),size = 7) # tighten y label
  ) +
  theme_combmatrix(
    combmatrix.panel.point.color.fill = "black",
    combmatrix.label.make_space = FALSE,   # <- stop reserving extra space
    combmatrix.label.text = element_text(color = "black", size = 6),
    #combmatrix.label.text = element_blank(),
    combmatrix.panel.line.color = "black",
    combmatrix.panel.striped_background = FALSE,
    combmatrix.panel.point.size = 2
  )

print(plot)

ggsave(path = "Biofilm Project Figures",file="Kraken2 CT 0.5 Bracken UpSet plot of pathogens detected by Genus-ONLY DATES SHARED WITH BIOFILM.png", plot=plot, width=7.5, height=5,dpi = 500)
ggsave(path = "Biofilm Project Figures",file="Kraken2 CT 0.5 Bracken UpSet plot of pathogens detected by Genus-ONLY DATES SHARED WITH BIOFILM.svg", plot=plot, width=7.5, height=5)



```

##Faceted Species-level Upset-ONLY SHARED DATES WITH BIOFILM
```{r}

#specify dates
dates = c("2024-09-09", "2024-09-11", "2024-09-12", "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09")


# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame


# RPIP taxa
rpip_targets <-  rpip_pathogen_species_list

# ---------- STEP 1: Extract and clean ----------
analysis_df <- kraken_merged_df_filtered %>%
  filter(
    taxonomy_lvl == "S",
    Location %in% locations,
    name %in% rpip_targets,
    Type %in% c("Drain", "Endcap", "H_WW")
  ) %>%
  mutate(
    rel_ab = as.numeric(fraction_total_reads),
    Date   = as.Date(Date),
    Type   = factor(Type, levels = c("Drain", "Endcap", "H_WW"))
  ) %>%
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>%
  drop_na(timepoint)


biofilm_dates <- analysis_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
  mutate(
    biofilm_ref_date = Date
  ) %>% 
    ungroup() %>% 
  dplyr::distinct(Location,timepoint,biofilm_ref_date)

filtered_df <- analysis_df %>%
  mutate(Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  ) %>% 
  
  dplyr::select(-best_date,-date_diff_days,-biofilm_ref_date)

df <- filtered_df %>%
  
  
  mutate(
    Type = factor(Type, levels = c("Drain", "Endcap", "H_WW")),
    TypeLabel = case_when(
      Type == "H_WW" ~ "Wastewater", 
      TRUE ~ as.character(Type)
    ),
    Month = format(as.Date(Date), "%m-%Y")
  )%>% 
  
  # now replace Site names 
  mutate(Location = case_when(Location == "MARIO" ~ "Site 1",
                              Location == "FIONA" ~ "Site 2",
                              Location == "LUIGI" ~ "Site 3",
                              Location == "SHREK" ~ "Site 4",
                              Location == "OSCAR" ~ "Site 5",
                              .default = Location))


# ---------- STEP: Rank species by sharedness across categories ----------

#  Count in how many unique sample Types each species appears
species_sharedness <- df %>%
  group_by(name) %>%
  summarise(
    n_types_detected = n_distinct(Type[rel_ab > 0]),
    mean_ab = mean(rel_ab, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_types_detected), desc(mean_ab))

n_species_shared  = 45

# Optional: extract top N
top_species <- species_sharedness %>%
  slice_head(n = n_species_shared) %>%
  pull(name)

# ---------- Filter for those species ----------
plot_df <- df %>%
  #inner_join(top_species_per_location %>% dplyr::select(Location, name),
             #by = c("Location", "name")) %>%
  # mutate(
  #   # Order species by abundance within each location for cleaner y-axis
  #   name = reorder_within(name, rel_ab, Location),
  #   FillValue = rel_ab
  # ) %>% 
  dplyr::select(name,SampleID, Date, Location, Type,rel_ab) %>% 
  
  
  #make presence/absence
  mutate(rel_ab = ifelse(rel_ab>0,1,0)) %>%
  mutate(
    Type  = dplyr::recode(Type,
                   "H_WW"   = "Wastewater",
                   "Endcap" = "Sewer Biofilm",
                   "Drain"  = "Sink Biofilm")
  )
  
  
wide_df <- plot_df %>%
  group_by(Location, name, Type) %>%
  summarise(detect = as.integer(sum(rel_ab) > 0), .groups = "drop") %>% 
  pivot_wider(
    names_from = Type,
    values_from = detect,
    values_fill = 0
  )



# install.packages("UpSetR")
# library(UpSetR)

library(ggupset)
library(forcats)


wide_long <- wide_df %>%
  pivot_longer(cols = c("Sink Biofilm","Sewer Biofilm","Wastewater"),
               names_to = "Type",
               values_to = "detect") %>%
  filter(detect == 1) %>%
  group_by(Location,name) %>%
  summarise(
    sets = list(unique(Type)),
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(
    name = ifelse(name %in% top_species, name, "Other"),
    name = factor(name, levels = c(sort(top_species), "Other"))
  )





# Get levels (after lumping)
species_levels <- levels(wide_long$name)

library(RColorBrewer)

# Assign colors: first N get a palette, "Other" gets grey
palette <- c(fortyfive_pal[1:n_species_shared], "grey70")
names(palette) <- c(setdiff(species_levels, "Other"), "Other")

# library(viridis)
# 
# cols <- viridis(12, option = "A")  # or "C", "B", "A", "E", "F"
# palette <- c(cols, "grey70")
# names(palette) <- c(setdiff(species_levels, "Other"), "Other")


# Plot

plot = ggplot(wide_long, aes(x = sets, fill = name)) +
  geom_bar() +
  scale_x_upset() +
  facet_wrap(~Location, ncol = 5) +
  labs(fill = "",
       x = "",
       y = "Number of Pathogen Species\n\ Shared Between Groups") +
  scale_fill_manual(values = palette) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    #strip.background = element_rect(colour = "black", fill = "white"),
    strip.background = element_rect(colour = "black"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),      # remove labels
    axis.ticks.x = element_blank(),     # remove ticks
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 5)
  )+
  theme_combmatrix(combmatrix.panel.point.color.fill = "black",
                   combmatrix.label.make_space = TRUE,
                   #combmatrix.label.text = element_blank(),  # << hide combmatrix labels
                   combmatrix.label.text = element_text(color = "black",size =7),
                   combmatrix.panel.line.color = "black",
                   combmatrix.panel.striped_background = FALSE,
                   combmatrix.panel.point.size = 2
                     #combmatrix.panel.line.size = 0,
                   ) 

print(plot)

ggsave(path = "Biofilm Project Figures",file="Kraken2 CT 0.5 Bracken UpSet plot of 45 pathogens detected by species (top shared across cats)-ONLY SHARED DATES WITH BIOFILM.png", plot=plot, width=7.5, height=5,dpi = 500)
ggsave(path = "Biofilm Project Figures",file="Kraken2 CT 0.5 Bracken UpSet plot of pathogens detected by species (top shared across cats)-ONLY SHARED DATES WITH BIOFILM.svg", plot=plot, width=7.5, height=5)

```



#UpSet by Pathogen Plots Relative Abundance and Presence Absence
##upset plot set themes
```{r}
font_add_google("Reem Kufi", "Reem Kufi")   # <-- REQUIRED
showtext_auto()                             # <-- activates showtext


## ggplot theme
theme_set(theme_minimal(base_family = "Reem Kufi", base_size = 12))

theme_update(
  plot.title = element_text(size = 27,
                            face = "bold",
                            hjust = .5,
                            margin = margin(10, 0, 30, 0)),
  plot.caption = element_text(size = 9,
                              color = "grey40",
                              hjust = .5,
                              margin = margin(20, 0, 5, 0)),
  axis.text.y = element_blank(),
  axis.title = element_blank(),
  plot.background = element_rect(fill = "grey88", color = NA),
  panel.background = element_rect(fill = NA, color = NA),
  panel.grid = element_blank(),
  panel.spacing.y = unit(0, "lines"),
  strip.text.y = element_blank(),
  legend.position = "bottom",
  legend.text = element_text(size = 9, color = "grey40"),
  legend.box.margin = margin(t = 30), 
  legend.background = element_rect(color = "grey40", 
                                   size = .3, 
                                   fill = "grey95"),
  legend.key.height = unit(.25, "lines"),
  legend.key.width = unit(2.5, "lines"),
  plot.margin = margin(rep(20, 4))
)

#second-------------------

font_add_google("Chivo", "Chivo")
font_add_google("Passion One", "Passion One")

## ggplot theme
theme_set(theme_minimal(base_family = "Chivo"))

theme_update(axis.text.x = element_text(size = 11, color = "grey20"),
             axis.text.y = element_text(size = 13, color = "black"),
             axis.ticks.x = element_line(color = "grey45"),
             axis.ticks.y = element_blank(),
             axis.ticks.length.x = unit(.4, "lines"),
             panel.grid = element_blank(),
             plot.background = element_rect(fill = "grey60", color = "grey60"))
```

##Upset Plot - All Pathogens Presence Absence- All Sites in One-RPIP only
```{r eval=FALSE, include=FALSE}
#theme set----------

## ggplot theme
font_add_google("Reem Kufi", "Reem Kufi")   # <-- REQUIRED
showtext_auto()                             # <-- activates showtext

theme_set(theme_minimal(base_family = "Reem Kufi", base_size = 12))

theme_update(
  plot.title = element_text(size = 27,
                            face = "bold",
                            hjust = .5,
                            margin = margin(10, 0, 30, 0)),
  plot.caption = element_text(size = 9,
                              color = "grey40",
                              hjust = .5,
                              margin = margin(20, 0, 5, 0)),
  axis.text.y = element_blank(),
  axis.title = element_blank(),
  plot.background = element_rect(fill = "grey88", color = NA),
  panel.background = element_rect(fill = NA, color = NA),
  panel.grid = element_blank(),
  panel.spacing.y = unit(0, "lines"),
  strip.text.y = element_blank(),
  legend.position = "bottom",
  legend.text = element_text(size = 9, color = "white"),
  legend.box.margin = margin(t = 30), 
  legend.background = element_rect(color = "white", 
                                   size = .3, 
                                   fill = "grey95"),
  legend.key.height = unit(.25, "lines"),
  legend.key.width = unit(2.5, "lines"),
  plot.margin = margin(rep(20, 4))
)
#--------
font_add_google("Chivo", "Chivo")
font_add_google("Passion One", "Passion One")

## ggplot theme
theme_set(theme_minimal(base_family = "Chivo"))

theme_update(
              #axis.text.x = element_text(size = 11, color = "grey20"),
             #axis.text.y = element_text(size = 13, color = "black", face = "bold"),
             axis.ticks.x = element_line(color = "white"),
             axis.ticks.y = element_blank(),
             axis.ticks.length.x = unit(.4, "lines"),
             panel.grid = element_blank(),
             plot.background = element_rect(fill = "white", color = "white"))



#R prep---------------------------------------
## packages
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggtext)
library(showtext)
library(sysfonts)
library(ggplot2)
library(ggstatsplot)
library(dplyr)
library(tidyverse)
library(readxl)
library(vegan)
library(OTUtable)
library(vegan)
library(MASS)
library(ggtext)
library(ggupset)
library(forcats)

font_add_google("Chivo", "Chivo")
font_add_google("Passion One", "Passion One")

## ggplot theme
theme_set(theme_minimal(base_family = "Chivo"))

theme_update(
              #axis.text.x = element_text(size = 11, color = "grey20"),
             #axis.text.y = element_text(size = 13, color = "black", face = "bold"),
             axis.ticks.x = element_line(color = "white"),
             axis.ticks.y = element_blank(),
             axis.ticks.length.x = unit(.4, "lines"),
             panel.grid = element_blank(),
             plot.background = element_rect(fill = "white", color = "white"))

#R data -------------------------
# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame



# RPIP taxa
rpip_targets <-  rpip_pathogen_species_list

df <- kraken_merged_df_filtered %>%
  filter(taxonomy_lvl == "S") %>% 
  mutate(rel_ab = as.numeric(fraction_total_reads)) %>% 
  filter(Location %in% locations) %>% 
  #filter(name %in% c(rpip_targets,detected_pathogens)) %>% 
  filter(name %in% c(rpip_targets)) %>%

  filter(Type %in% c("Drain", "Endcap", "H_WW")) %>% 
  mutate(
    Type = factor(Type, levels = c("Drain", "Endcap", "H_WW")),
    TypeLabel = case_when(Type == "H_WW" ~ "Hospital Wastewater", 
                          Type == "Drain" ~ "Sink Biofilm", 
                          Type == "Endcap" ~ "Sewer Biofilm", 
                          TRUE ~ as.character(Type)),
    Month = format(as.Date(Date), "%m-%Y")
  ) 


# Count detects per Species × Type
df_detects <- df %>%
  group_by(name, SampleID) %>% 
  mutate(sum_rel_ab = sum(rel_ab,na.rm = TRUE)) %>% 
  ungroup() %>% 
  
  dplyr::select(name,Type,Location,SampleID,sum_rel_ab) %>% distinct() %>%
  
  mutate(detect = ifelse(sum_rel_ab>0,1,0)) %>% 
  
  group_by(name, Type) %>%
  summarise(
    n_detects_type = sum(detect, na.rm = TRUE),   # how many samples detected?
    .groups = "drop"
  ) %>%
  ungroup() %>% 
  
  filter(n_detects_type >0)

# Total detects per species (for ordering)
detect_order <- df_detects %>%
  group_by(name) %>%
  summarise(n_detects = sum(n_detects_type)
         #, .groups = "drop"
         ) %>%
  arrange(desc(n_detects))

# compute detection + average rel_ab
df_ranks_type <- df_detects %>%
  left_join(detect_order, by = "name") %>%
  
mutate(ID = rank(-n_detects_type, ties.method = "first")) %>%
  
  mutate(
    name = factor(
      name,
      levels = unique(name[order(n_detects, -n_detects_type)])
    )
  ) %>% 
  
  mutate(
    TypeLabel = case_when(
      Type == "H_WW" ~ "Hospital Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain" ~ "Sink Biofilm",
      .default = Type
    ),
    TypeLabel = factor(TypeLabel, 
                       levels = c("Hospital Wastewater", 
                                  "Sewer Biofilm", 
                                  "Sink Biofilm")),
    type_num = as.numeric(TypeLabel)
  ) %>% 
  
 mutate(Species_label_italic = paste0("<i>", name, "</i>")) %>%  # Convert Species_label to italicized labels

  
  mutate(
    Species_label_italic = fct_reorder(Species_label_italic, n_detects, .desc = FALSE)
  ) %>% 
  
  filter(n_detects >0)
  




#Plot-------------
max_detect <- max(df_ranks_type$n_detects)

# candidate step
raw_step <- max_detect / 4

# function to snap to nice values (2, 5, 10, etc.)
nice_step <- function(x) {
  base <- 10^floor(log10(x))
  if (x/base <= 2) return(2 * base)
  if (x/base <= 5) return(5 * base)
  return(10 * base)
}

break_step <- 50

bars <- 
  df_ranks_type %>% 
  ggplot(aes(Species_label_italic, -n_detects_type)) +
    geom_col(aes(fill = TypeLabel),
             #color = "white",
             size = .5,
             width = 1.02) +
    # geom_curve(aes(x = 51.2, xend = 47, 
    #                y = -148, yend = -166),
    #            curvature = -.4) +
    # annotate("text", x = 47, y = -300, 
    #          label = "Each rectangle represents\none song included in the\nBBC ranking, its length\n the total points and the\ncolor indicates the rank",
    #          family = "Chivo",
    #          size = 2,
    #          lineheight = .9) +
    # annotate("text", x = 21.5, y = -120, 
    #          label = 'The Top Artists featured in the BBC´s\n"Greatest Hip-Hop Songs of All Time"',
    #          family = "Chivo",
    #          fontface = "bold",
    #          size = 4,
    #          lineheight = .9) +
    # annotate("text", x = 17, y = -120,
    #          label = 'In Autumn 2019, 108 hip-hop and music experts ranked their 5 favorites out of\n311 nominated songs in an online survey by the BBC. The graphic shows points\nscored in total and per song for the top ranked artists and broken down by era.',
    #          family = "Chivo",
    #          fontface = "bold",
    #          color = "grey30",
    #          size = 2.5,
    #          lineheight = .9) +
    coord_flip(clip = "off") +
    scale_x_discrete(position = "top") +
    scale_y_continuous(
      expand = c(.02, .02),
      limits = c(
        -max_detect - max_detect/4,
        0
      ),
      breaks = seq(
        -ceiling(max_detect / break_step) * break_step,
        0,
        by = break_step
      ),
      labels = function(x) abs(as.integer(x)),
      position = "right"
    )+
    # nord::scale_fill_nord(palette = "halifax_harbor",
    #                       discrete = T,
    #                       reverse = F
    #                       #guide = F
    #                       ) +
    scale_fill_manual(values = c("#6d7ecd", "#9fdeca", "#c86c69"))+
    theme(
          axis.text.y.right = element_markdown(hjust = .5,size = 7),
          axis.text.x = element_text(size = 12),
          axis.ticks.x = element_line(color = "black"),
          axis.ticks.length = unit(0.5, "pt"),
          #axis.text.y = element_blank(),
          #axis.text.x = element_markdown(size = 35),
          axis.title = element_text(size = 14,face = "bold"),
          legend.position = "left",
          legend.text = element_text(size = 12),
          plot.margin = margin(5, 0, 5, 5)) +
    labs(x = NULL, 
         fill = "",
         y = "Number of Detections",
         title = "")


dots <- df_ranks_type %>% 

  ggplot(aes(Species_label_italic, type_num, group = Species_label_italic)) +
    geom_point(aes(Species_label_italic, 1), color = "grey75", size = 2) +
    geom_point(aes(Species_label_italic, 2), color = "grey75", size = 2) +
    geom_point(aes(Species_label_italic, 3), color = "grey75", size = 2) +
  
    geom_segment(aes(x = Species_label_italic, xend = Species_label_italic, 
                     y = 1, yend = 3), 
                 color = "grey75",
                 size = .3) +
    geom_line(color = "black",
              size = .09) +
    geom_point(aes(fill = TypeLabel, size = n_detects_type),
               shape = 21,
               color = "black",
               stroke = 1.2) +
  # geom_point(aes(fill = TypeLabel, size = n_detects_type), 
  #          shape = 21, 
  #          color = "black", 
  #          stroke = 1) +

    #geom_curve(aes(x = 47, xend = 51, 
    # geom_curve(aes(x = 45, xend = 49, 
    #                y = 6.1, yend = 4.3),
    #            curvature = .4) +
    #annotate("text", x = 45.1, y = 6.1, 
    # annotate("text", x = 43.1, y = 6.1, 
    #          label = "The dot size indicates\nthe number of songs,\nthe dot color the best\nrank in each era",
    #          family = "Chivo",
    #          size = 3.8, 
    #          lineheight = .9) +
    coord_flip() +
    scale_y_continuous(limits = c(.5, 7.3),
                       breaks = 1:3,
                       labels = c("Hospital Wastewater", 
                                  "Sewer Biofilm", 
                                  "Sink Biofilm"), 
                       position = "right") +
  scale_size(range = c(0.6, 2.7), guide = FALSE)+
  #scale_size(range = c(2, 5.5), guide = F) +
  scale_fill_manual(values = c("#6d7ecd", "#9fdeca", "#c86c69"))+
    # nord::scale_fill_nord(palette = "halifax_harbor", 
    #                       discrete = T, 
    #                       reverse = F, 
    #                       guide = F) +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 7, hjust = .01, vjust = 0.1, angle = 35),
          #axis.text.x = element_text(size = 35, face = "bold", angle = 35),
          #axis.text.y = element_blank(),
          #axis.text.y = element_markdown(size = 20,face = "bold"),
          #axis.text.y = element_text(size = 20,face = "bold"),
          axis.text.y = element_blank(),
          legend.position = "none",
          plot.margin = margin(5, 5, 5, 0),
          plot.caption = element_text(face = "bold", color = "grey30", 
                                      size = 15, margin = margin(t = 15))) +
    labs(x = NULL, y = NULL,
         caption = "")


plot = bars + dots +  plot_layout(widths = c(1, .35))

showtext_opts(dpi = 300)

ggsave(path = "Biofilm Project Figures",file="Kraken2CT0.5-Upset Plot Detections for Pathogens.png", plot=plot, width=8, height= 9,dpi = 300)
ggsave(path = "Biofilm Project Figures",file="Kraken2CT0.5-Upset Plot Detections for Pathogens.svg", plot=plot, width=8, height=9)

print(plot)



#------------ Get counts of pathogens for each environment --------

presence_df <- df %>%
  mutate(
    TypeLabel = case_when(
      Type == "H_WW"   ~ "Hospital Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain"  ~ "Sink Biofilm"
    )
  ) %>%
  group_by(name, TypeLabel) %>%
  summarise(present = any(rel_ab > 0), .groups = "drop") %>%
  filter(present)


presence_wide <- presence_df %>%
  mutate(present = 1) %>%
  tidyr::pivot_wider(
    names_from  = TypeLabel,
    values_from = present,
    values_fill = 0
  )


#only ww
ww_only <- presence_wide %>%
  filter(
    `Hospital Wastewater` == 1,
    `Sewer Biofilm` == 0,
    `Sink Biofilm` == 0
  ) %>%
  pull(name)


# sewer biofilm and wastewater
sewer_ww <- presence_wide %>%
  filter(
    `Hospital Wastewater` == 1,
    `Sewer Biofilm` == 1,
    `Sink Biofilm` == 0 |`Sink Biofilm` == 1
  ) %>%
  pull(name)


# all three

all_three <- presence_wide %>%
  filter(
    `Hospital Wastewater` == 1,
    `Sewer Biofilm` == 1,
    `Sink Biofilm` == 1
  ) %>%
  pull(name)

#sink biofilm and ww
sink_ww <- presence_wide %>%
  filter(
    `Hospital Wastewater` == 1,
    `Sewer Biofilm` == 0,
    `Sink Biofilm` == 1
  ) %>%
  pull(name)

#------------ Get summed counts of pathogens for each environment --------

presence_filtered_df <- df %>%
  mutate(
    TypeLabel = case_when(
      Type == "H_WW"   ~ "Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain"  ~ "Sink Biofilm"
    )
  ) %>%
  dplyr::select(SampleID,name,rel_ab,TypeLabel) %>% distinct() %>% 
  filter(rel_ab>0) %>% 
  group_by(name, TypeLabel) %>%
  summarise(present_counts = n(), .groups = "drop") %>% 
  ungroup()


presence_wide <- presence_filtered_df %>%
  tidyr::pivot_wider(
    names_from  = TypeLabel,
    values_from = present_counts,
    values_fill = 0
  )

presence_wide_filtered = presence_wide %>% 
  filter(`Wastewater`>1 & `Sewer Biofilm`>1 & `Sink Biofilm`>1)
```

##Upset Plot - All Pathogens Presence Absence- All Sites in One-0NLY SHARED DATES WITH BIOFILM-RPIP only
```{r eval=FALSE, include=FALSE}
#theme set----------

## ggplot theme
font_add_google("Reem Kufi", "Reem Kufi")   # <-- REQUIRED
showtext_auto()                             # <-- activates showtext

theme_set(theme_minimal(base_family = "Reem Kufi", base_size = 12))

theme_update(
  plot.title = element_text(size = 27,
                            face = "bold",
                            hjust = .5,
                            margin = margin(10, 0, 30, 0)),
  plot.caption = element_text(size = 9,
                              color = "grey40",
                              hjust = .5,
                              margin = margin(20, 0, 5, 0)),
  axis.text.y = element_blank(),
  axis.title = element_blank(),
  plot.background = element_rect(fill = "grey88", color = NA),
  panel.background = element_rect(fill = NA, color = NA),
  panel.grid = element_blank(),
  panel.spacing.y = unit(0, "lines"),
  strip.text.y = element_blank(),
  legend.position = "bottom",
  legend.text = element_text(size = 9, color = "grey40"),
  legend.box.margin = margin(t = 30), 
  legend.background = element_rect(color = "grey40", 
                                   size = .3, 
                                   fill = "grey95"),
  legend.key.height = unit(.25, "lines"),
  legend.key.width = unit(2.5, "lines"),
  plot.margin = margin(rep(20, 4))
)
#--------
font_add_google("Chivo", "Chivo")
font_add_google("Passion One", "Passion One")

## ggplot theme
theme_set(theme_minimal(base_family = "Chivo"))

theme_update(axis.text.x = element_text(color = "grey12"),
             #axis.text.y = element_text(color = "black", face = "bold"),
             axis.ticks.x = element_line(color = "white"),
             axis.ticks.y = element_blank(),
             axis.ticks.length.x = unit(.4, "lines"),
             panel.grid = element_blank(),
             plot.background = element_rect(fill = "white", color = "white"))



#R prep---------------------------------------
## packages
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggtext)
library(showtext)
library(sysfonts)
library(ggplot2)
library(ggstatsplot)
library(dplyr)
library(tidyverse)
library(readxl)
library(vegan)
library(OTUtable)
library(vegan)
library(MASS)
library(ggtext)
library(ggupset)
library(forcats)

font_add_google("Chivo", "Chivo")
font_add_google("Passion One", "Passion One")

## ggplot theme
theme_set(theme_minimal(base_family = "Chivo"))

theme_update(axis.text.x = element_text(color = "grey20"),
             axis.text.y = element_text(color = "black", face = "bold"),
             axis.ticks.x = element_line(color = "white"),
             axis.ticks.y = element_blank(),
             axis.ticks.length.x = unit(.4, "lines"),
             panel.grid = element_blank(),
             plot.background = element_rect(fill = "white", color = "white"))

#R data -------------------------
dates = c("2024-09-09", "2024-09-11", "2024-09-12", "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09")

# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame



analysis_df <- kraken_merged_df_filtered %>%
  filter(
    taxonomy_lvl == "S",
    Location %in% locations,
    name %in% rpip_pathogen_species_list,
    #name %in% c(rpip_pathogen_species_list,detected_pathogens),
    Type %in% c("Drain", "Endcap", "H_WW")
  ) %>% 
  mutate(
    rel_ab = as.numeric(fraction_total_reads),
    Date   = as.Date(Date),
    Type   = factor(Type, levels = c("Drain", "Endcap", "H_WW"))
  ) %>%
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>% 
  drop_na(timepoint)


biofilm_dates <- analysis_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
  mutate(
      biofilm_ref_date = Date
    ) %>% 
  ungroup() %>% 
  dplyr::distinct(Location,timepoint,biofilm_ref_date)

filtered_df <- analysis_df %>%
  mutate(Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  ) %>% 
  
  dplyr::select(-best_date,-date_diff_days,-biofilm_ref_date)

df <- filtered_df %>%
  
  mutate(
    Type = factor(Type, levels = c("Drain", "Endcap", "H_WW")),
    TypeLabel = case_when(Type == "H_WW" ~ "Wastewater", 
                          Type == "Drain" ~ "Sink Biofilm", 
                          Type == "Endcap" ~ "Sewer Biofilm", 
                          TRUE ~ as.character(Type)),
    Month = format(as.Date(Date), "%m-%Y")
  )%>% 
  
  # now replace Site names 
  mutate(Location = case_when(Location == "MARIO" ~ "Site 1",
                              Location == "FIONA" ~ "Site 2",
                              Location == "LUIGI" ~ "Site 3",
                              Location == "SHREK" ~ "Site 4",
                              Location == "OSCAR" ~ "Site 5",
                              .default = Location))


# Count detects per gene × Type
df_detects <- df %>%
  group_by(name, SampleID) %>% 
  mutate(sum_rel_ab = sum(rel_ab,na.rm = TRUE)) %>% 
  ungroup() %>% 
  
  dplyr::select(name,Type,Location,SampleID,sum_rel_ab) %>% distinct() %>%
  
  mutate(detect = ifelse(sum_rel_ab>0,1,0)) %>% 
  
  group_by(name, Type) %>%
  summarise(
    n_detects_type = sum(detect, na.rm = TRUE),   # how many samples detected?
    .groups = "drop"
  ) %>%
  ungroup() %>% 
  
  filter(n_detects_type >0)

# Total detects per species (for ordering)
detect_order <- df_detects %>%
  group_by(name) %>%
  summarise(n_detects = sum(n_detects_type)
         #, .groups = "drop"
         ) %>%
  arrange(desc(n_detects))

# compute detection + average rel_ab
df_ranks_type <- df_detects %>%
  left_join(detect_order, by = "name") %>%
  
mutate(ID = rank(-n_detects_type, ties.method = "first")) %>%
  
  mutate(
    name = factor(
      name,
      levels = unique(name[order(n_detects, -n_detects_type)])
    )
  ) %>% 
  
  mutate(
    TypeLabel = case_when(
      Type == "H_WW" ~ "Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain" ~ "Sink Biofilm",
      .default = Type
    ),
    TypeLabel = factor(TypeLabel, 
                       levels = c("Wastewater", 
                                  "Sewer Biofilm", 
                                  "Sink Biofilm")),
    type_num = as.numeric(TypeLabel)
  ) %>% 
  
  mutate(Species_label_italic = paste0("<i>", name, "</i>")) %>%  # Convert Species_label to italicized labels
  
  mutate(
    Species_label_italic = fct_reorder(Species_label_italic, n_detects, .desc = FALSE)
  ) %>% 
  
  filter(n_detects >2)

list = df_ranks_type %>% dplyr::select(name) %>% distinct() %>% pull()  


#Plot-------------
max_detect <- max(df_ranks_type$n_detects)

# candidate step
raw_step <- max_detect / 4

# function to snap to nice values (2, 5, 10, etc.)
nice_step <- function(x) {
  base <- 10^floor(log10(x))
  if (x/base <= 2) return(2 * base)
  if (x/base <= 5) return(5 * base)
  return(10 * base)
}

break_step <- 5

bars <- 
  df_ranks_type %>% 
  ggplot(aes(Species_label_italic, -n_detects_type)) +
    geom_col(aes(fill = TypeLabel),
             #color = "white",
             size = .5,
             width = 1.02) +
    # geom_curve(aes(x = 51.2, xend = 47, 
    #                y = -148, yend = -166),
    #            curvature = -.4) +
    # annotate("text", x = 47, y = -300, 
    #          label = "Each rectangle represents\none song included in the\nBBC ranking, its length\n the total points and the\ncolor indicates the rank",
    #          family = "Chivo",
    #          size = 2,
    #          lineheight = .9) +
    # annotate("text", x = 21.5, y = -120, 
    #          label = 'The Top Artists featured in the BBC´s\n"Greatest Hip-Hop Songs of All Time"',
    #          family = "Chivo",
    #          fontface = "bold",
    #          size = 4,
    #          lineheight = .9) +
    # annotate("text", x = 17, y = -120,
    #          label = 'In Autumn 2019, 108 hip-hop and music experts ranked their 5 favorites out of\n311 nominated songs in an online survey by the BBC. The graphic shows points\nscored in total and per song for the top ranked artists and broken down by era.',
    #          family = "Chivo",
    #          fontface = "bold",
    #          color = "grey30",
    #          size = 2.5,
    #          lineheight = .9) +
    coord_flip(clip = "off") +
    scale_x_discrete(position = "top") +
    scale_y_continuous(
      expand = c(.02, .02),
      limits = c(
        -max_detect - max_detect/4,
        0
      ),
      breaks = seq(
        -ceiling(max_detect / break_step) * break_step,
        0,
        by = break_step
      ),
      labels = function(x) abs(as.integer(x)),
      position = "right"
    )+
  
    # scale_fill_viridis_d(
    #                   #discrete = T, 
    #                    option = "D"
    #                       #reverse = F, 
    #                       #guide = F
    #                       )+
    # nord::scale_fill_nord(palette = "halifax_harbor", 
    #                       discrete = T, 
    #                       reverse = F, 
    #                       #guide = F
    #) +
    scale_fill_manual(values = c("#6d7ecd", "#9fdeca", "#c86c69"))+
    theme(
          #axis.text.y.right = element_blank(),
          axis.text.y.right = element_markdown(hjust = .5,size = 10), #pathogen names
          #axis.text.y.right = element_markdown(hjust = 0, size = 12),  # was hjust = .5
          axis.ticks.x = element_line(color = "black"),
          axis.ticks.length = unit(0.5, "pt"),
          axis.text.x = element_text(size = 16),
          #plot.margin = margin(5, 120, 5, 5),  # was margin(5, 0, 5, 5)
          axis.title = element_text(size = 20,face = "bold"),
          legend.position = "left",
          plot.margin = margin(5, 0, 5, 5),
          legend.text = element_text(size = 23, face = "bold")
          ) +
    labs(x = NULL, 
         fill = "",
         y = "Number of detections")

dots <- df_ranks_type %>% 

  ggplot(aes(Species_label_italic, type_num, group = Species_label_italic)) +
    geom_point(aes(Species_label_italic, 1), color = "grey75", size = 2) +
    geom_point(aes(Species_label_italic, 2), color = "grey75", size = 2) +
    geom_point(aes(Species_label_italic, 3), color = "grey75", size = 2) +
  
    geom_segment(aes(x = Species_label_italic, xend = Species_label_italic, 
                     y = 1, yend = 3), 
                 color = "grey75",
                 size = .3) +
    geom_line(color = "black",
              size = .09) +
    geom_point(aes(fill = TypeLabel, size = n_detects_type),
               shape = 21,
               color = "black",
               stroke = 1.2) +
  # geom_point(aes(fill = TypeLabel, size = n_detects_type), 
  #          shape = 21, 
  #          color = "black", 
  #          stroke = 1) +

    #geom_curve(aes(x = 47, xend = 51, 
    # geom_curve(aes(x = 45, xend = 49, 
    #                y = 6.1, yend = 4.3),
    #            curvature = .4) +
    #annotate("text", x = 45.1, y = 6.1, 
    # annotate("text", x = 43.1, y = 6.1, 
    #          label = "The dot size indicates\nthe number of songs,\nthe dot color the best\nrank in each era",
    #          family = "Chivo",
    #          size = 3.8, 
    #          lineheight = .9) +
    coord_flip() +
    scale_y_continuous(limits = c(.5, 7.3),
                       breaks = 1:3,
                       labels = c("Wastewater", 
                                  "Sewer Biofilm", 
                                  "Sink Biofilm"), 
                       position = "right") +
  scale_size(range = c(0.5, 2.6), guide = FALSE)+
  #scale_size(range = c(2, 5.5), guide = F) +
  # scale_fill_viridis_d(
  #                     #discrete = T, 
  #                      option = "D"
  #                         #reverse = F, 
  #                         #guide = F
  #                         )+
  scale_fill_manual(values = c("#6d7ecd", "#9fdeca", "#c86c69"))+
    # nord::scale_fill_nord(palette = "halifax_harbor", 
    #                       discrete = T, 
    #                       reverse = F, 
    #                       guide = F) +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 12, 
                                     #face = "bold", 
                                     hjust = .1, vjust = 0, angle = 35),
          legend.position = "none",
          axis.text.y = element_blank(),
          plot.margin = margin(5, 5, 5, 0),
          plot.caption = element_text(face = "bold", color = "grey30", 
                                      size = 15, margin = margin(t = 15))) +
    labs(x = NULL, y = NULL,
         caption = "")



plot = bars + dots + plot_layout(widths = c(1, .35))

#plot_spacer() + bars + dots + plot_layout(widths = c(1, .35), heights = c(.1, 1))


# 4. Critical step for high-res saving
showtext_opts(dpi = 300) # Synchronize showtext with ggsave

ggsave(path = "Biofilm Project Figures",file="Kraken2 CT 0.5-Upset Plot Detections for RPIP Pathogens-ONLY SHARED DATES WITH BIOFILM.png", plot=plot, width=12, height= 12,dpi = 300)
ggsave(path = "Biofilm Project Figures",file="Kraken2 CT 0.5-Upset Plot Detections for RPIP Pathogens-ONLY SHARED DATES WITH BIOFILM.svg", plot=plot, width=12, height=12,dpi = 300)

# ggsave(
#   path = "Biofilm Project Figures",
#   "Kraken2 CT 0.5-Upset Plot Detections for Pathogens-ONLY SHARED DATES WITH BIOFILM.svg",
#   plot = plot,
#   device = svglite::svglite,
#   width = 12,
#   height = 12,
#   system_fonts = list(
#     sans = "Chivo"
#   )
# )

print(plot)

#notes for figure caption: we filtered to detects to greater than 1 across all Sites. Combined detects across all Sites and reads were included that were greater than 0.0001 relative abundance

#sysfonts::font_add_google("Chivo", "Chivo", db_cache = TRUE)

#------------ Get counts of pathogens for each environment --------

presence_df <- df %>%
  mutate(
    TypeLabel = case_when(
      Type == "H_WW"   ~ "Hospital Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain"  ~ "Sink Biofilm"
    )
  ) %>%
  group_by(name, TypeLabel) %>%
  summarise(present = any(rel_ab > 0), .groups = "drop") %>%
  filter(present)

presence_wide <- presence_df %>%
  mutate(present = 1) %>%
  tidyr::pivot_wider(
    names_from  = TypeLabel,
    values_from = present,
    values_fill = 0
  )


#only ww
ww_only <- presence_wide %>%
  filter(
    `Hospital Wastewater` == 1,
    `Sewer Biofilm` == 0,
    `Sink Biofilm` == 0
  ) %>%
  pull(name)


# sewer biofilm and wastewater
sewer_ww <- presence_wide %>%
  filter(
    `Hospital Wastewater` == 1,
    `Sewer Biofilm` == 1,
    #`Sink Biofilm` == 0 |`Sink Biofilm` == 1
  ) %>%
  pull(name)


# all three

all_three <- presence_wide %>%
  filter(
    `Hospital Wastewater` == 1,
    `Sewer Biofilm` == 1,
    `Sink Biofilm` == 1
  ) %>%
  pull(name)

#sink biofilm and ww
sink_ww <- presence_wide %>%
  filter(
    `Hospital Wastewater` == 1,
    `Sewer Biofilm` == 0,
    `Sink Biofilm` == 1
  ) %>%
  pull(name)

#sink biofilm and ww
sink_only <- presence_wide %>%
  filter(
    `Hospital Wastewater` == 0,
    `Sewer Biofilm` == 0,
    `Sink Biofilm` == 1
  ) %>%
  pull(name)

#------------ Get summed counts of pathogens for each environment --------

presence_filtered_df <- df %>%
  mutate(
    TypeLabel = case_when(
      Type == "H_WW"   ~ "Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain"  ~ "Sink Biofilm"
    )
  ) %>%
  dplyr::select(SampleID,name,rel_ab,TypeLabel) %>% distinct() %>% 
  filter(rel_ab>0) %>% 
  group_by(name, TypeLabel) %>%
  summarise(present_counts = n(), .groups = "drop") %>% 
  ungroup()


presence_wide <- presence_filtered_df %>%
  tidyr::pivot_wider(
    names_from  = TypeLabel,
    values_from = present_counts,
    values_fill = 0
  )

presence_wide_filtered = presence_wide %>% 
  filter(`Wastewater`>1 & `Sewer Biofilm`>1 & `Sink Biofilm`>1)
```

#UpSet Plots (Pathogens detected in >3 sewer biofilm samples)
##upset plot set themes
```{r}
font_add_google("Reem Kufi", "Reem Kufi")   # <-- REQUIRED
showtext_auto()                             # <-- activates showtext


## ggplot theme
theme_set(theme_minimal(base_family = "Reem Kufi", base_size = 12))

theme_update(
  plot.title = element_text(size = 27,
                            face = "bold",
                            hjust = .5,
                            margin = margin(10, 0, 30, 0)),
  plot.caption = element_text(size = 9,
                              color = "grey40",
                              hjust = .5,
                              margin = margin(20, 0, 5, 0)),
  axis.text.y = element_blank(),
  axis.title = element_blank(),
  plot.background = element_rect(fill = "grey88", color = NA),
  panel.background = element_rect(fill = NA, color = NA),
  panel.grid = element_blank(),
  panel.spacing.y = unit(0, "lines"),
  strip.text.y = element_blank(),
  legend.position = "bottom",
  legend.text = element_text(size = 9, color = "grey40"),
  legend.box.margin = margin(t = 30), 
  legend.background = element_rect(color = "grey40", 
                                   size = .3, 
                                   fill = "grey95"),
  legend.key.height = unit(.25, "lines"),
  legend.key.width = unit(2.5, "lines"),
  plot.margin = margin(rep(20, 4))
)

#second-------------------

font_add_google("Chivo", "Chivo")
font_add_google("Passion One", "Passion One")

## ggplot theme
theme_set(theme_minimal(base_family = "Chivo"))

theme_update(axis.text.x = element_text(size = 11, color = "grey20"),
             axis.text.y = element_text(size = 13, color = "black", face = "bold"),
             axis.ticks.x = element_line(color = "grey45"),
             axis.ticks.y = element_blank(),
             axis.ticks.length.x = unit(.4, "lines"),
             panel.grid = element_blank(),
             plot.background = element_rect(fill = "grey60", color = "grey60"))
```
##Upset Plot - All Pathogens Presence Absence- All Sites in One-RPIP only - Narrowed down to Pathogens Detected in at least 3 biofilm samples
```{r eval=FALSE, include=FALSE}
#theme set----------

## ggplot theme
font_add_google("Reem Kufi", "Reem Kufi")   # <-- REQUIRED
showtext_auto()                             # <-- activates showtext

theme_set(theme_minimal(base_family = "Reem Kufi", base_size = 12))

theme_update(
  plot.title = element_text(size = 27,
                            face = "bold",
                            hjust = .5,
                            margin = margin(10, 0, 30, 0)),
  plot.caption = element_text(size = 9,
                              color = "grey40",
                              hjust = .5,
                              margin = margin(20, 0, 5, 0)),
  axis.text.y = element_blank(),
  axis.title = element_blank(),
  plot.background = element_rect(fill = "grey88", color = NA),
  panel.background = element_rect(fill = NA, color = NA),
  panel.grid = element_blank(),
  panel.spacing.y = unit(0, "lines"),
  strip.text.y = element_blank(),
  legend.position = "bottom",
  legend.text = element_text(size = 9, color = "white"),
  legend.box.margin = margin(t = 30), 
  legend.background = element_rect(color = "white", 
                                   size = .3, 
                                   fill = "grey95"),
  legend.key.height = unit(.25, "lines"),
  legend.key.width = unit(2.5, "lines"),
  plot.margin = margin(rep(20, 4))
)
#--------
font_add_google("Chivo", "Chivo")
font_add_google("Passion One", "Passion One")

## ggplot theme
theme_set(theme_minimal(base_family = "Chivo"))

theme_update(
              #axis.text.x = element_text(size = 11, color = "grey20"),
             #axis.text.y = element_text(size = 13, color = "black", face = "bold"),
             axis.ticks.x = element_line(color = "white"),
             axis.ticks.y = element_blank(),
             axis.ticks.length.x = unit(.4, "lines"),
             panel.grid = element_blank(),
             plot.background = element_rect(fill = "white", color = "white"))



#R prep---------------------------------------
## packages
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggtext)
library(showtext)
library(sysfonts)
library(ggplot2)
library(ggstatsplot)
library(dplyr)
library(tidyverse)
library(readxl)
library(vegan)
library(OTUtable)
library(vegan)
library(MASS)
library(ggtext)
library(ggupset)
library(forcats)

font_add_google("Chivo", "Chivo")
font_add_google("Passion One", "Passion One")

## ggplot theme
theme_set(theme_minimal(base_family = "Chivo"))

theme_update(
              #axis.text.x = element_text(size = 11, color = "grey20"),
             #axis.text.y = element_text(size = 13, color = "black", face = "bold"),
             axis.ticks.x = element_line(color = "white"),
             axis.ticks.y = element_blank(),
             axis.ticks.length.x = unit(.4, "lines"),
             panel.grid = element_blank(),
             plot.background = element_rect(fill = "white", color = "white"))

#R data -------------------------
# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame



# RPIP taxa
rpip_targets <-  rpip_pathogen_species_list

df <- kraken_merged_df_filtered %>%
  filter(taxonomy_lvl == "S") %>% 
  mutate(rel_ab = as.numeric(fraction_total_reads)) %>% 
  filter(Location %in% locations) %>% 
  filter(name %in% rpip_targets) %>%
  filter(Type %in% c("Drain", "Endcap", "H_WW")) %>% 
  mutate(
    Type = factor(Type, levels = c("Drain", "Endcap", "H_WW")),
    TypeLabel = case_when(
      Type == "H_WW" ~ "Hospital Wastewater", 
      Type == "Drain" ~ "Sink Biofilm", 
      Type == "Endcap" ~ "Sewer Biofilm", 
      TRUE ~ as.character(Type)
    ),
    Month = format(as.Date(Date), "%m-%Y")
  )

# species detected in at least 3 sewer biofilm samples
keep_pathogens <- df %>%
  group_by(name, SampleID, Type) %>%
  summarise(sum_rel_ab = sum(rel_ab, na.rm = TRUE), .groups = "drop") %>%
  mutate(detect = sum_rel_ab > 0) %>%
  filter(Type == "Endcap", detect) %>%
  group_by(name) %>%
  summarise(n_sewer_detects = n_distinct(SampleID), .groups = "drop") %>%
  filter(n_sewer_detects >= 3) %>%
  pull(name)

# keep only those pathogens everywhere downstream
df <- df %>%
  filter(name %in% keep_pathogens)


# Count detects per Species × Type
df_detects <- df %>%
  group_by(name, SampleID) %>% 
  mutate(sum_rel_ab = sum(rel_ab,na.rm = TRUE)) %>% 
  ungroup() %>% 
  
  dplyr::select(name,Type,Location,SampleID,sum_rel_ab) %>% distinct() %>%
  
  mutate(detect = ifelse(sum_rel_ab>0,1,0)) %>% 
  
  group_by(name, Type) %>%
  summarise(
    n_detects_type = sum(detect, na.rm = TRUE),   # how many samples detected?
    .groups = "drop"
  ) %>%
  ungroup() 


# Total detects per species (for ordering)
detect_order <- df_detects %>%
  group_by(name) %>%
  summarise(n_detects = sum(n_detects_type)
         #, .groups = "drop"
         ) %>%
  arrange(desc(n_detects))

# compute detection + average rel_ab
df_ranks_type <- df_detects %>%
  left_join(detect_order, by = "name") %>%
  
mutate(ID = rank(-n_detects_type, ties.method = "first")) %>%
  
  mutate(
    name = factor(
      name,
      levels = unique(name[order(n_detects, -n_detects_type)])
    )
  ) %>% 
  
  mutate(
    TypeLabel = case_when(
      Type == "H_WW" ~ "Hospital Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain" ~ "Sink Biofilm",
      .default = Type
    ),
    TypeLabel = factor(TypeLabel, 
                       levels = c("Hospital Wastewater", 
                                  "Sewer Biofilm", 
                                  "Sink Biofilm")),
    type_num = as.numeric(TypeLabel)
  ) %>% 
  
 mutate(Species_label_italic = paste0("<i>", name, "</i>")) %>%  # Convert Species_label to italicized labels

  
  mutate(
    Species_label_italic = fct_reorder(Species_label_italic, n_detects, .desc = FALSE)
  ) %>% 
  
  filter(n_detects >0)
  




#Plot-------------
max_detect <- max(df_ranks_type$n_detects)

# candidate step
raw_step <- max_detect / 4

# function to snap to nice values (2, 5, 10, etc.)
nice_step <- function(x) {
  base <- 10^floor(log10(x))
  if (x/base <= 2) return(2 * base)
  if (x/base <= 5) return(5 * base)
  return(10 * base)
}

break_step <- 25

bars <- 
  df_ranks_type %>% 
  ggplot(aes(Species_label_italic, -n_detects_type)) +
    geom_col(aes(fill = TypeLabel),
             #color = "white",
             size = .5,
             width = 1.02) +
    # geom_curve(aes(x = 51.2, xend = 47, 
    #                y = -148, yend = -166),
    #            curvature = -.4) +
    # annotate("text", x = 47, y = -300, 
    #          label = "Each rectangle represents\none song included in the\nBBC ranking, its length\n the total points and the\ncolor indicates the rank",
    #          family = "Chivo",
    #          size = 2,
    #          lineheight = .9) +
    # annotate("text", x = 21.5, y = -120, 
    #          label = 'The Top Artists featured in the BBC´s\n"Greatest Hip-Hop Songs of All Time"',
    #          family = "Chivo",
    #          fontface = "bold",
    #          size = 4,
    #          lineheight = .9) +
    # annotate("text", x = 17, y = -120,
    #          label = 'In Autumn 2019, 108 hip-hop and music experts ranked their 5 favorites out of\n311 nominated songs in an online survey by the BBC. The graphic shows points\nscored in total and per song for the top ranked artists and broken down by era.',
    #          family = "Chivo",
    #          fontface = "bold",
    #          color = "grey30",
    #          size = 2.5,
    #          lineheight = .9) +
    coord_flip(clip = "off") +
    scale_x_discrete(position = "top") +
    scale_y_continuous(
      expand = c(.02, .02),
      limits = c(
        -max_detect - max_detect/4,
        0
      ),
      breaks = seq(
        -ceiling(max_detect / break_step) * break_step,
        0,
        by = break_step
      ),
      labels = function(x) abs(as.integer(x)),
      position = "right"
    )+
    # nord::scale_fill_nord(palette = "halifax_harbor",
    #                       discrete = T,
    #                       reverse = F
    #                       #guide = F
    #                       ) +
    scale_fill_manual(values = c("#6d7ecd", "#9fdeca", "#c86c69"))+
    theme(
          axis.text.y.right = element_markdown(hjust = .5,size = 8),
          axis.text.x = element_text(size = 9),
          axis.ticks.x = element_line(color = "black"),
          axis.ticks.length = unit(0.5, "pt"),
          #axis.text.y = element_blank(),
          #axis.text.x = element_markdown(size = 35),
          axis.title = element_text(size = 14,face = "bold"),
          legend.position = "left",
          legend.text = element_text(size = 12),
          plot.margin = margin(5, 0, 5, 5)) +
    labs(x = NULL, 
         fill = "",
         y = "Number of Detections",
         title = "")


dots <- df_ranks_type %>% 

  ggplot(aes(Species_label_italic, type_num, group = Species_label_italic)) +
    geom_point(aes(Species_label_italic, 1), color = "grey75", size = 2) +
    geom_point(aes(Species_label_italic, 2), color = "grey75", size = 2) +
    geom_point(aes(Species_label_italic, 3), color = "grey75", size = 2) +
  
    geom_segment(aes(x = Species_label_italic, xend = Species_label_italic, 
                     y = 1, yend = 3), 
                 color = "grey75",
                 size = .3) +
    geom_line(color = "black",
              size = .09) +
    geom_point(aes(fill = TypeLabel, size = n_detects_type),
               shape = 21,
               color = "black",
               stroke = 1.2) +
  # geom_point(aes(fill = TypeLabel, size = n_detects_type), 
  #          shape = 21, 
  #          color = "black", 
  #          stroke = 1) +

    #geom_curve(aes(x = 47, xend = 51, 
    # geom_curve(aes(x = 45, xend = 49, 
    #                y = 6.1, yend = 4.3),
    #            curvature = .4) +
    #annotate("text", x = 45.1, y = 6.1, 
    # annotate("text", x = 43.1, y = 6.1, 
    #          label = "The dot size indicates\nthe number of songs,\nthe dot color the best\nrank in each era",
    #          family = "Chivo",
    #          size = 3.8, 
    #          lineheight = .9) +
    coord_flip() +
    scale_y_continuous(limits = c(.5, 7.3),
                       breaks = 1:3,
                       labels = c("Hospital Wastewater", 
                                  "Sewer Biofilm", 
                                  "Sink Biofilm"), 
                       position = "right") +
  scale_size(range = c(0.6, 2.7), guide = FALSE)+
  #scale_size(range = c(2, 5.5), guide = F) +
    # nord::scale_fill_nord(palette = "halifax_harbor", 
    #                       discrete = T, 
    #                       reverse = F, 
    #                       guide = F) +
  scale_fill_manual(values = c("#6d7ecd", "#9fdeca", "#c86c69"))+
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 10, hjust = .01, vjust = 0.1, angle = 45),
          #axis.text.x = element_text(size = 35, face = "bold", angle = 35),
          #axis.text.y = element_blank(),
          #axis.text.y = element_markdown(size = 20,face = "bold"),
          #axis.text.y = element_text(size = 20,face = "bold"),
          axis.text.y = element_blank(),
          plot.margin = margin(5, 5, 5, 0),
          legend.position = "none",
          plot.caption = element_text(face = "bold", color = "grey30", 
                                      size = 15, margin = margin(t = 15))) +
    labs(x = NULL, y = NULL,fill = NULL,
         caption = "")


plot = bars + dots +  plot_layout(widths = c(1, .55))

showtext_opts(dpi = 300)

ggsave(path = "Biofilm Project Figures",file="Kraken2CT0.5-Upset Plot Detections for Pathogens detected in >3 Sewer biofilm.png", plot=plot, width=10, height= 7,dpi = 300)
ggsave(path = "Biofilm Project Figures",file="Kraken2CT0.5-Upset Plot Detections for Pathogens detected in >3 Sewer biofilm.svg", plot=plot, width=10, height=7)

print(plot)



#------------ Get counts of pathogens for each environment --------

presence_df <- df %>%
  mutate(
    TypeLabel = case_when(
      Type == "H_WW"   ~ "Hospital Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain"  ~ "Sink Biofilm"
    )
  ) %>%
  group_by(name, TypeLabel) %>%
  summarise(present = any(rel_ab > 0), .groups = "drop") %>%
  filter(present)

presence_wide <- presence_df %>%
  mutate(present = 1) %>%
  tidyr::pivot_wider(
    names_from  = TypeLabel,
    values_from = present,
    values_fill = 0
  )


#only ww
ww_only <- presence_wide %>%
  filter(
    `Hospital Wastewater` == 1,
    `Sewer Biofilm` == 0,
    `Sink Biofilm` == 0
  ) %>%
  pull(name)


# sewer biofilm and wastewater
sewer_ww <- presence_wide %>%
  filter(
    `Hospital Wastewater` == 1,
    `Sewer Biofilm` == 1,
    `Sink Biofilm` == 0 |`Sink Biofilm` == 1
  ) %>%
  pull(name)

# sewer total
sewer_ww <- presence_wide %>%
  filter(
    `Sewer Biofilm` == 1
  ) %>%
  pull(name)

#only ww
ww_total<- presence_wide %>%
  filter(
    `Hospital Wastewater` == 1
  ) %>%
  pull(name)


# all three

all_three <- presence_wide %>%
  filter(
    `Hospital Wastewater` == 1,
    `Sewer Biofilm` == 1,
    `Sink Biofilm` == 1
  ) %>%
  pull(name)

#sink biofilm and ww
sink_ww <- presence_wide %>%
  filter(
    `Hospital Wastewater` == 1,
    `Sewer Biofilm` == 0,
    `Sink Biofilm` == 1
  ) %>%
  pull(name)

#------------ Get summed counts of pathogens for each environment --------

presence_filtered_df <- df %>%
  mutate(
    TypeLabel = case_when(
      Type == "H_WW"   ~ "Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain"  ~ "Sink Biofilm"
    )
  ) %>%
  dplyr::select(SampleID,name,rel_ab,TypeLabel) %>% distinct() %>% 
  filter(rel_ab>0) %>% 
  group_by(name, TypeLabel) %>%
  summarise(present_counts = n(), .groups = "drop") %>% 
  ungroup()


presence_wide <- presence_filtered_df %>%
  tidyr::pivot_wider(
    names_from  = TypeLabel,
    values_from = present_counts,
    values_fill = 0
  )

presence_wide_filtered = presence_wide %>% 
  filter(`Wastewater`>1 & `Sewer Biofilm`>1 & `Sink Biofilm`>1)
```


#Stacked Bars
##Stacked bars (facet nested by Location and Type) top 44 species by mean relative abundance colored) -Real relative abundance- x axis is unique sample type and date combination as row numbers-ONLY SHARED DATES WITH BIOFILM WITH UPDATED DATE CLASSIFICATION BY TIMEPOINT
```{r}
#specify dates
dates = c("2024-09-09", "2024-09-11", "2024-09-12", "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09")


# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame

# RPIP taxa
rpip_targets <-  rpip_pathogen_species_list


library(dplyr)
library(ggplot2)
library(forcats)
library(viridis)
library(ggh4x)

# ---------- STEP 1: Extract and clean ----------

# ---- Parameters ----
INPUT  <- "Input Data/kraken_ct_0_15_min_hit_2_all_samples.tsv"  # path to your combined k2report TSV
OUTPUT <- "Output Spreadsheets/kraken_relative_abundance.tsv"             # output file
MIN_READS <- 10                                        # minimum reads_clade to keep a species (filters noise)

# ---- Load data ----
kraken_stacked_bar <- read_tsv(INPUT, col_types = cols(
  SampleID       = col_character(),
  pct_reads      = col_double(),
  reads_clade    = col_double(),
  reads_direct   = col_double(),
  rank           = col_character(),
  taxid          = col_double(),
  name           = col_character()
))



# ---- Filter to species level and apply minimum read threshold ----

kraken_df_stacked_bar <- kraken_stacked_bar %>%
  filter(rank == "S")

# ---- Denominator = total species-level reads per sample ----
species_totals <- kraken_df_stacked_bar %>%
  group_by(SampleID) %>%
  summarise(species_total = sum(reads_clade))

# ---- Calculate relative abundance and filter noise ----
kraken_df_stacked_bar <- kraken_df_stacked_bar %>%
  left_join(species_totals, by = "SampleID") %>%
  mutate(fraction_total_reads = reads_clade / species_total) 

#get kracken_merged_df ----------------
kraken_merged_df_stacked_bar <- full_join(kraken_df_stacked_bar, metadata_df, by = 'SampleID') %>% 
  distinct() %>% 
  ungroup() %>% 
  
  dplyr::rename("taxonomy_lvl" = "rank",
                "taxonomy_id" = "taxid")


check_counts_merged = kraken_merged_df_stacked_bar %>% 
  
  filter(Type %in% c("Endcap", "H_WW","Drain")) %>% 
  
  dplyr::select("SampleID") %>% distinct()


analysis_df <- kraken_merged_df_stacked_bar %>%
  filter(
    taxonomy_lvl == "S",
    Location %in% locations,
    #name %in% rpip_targets,
    Type %in% c("Drain", "Endcap", "H_WW")
  ) %>%
  mutate(
    rel_ab = as.numeric(fraction_total_reads),
    Date   = as.Date(Date),
    Type   = factor(Type, levels = c("Drain", "Endcap", "H_WW"))
  ) %>%
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>%
  drop_na(timepoint)


biofilm_dates <- analysis_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
  mutate(
    biofilm_ref_date = Date
  ) %>% 
    ungroup() %>% 
  dplyr::distinct(Location,timepoint,biofilm_ref_date)

filtered_df <- analysis_df %>%
  mutate(Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  ) %>% 
  
  dplyr::select(-best_date,-date_diff_days,-biofilm_ref_date)

df <- filtered_df %>%
  
  
  mutate(
    Type = factor(Type, levels = c("H_WW", "Endcap","Drain")),
    TypeLabel = case_when(
      Type == "H_WW" ~ "Hospital Wastewater", 
      TRUE ~ as.character(Type)
    ),
    # Make Date a string in MM/DD/YYYY format
    DateString = format(as.Date(Date), "%m/%d/%Y"),
    Month = format(as.Date(Date), "%m-%Y")
  ) 

# ---------- New Way: Assign letters time points ----------
drain_samples <- df %>%
  filter(Type == "Drain") %>%
  distinct(Location, DateString, Type,SampleID,timepoint) %>% 
  group_by(Location, DateString) %>%
  mutate(SampleLetter = letters[seq_len(n())],
         NewDateString = paste0(timepoint, "-", SampleLetter)) %>%
  ungroup() %>%
  dplyr::select(SampleID,Location, DateString, Type, NewDateString)

# ---------- Assign letters for sewer biofilm replicates ----------
sewer_samples  <- df %>%
  filter(Type == "Endcap") %>%
  distinct(Location, DateString, Type,SampleID,timepoint) %>% 
  group_by(Location, DateString) %>%
  mutate(SampleLetter = letters[seq_len(n())],
         #NewDateString = paste0(timepoint, "-", SampleLetter)) %>%
         NewDateString = timepoint) %>%
  ungroup() %>%
  dplyr::select(SampleID,Location, DateString, Type, NewDateString)

# ---------- Assign letters for ww replicates ----------
ww_samples  <- df %>%
  filter(Type == "H_WW") %>%
  distinct(Location, DateString, Type,SampleID,timepoint) %>% 
  group_by(Location, DateString) %>%
  mutate(SampleLetter = letters[seq_len(n())],
         #NewDateString = paste0(timepoint, "-", SampleLetter)) %>%
          NewDateString = timepoint) %>%
  ungroup() %>%
  dplyr::select(SampleID,Location, DateString, Type, NewDateString)

# Combine drain and sewer samples
replicate_samples <- bind_rows(
  drain_samples,
  ww_samples,
  sewer_samples %>% dplyr::select(SampleID, Location, DateString, Type, NewDateString)
)

# ---------- Join back to main dataframe ----------
df_update <- df %>%
  left_join(replicate_samples, by = c("SampleID", "Location", "DateString", "Type")) %>%
  mutate(
    NewDateString = ifelse(is.na(NewDateString) == TRUE, DateString, NewDateString)
  )


# ---------- Collapse into top species ----------

n_species = 42

# get top 25 globally by total rel_ab
top_species <- df_update %>%
  filter(name %in% rpip_pathogen_species_list) %>% 
  group_by(name) %>%
  summarise(total_ab = sum(rel_ab, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_ab)) %>%
  slice_head(n = n_species) %>%
  pull(name)

# assign species or "Other"
plot_df <- df_update %>%
  dplyr::select(name, SampleID, Date, Location, Type, rel_ab,NewDateString) %>%
  mutate(
    Type = case_when(
      Type == "H_WW"   ~ "Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain"  ~ "Sink Biofilm",
      .default = Type
    ),
    class = case_when(
    name %in% top_species ~ name,
    name %in% detected_pathogens & !(name %in% top_species) ~ "Pathogens not Enriched",
    !(name %in% rpip_pathogen_species_list) & !(name %in% detected_pathogens) ~ "Non Pathogenic Taxa",
    name %in% rpip_pathogen_species_list & !(name %in% top_species) ~ "Other Enriched Pathogens"
  ),
  class = factor(class),
    Type  = factor(Type, levels = c("Wastewater", "Sewer Biofilm","Sink Biofilm"))
  ) %>% 
  
  # now replace Site names 
  mutate(Location = case_when(Location == "MARIO" ~ "Site 1",
                              Location == "FIONA" ~ "Site 2",
                              Location == "LUIGI" ~ "Site 3",
                              Location == "SHREK" ~ "Site 4",
                              Location == "OSCAR" ~ "Site 5",
                              .default = Location)) 

# ---------- Colors for Plot ----------

# # Assign colors: first 15 get a palette, "Other" gets grey
# ref = "Other Enriched Pathogens"
# ref2 = "Pathogens not Enriched"
# ref3 = "Non Pathogenic Taxa"
# myColors <- fortyfive_pal[1:length(levels(plot_df$class))]
# names(myColors) <- levels(plot_df$class)
# myColors[names(myColors)==ref] <- "grey"
# myColors[names(myColors)==ref2] <- "grey30"
# myColors[names(myColors)==ref3] <- "black"


##actual plotting-----------
library(dplyr)
library(ggplot2)
library(ggh4x)
library(ggtext)   # for element_markdown

# # 1) Build a per-facet ordering of samples (one row per unique sample)
# sample_order <- plot_df %>%
#   distinct(Location, Type, SampleID, Date) %>%
#   group_by(Location, Type) %>%
#   arrange(Date, SampleID, .by_group = TRUE) %>%
#   mutate(
#     Event = row_number(),                                # 1..N per facet (no gaps)
#     EventLabel = format(as.Date(Date), "%Y-%m-%d")       # optional: a label you can show
#   ) %>%
#   ungroup()
# 
# # 2) Join that back so every row has its facet-local Event index
# plot_df2 <- plot_df %>%
#   left_join(sample_order, by = c("Location","Type","SampleID","Date"))

category_levels <- c("Other Enriched Pathogens",
                     "Pathogens not Enriched",
                     "Non Pathogenic Taxa")







# Plot using Event (factor) on x → no gaps within each facet


library(scales)  # for percent_format()

plot_df_updated <- plot_df %>%
  mutate(
    Type = dplyr::recode(Type,
                  "Sink Biofilm" = "Sink\nBiofilm",
                  "Sewer Biofilm" = "Sewer\nBiofilm"),
    class = case_when(
      class %in% c("Pathogens not Enriched", "Non Pathogenic Taxa") ~ "Non Enriched Taxa",
      TRUE ~ as.character(class)
    ),
    class = factor(class)
  )%>% 
  
  filter(class != "Non Enriched Taxa", class != "Other Enriched Pathogens")

# ---------- Colors for Plot ----------

# Assign colors: first 15 get a palette, "Other" gets grey
ref = "Other Enriched Pathogens"
ref2 = "Non Enriched Taxa"
myColors <- fortyfive_pal[1:length(levels(plot_df_updated$class))]
names(myColors) <- levels(plot_df_updated$class)
myColors[names(myColors)==ref] <- "grey"
myColors[names(myColors)==ref2] <- "grey20"


  
# italicize legend labels
species_labels <- setNames(
  ifelse(levels(plot_df_updated$class) %in% category_levels,
         levels(plot_df_updated$class),                       # plain text for categories
         paste0("*", levels(plot_df_updated$class), "*")),    # markdown italics for species
  levels(plot_df_updated$class)
)


plot = ggplot(plot_df_updated, aes(x = NewDateString, y = rel_ab, fill = class)) +
  geom_col(position = "stack") +
  facet_nested(Location ~ Type, scales = "free", space = "free_x", switch = "y") +
  scale_fill_manual(values = myColors,
                    labels = species_labels) +
  # scale_y_continuous(
  #   breaks = c(0, 0.25, 0.50, 0.75),   # tick marks at 25, 50, 75
  #   labels = c("0","25", "50", "75"),   # plain numbers, no %
  #   limits = c(0, 1)                 # optional: keep to 100%
  # ) +
  labs(x = NULL, y = "Fraction of Annotated RPIP Reads", fill = "", title = "") +
theme_minimal(base_size = 11) +
  guides(fill  = guide_legend(ncol = 6))+
theme(
  panel.background = element_rect(fill = "white", colour = "black"),
  strip.background = element_rect(fill = "white", colour = "black"),
  strip.text = element_text(
    face = "bold", 
    hjust = 0.5,          # center horizontally
    vjust = 0.5,          # center vertically
    lineheight = 0.9      # tighter line spacing if wrapped
  ),
  axis.text.x        = element_blank(),
  panel.spacing      = unit(0.03, "lines"),
  legend.position    = "bottom",
  legend.key.size    = unit(0.3, "cm"),
  legend.text        = ggtext::element_markdown(size = 5.5),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  #panel.grid.major.y = element_line(color = "grey80"),
  panel.grid.minor.y = element_blank(),
  panel.ontop        = FALSE   # <- ensures geoms (bars) are drawn over grid lines
)


print(plot)

  ggsave(path = "Biofilm Project Figures",file="Kraken CT 0.5 stacked bar plot of pathogens detected Overall-ONLY SHARED DATES WITH BIOFILM .png", plot=plot, width=8.5, height=7.5,dpi = 400)
  ggsave(path = "Biofilm Project Figures",file="Kraken CT 0.5 stacked bar plot of pathogens detected Overall-ONLY SHARED DATES WITH BIOFILM .svg", plot=plot, width=8.5, height=7.5)
  
#for slides
  
plot = ggplot(plot_df, aes(x = NewDateString, y = rel_ab, fill = class)) +
  geom_col(position = "stack") +
  facet_nested(Location ~ Type, scales = "free", space = "free_x", switch = "y") +
  scale_fill_manual(values = myColors,
                    labels = species_labels) +
  # scale_y_continuous(
  #   breaks = c(0, 0.25, 0.50, 0.75),   # tick marks at 25, 50, 75
  #   labels = c("0","25", "50", "75"),   # plain numbers, no %
  #   limits = c(0, 1)                 # optional: keep to 100%
  # ) +
  labs(x = NULL, y = "Fraction of Annotated RPIP Reads", fill = "", title = "Kraken Species Detected") +
theme_minimal(base_size = 11) +
theme(
  panel.background = element_rect(fill = "white", colour = "black"),
  strip.background = element_rect(fill = "white", colour = "black"),
  strip.text = element_text(
    face = "bold", 
    hjust = 0.5,          # center horizontally
    vjust = 0.5,          # center vertically
    lineheight = 0.9      # tighter line spacing if wrapped
  ),
  axis.text.x        = element_blank(),
  panel.spacing      = unit(0.05, "lines"),
  legend.position    = "bottom",
  legend.key.size    = unit(0.1, "cm"),
  legend.spacing.x = unit(0.2, "cm"),   # space between items horizontally
  legend.spacing.y = unit(0.2, "cm"),    # space between rows if wrapped
  legend.text        = ggtext::element_markdown(size = 5),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  legend.key.width  = unit(1, "cm"),
  legend.key.height = unit(0.5, "cm"),
  #panel.grid.major.y = element_line(color = "grey80"),
  panel.grid.minor.y = element_blank(),
  panel.ontop        = FALSE   # <- ensures geoms (bars) are drawn over grid lines
)+
  guides(fill = guide_legend(ncol = 4))


print(plot)

  ggsave(path = "Biofilm Project Figures",file="Kraken CT 0.5 stacked bar plot of pathogens detected-slides-ONLY SHARED DATES WITH BIOFILM .png", plot=plot, width=9.5, height=8,dpi = 600)
  ggsave(path = "Biofilm Project Figures",file="Kraken CT 0.5 stacked bar plot of pathogens detected-slides-ONLY SHARED DATES WITH BIOFILM .svg", plot=plot, width=9.5, height=8)
  
  
# now vertical:
  
  plot = ggplot(plot_df, aes(y = NewDateString, x = rel_ab, fill = class)) +
  geom_col(position = "stack") +
  facet_nested(Type ~ Location, scales = "free", space = "free_y",switch = "y") +
  scale_fill_manual(values = myColors,
                    labels = species_labels,
                    guide = guide_legend(
    nrow = 4,              # number of rows in the legend
    byrow = TRUE,          # fill across rows instead of down columns
    title.position = "top" # keep the title above
  )
  ) +
  # scale_x_continuous(
  #   breaks = c(0, 0.25, 0.50, 0.75),   # tick marks at 25, 50, 75
  #   labels = c("0","25", "50", "75"),  # plain numbers, no %
  #   limits = c(0, 1)
  # ) +
  labs(y = NULL, x = "Fraction of Annotated RPIP Reads", fill = "") +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(
      face = "bold", 
      hjust = 0.5,
      vjust = 0.5,
      lineheight = 0.9
    ),
    axis.text.y        = element_blank(),  # hide if crowded
    panel.spacing      = unit(0.05, "lines"),
    legend.position    = "bottom",
    legend.key.size    = unit(0.5, "cm"),
    legend.text        = ggtext::element_markdown(size = 6),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.ontop        = FALSE
  )

  
  print(plot)

  ggsave(path = "Biofilm Project Figures",file="Kraken CT 0.5 Bracken stacked bar plot of pathogens detected Overall -vertical-ONLY SHARED DATES WITH BIOFILM .png", plot=plot, width=8, height=10,dpi = 500)
  ggsave(path = "Biofilm Project Figures",file="Kraken CT 0.5 Bracken stacked bar plot of pathogens detected in Overall-vertical-ONLY SHARED DATES WITH BIOFILM .svg", plot=plot, width=8, height=9)


library(dplyr)

# Create summary table of percent contribution by pathogen and sample type
pathogen_summary <- plot_df_updated %>%
  group_by(class, Type) %>%
  summarise(
    MinPercent = min(rel_ab, na.rm = TRUE) * 100,
    MaxPercent = max(rel_ab, na.rm = TRUE) * 100,
    MeanPercent = mean(rel_ab, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  dplyr::rename(Pathogen = class, SampleType = Type) %>%
  arrange(SampleType, desc(MeanPercent))

# View the table
pathogen_summary


```





#MaAsLin Pathogens
##Maaslin3- Pathogens 
```{r}

#specify dates-------------
dates <- as.Date(c(
  "2024-09-09", "2024-09-11", "2024-09-12",
  "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09"
))

# Specify locations to include--------
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame

# RPIP taxa
rpip_targets <-  rpip_pathogen_species_list

#libraries--------
library(ggplot2)
library(readxl)
library(here)
library(cowplot)
library(tidyverse)
library(patchwork)
library(readr)
library(phyloseq)
library(microshades)
library(ape)
library(vegan)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(maaslin3)
library(writexl)
library(ggrepel)

#maaslin3----------------
df <- kraken_merged_df_filtered %>%
  filter(
    taxonomy_lvl == "S",
    Location %in% locations,
    name %in% rpip_targets,
    Type %in% c("Drain", "Endcap", "H_WW")
  ) %>%
  mutate(
    rel_ab = as.numeric(fraction_total_reads),
    Date   = as.Date(Date),
    Type   = factor(Type, levels = c("H_WW", "Endcap" ,"Drain"))
  ) %>%
  
  #clean ontology df
  dplyr::select(name,SampleID, Date, Location, Type,rel_ab) %>% 

  #filter for sample type
  filter(Type %in% c("Endcap", "H_WW")) 
  # 
  # #filter specific dates shared with biofilm
  # filter(
  #   (Type == "H_WW" & Date %in% dates) |
  #   (Type != "H_WW")
  # ) 
  

df_sum <- df %>%
  group_by(SampleID, name) %>%
  summarise(
    rel_ab = sum(rel_ab, na.rm = TRUE),
    .groups = "drop"
  )

otu <- df_sum %>%
  dplyr::select(SampleID, name, rel_ab) %>%
  tidyr::pivot_wider(
    names_from = name,
    values_from = rel_ab,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("SampleID")

meta <- df %>%
  dplyr::select(
    SampleID,
    Date,
    Location,
    Type
  ) %>%
  distinct() %>%
  tibble::column_to_rownames("SampleID")


#MaAsLin3------------------------

fit_data <- maaslin3(
  input_data = otu,
  input_metadata = meta,
  output = "maaslin3_output",
  fixed_effects = c("Type"),
  random_effects = c("Location"),
  normalization = "NONE",
  transform = "LOG",
  standardize = TRUE,
  min_prevalence = 0.1,
  verbosity = "ERROR"
)



#plotting it-------------
plot_df <- fit_data$normalized %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "feature",
    values_to = "abundance"
  ) %>%
  left_join(
    fit_data$metadata %>%
      rownames_to_column("SampleID"),
    by = "SampleID"
  )


#volcano plot----------

#get mean abundance directly from MaAslin3 normalized data
# Significant features
sig_abund_joint <- subset(fit_data$fit_data_abundance$results, qval_joint < 0.05)
sig_prev_joint <- subset(fit_data$fit_data_prevalence$results, qval_joint < 0.05)

sig_feats <- union(sig_abund_joint$feature, sig_prev_joint$feature)

top_feat <- sig_abund_joint %>%
  filter(feature %in% sig_feats) %>%
  arrange(qval_joint)

# Compute mean abundance directly from MaAsLin3 normalized table
mean_abund <- colMeans(fit_data$normalized[, top_feat$feature, drop = FALSE], na.rm = TRUE)
top_feat$mean_abund <- mean_abund

# Optionally join extra info if you have a taxonomy table
# top_feat <- left_join(top_feat, top_feat_unique, by = "feature")

# Example: volcano for a subset of features (adjust as needed)
b_rg <- top_feat %>%
  # filter(name.x == "growthregrowth") %>% # if you have a variable for before/regrowth
  filter(!is.na(feature))  # just keep valid features

# Label a few genera if you have that info
# b_rg$label_it <- ifelse(b_rg$genus %in% c("Cupriavidus","Ralstonia","Burkholderia","Pseudomonas","Alcaligenes","Stenotrophomonas"), b_rg$genus, NA)

poster_colors = c("#dc9298","#a51c2f","#3d3d3d","#f3e3e5","#417c8c","#e57262","#915493")

poster_colors = c("#ce2e2d","#c00300","#88b8c5","#dd9298","#ab74a8","#fc3b3b","#a61d30","#eededf","#dc6c64","#467d91")


volc <- ggplot(b_rg, aes(
    x = coef,
    y = -log10(qval_joint),
    color = coef > 0,
    size = mean_abund
)) +
  geom_point() +
  # geom_text(aes(label = label_it), vjust = -1, size = 3, na.rm = TRUE, show.legend = FALSE, fontface = "italic") +
  scale_color_manual(values =c("#6d7ecd", "#c86c69"),labels = c("Biofilm","Wastewater")) +
  scale_size_continuous(range = c(1, 5)) +
  labs(
    title = "",
    x = "Effect Size (MaAsLin3 coefficient)",
    y = "-log10(q-value)",
    color = "",
    size = "Mean Rel. Abundance"
  ) +
  theme_minimal() +
  coord_cartesian(clip = "off")

volc

#get top pahogens-------

#label top args in new volcano plot----------
# Get top pathogens (e.g., top 5 by q-value)

top_pathogens = b_rg %>% 
  #filter(-log10(qval_joint)>150) %>% 
  filter(-log10(qval_joint)>8 | (-log10(qval_joint)>5 & coef > 3) | mean_abund > 0.06) %>% 
  dplyr::select(feature) %>% distinct() %>% pull(); top_pathogens

top_pathogens <- b_rg %>%
  arrange(qval_joint) %>%
  #slice(1:15) %>%        # top 15 most significant
  pull(feature)


# Add a column for labels in the volcano data
b_rg$label_it <- ifelse(b_rg$feature %in% top_pathogens, b_rg$feature, NA)

# Volcano plot with top ARG labels
volc <- ggplot(b_rg, aes(
    x = coef,
    y = -log10(qval_joint),
    color = coef > 0,
    size = mean_abund
)) +
  geom_point() +
    geom_text_repel(
    aes(label = label_it),
    size = 3,                 # smaller text
    max.overlaps = Inf,         # don’t drop labels
    box.padding = 0.3,          # space around text
    point.padding = 0.3,        # space from points
    min.segment.length = 0,     # always draw connecting line
    segment.size = 0.1,
    na.rm = TRUE,
    show.legend = FALSE,
    fontface = "italic"
  ) +
  #geom_text(aes(label = label_it), vjust = -1, size = 3, na.rm = TRUE, show.legend = FALSE, fontface = "italic") +
  scale_color_manual(values = c("#6d7ecd", "#c86c69"),labels = c("Biofilm","Wastewater")) +
  scale_size_continuous(range = c(1, 5)) +
  labs(
    title = "",
    x = "Effect Size (MaAsLin3 coefficient)",
    y = "-log10(q-value)",
    color = "",
    size = "Mean Rel. Abundance"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey88"),  # <--- major grid color
    panel.grid.minor = element_line(color = "grey88")  # <--- minor grid color
    #panel.border = element_rect(color = "", fill = NA)  # optional border
  ) +
  coord_cartesian(clip = "off")

volc

# Save final heatmap
ggsave(
  plot = volc,
  filename = "Biofilm Project Figures/Kraken 0.5 - MaAsLin3 Volcano Plot by Pathogen.png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)

# Save final heatmap
ggsave(
  plot = volc,
  filename = "Biofilm Project Figures/Kraken 0.5 -MaAsLin3 Volcano Plot by Pathogen.svg",
  width = 8,
  height = 6,
  units = "in",
  dpi = 500
)

# Get plot with all pathogens for visualizing---

all_pathogens <- b_rg %>%
  arrange(qval_joint) %>%
  #slice(1:15) %>%        # top 15 most significant
  pull(feature)


# Add a column for labels in the volcano data
b_rg$label_it <- ifelse(b_rg$feature %in% all_pathogens, b_rg$feature, NA)

# Volcano plot with top ARG labels
volc <- ggplot(b_rg, aes(
    x = coef,
    y = -log10(qval_joint),
    color = coef > 0,
    size = mean_abund
)) +
  geom_point() +
    geom_text_repel(
    aes(label = label_it),
    size = 3,                 # smaller text
    max.overlaps = Inf,         # don’t drop labels
    box.padding = 0.3,          # space around text
    point.padding = 0.3,        # space from points
    min.segment.length = 0,     # always draw connecting line
    segment.size = 0.1,
    na.rm = TRUE,
    show.legend = FALSE,
    fontface = "italic"
  ) +
  #geom_text(aes(label = label_it), vjust = -1, size = 3, na.rm = TRUE, show.legend = FALSE, fontface = "italic") +
  scale_color_manual(values = c("#6d7ecd", "#c86c69"),labels = c("Biofilm","Wastewater")) +
  scale_size_continuous(range = c(1, 5)) +
  labs(
    title = "",
    x = "Effect Size (MaAsLin3 coefficient)",
    y = "-log10(q-value)",
    color = "",
    size = "Mean Rel. Abundance"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey88"),  # <--- major grid color
    panel.grid.minor = element_line(color = "grey88")  # <--- minor grid color
    #panel.border = element_rect(color = "", fill = NA)  # optional border
  ) +
  coord_cartesian(clip = "off")

volc

# Save final heatmap
ggsave(
  plot = volc,
  filename = "Biofilm Project Figures/CARD-MaAsLin3 Volcano Plot by Pathogen All Pathogens Labeled.png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)



#filtering for boxplots-----

sig_feats <- union(
  sig_abund_joint$feature,
  sig_prev_joint$feature
)

boxplot_df <- plot_df %>%
  filter(feature %in% sig_feats)

#join MaAsLin3 stats for box plots (effect sizes, q-values) -------
top_feat <- sig_abund_joint %>%
  filter(feature %in% sig_feats) %>%
  arrange(qval_joint)

boxplot_joined_df <- boxplot_df %>%
  left_join(
    top_feat %>% select(feature, coef, qval_joint),
    by = "feature"
  )


#loop through top pathogens------------

for (pathogen_var in top_pathogens) {


plot <- boxplot_joined_df %>% 
  filter(feature == pathogen_var) %>% 
  ggplot(., aes(x = Type, y = abundance, colour = Type)) +
  theme_minimal() +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(y = paste0(pathogen_var," Abundance"))

print(plot)

}

#create arg box plots using faceting approach

library(ggpubr)

plot <- boxplot_joined_df %>% 
  filter(feature %in% top_pathogens) %>% 
  ggplot(aes(x = Type, y = abundance)) +
  
  geom_boxplot(
    outlier.shape = NA,
    color = "black",
    width = 0.6
  ) +
  
  geom_jitter(
    aes(color = Type),
    width = 0.15,
    size = 1.2,
    alpha = 0.6
  ) +
  
  stat_compare_means(
    method = "kruskal.test",
    label = "p.format",
    size = 3
  ) +
  
  facet_wrap(
    ~ feature,
    scales = "free_y",
    ncol = 4
  ) +
  
  labs(
    x = "",
    y = "Abundance"
  ) +
  
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 5),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "italic", size = 9),
    legend.position = "bottom"
  )

plot

# Save final heatmap
ggsave(
  plot = plot,
  filename = "Biofilm Project Figures/Kraken 0.5 --Box Plots Faceted of Top Pathogens.png",
  width = 8,
  height = 12,
  units = "in",
  dpi = 300
)




#stats-------------
# 
# # Cupriavidus
# stat.test <- plot_df %>%
#   filter(genus == "Cupriavidus") %>% 
#   pairwise_wilcox_test(
#     abundance ~ growth,
#     p.adjust.method = "fdr")
# 
# c_plot <- cup_plot + stat_pvalue_manual(stat.test, label = "p.adj", y.position = c(0.55, 0.65, 0.75))


#MaAsLin3 table------

# mas_table <- test %>% 
#   select(name.x, qval_joint.x, mean_abund, genus) %>% 
#   filter(qval_joint.x > 0) %>% 
#   arrange(qval_joint.x) %>%  
#   slice(1:21)
# 
# write_xlsx(mas_table, "maaslin3_table.xlsx")



#------Extract All ARGs in Table Format by Mean Abundance, Coefficient, 
ontologies = df %>% 
  
  dplyr::select(name)

pathogens = b_rg %>% 
  dplyr::select(feature, coef,N, mean_abund,qval_joint) %>% 
  
  mutate(enriched = ifelse(coef>0, "Wastewater", "Biofilm")) %>% 
  
  dplyr::rename("name" = "feature") %>% 
  
  left_join(ontologies,by = c("name")) %>% distinct()

#------Stacked Bar of What's Enriched in WW vs Biofilm----------------

# stacked_bar =  pathogens %>% 
#   
#   group_by(enriched) %>% 
#   mutate(total = n()) %>% 
#   ungroup() %>% 
#   
#   
#   group_by(enriched,card_edited_drug_class) %>% 
#   mutate(count = n()) %>% 
#   dplyr::select(enriched,card_edited_drug_class,total,count) %>% distinct() %>% 
#   summarise(percent = count/total * 100) %>% 
#   ungroup() 
#   
# 
#   
#   
#   
# ggplot(stacked_bar, aes(x = factor(enriched), y = percent, fill = card_edited_drug_class)) +
#   geom_col(position = "stack") 
  


```
 
 




#  SET WORKING DIRECTORY #######################################

```{r}
#Setting working directory for importing data from excel
setwd("/Users/amd689/Documents/Healy Lab/Biofilm Project/Biofilm YNHH Project R")
```


# SET NAME OF PATH TO DATA
```{r}
data_path = "Input Data/"

```
 

#LOAD LIBRARIES
```{r}
rm(list = ls())
library(ggplot2)
library(ggtext)
library(ggstatsplot)
library(dplyr)
library(tidyverse)
library(readxl)
library(vegan)
library(OTUtable)
library(vegan)
library(ggupset)
library(MASS)
library(ggh4x)
library(ggplot2)
library(ggstatsplot)
library(dplyr)
library(tidyverse)
library(readxl)
library(vegan)
library(OTUtable)
library(vegan)
library(MASS)
library(sysfonts)
library(showtext)
library(tidytext)
library(ggtext)
```

#INSTALL FUNCTION FROM GITHUB
```{r eval=FALSE, include=FALSE}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biobakery/maaslin3")

Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)


install.packages("BiocManager")
#BiocManager::install("maaslin3")
BiocManager::install(c("phyloseq", "vegan", "Biostrings"))

remotes::install_github("KarstensLab/microshades", dependencies = TRUE)

#function pulled from github

#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(devtools)

#devtools::install_github("gauravsk/ranacapa", force = TRUE)

# library(ranacapa)
# 
# data(iris)
# ranacapa::pairwise_adonis(iris[, 1:4], iris$Species)
# 
# ranacapa::pairwise_adonis(iris[, 1:4], iris$Species, reduce = 'setosa')
# 
# # similarity euclidean from vegdist and holm correction
# pairwise_adonis(x = iris[, 1:4], factors = iris$Species,
# sim_method = 'euclidian', p_adjust_m = 'holm')
# 
# #similarity manhattan from daisy and bonferroni correction
# pairwise_adonis(x = iris[, 1:4], factors = iris$Species,
# # sim_method = 'manhattan', p_adjust_m = 'bonferroni')
# 
# 
# install.packages("spaa")
# library(spaa)
# 
# devtools::install_github("GuillemSalazar/EcolUtils")
# 
# library(EcolUtils)
# 
# data(dune)
# data(dune.env)
# adonis.pair(vegdist(dune),dune.env$Management)
```

#COLORS
```{r}

colors = c("red",  "blue",  "green", "purple", "brown", "gold1", "orange", "hotpink", "black", "seagreen", "lightblue", "pink", "goldenrod4")

colors = c("#9FB498", "#8488AC","#BA6646","#C8DFEA","#565656","#865856","#E8DCDD","#F8E5DC","#E8DCDC","#B99796","black","seagreen")



resistome_colors = c("#8488AC", "#9FB498", "#BA6646", "#565656", "#865856", "#B99796", "#F8E5DC")

posterCol1=  "#8488AC"
posterCol2 = "#9FB498"
posterCol3 = "#BA6646"
posterCol4 = "#565656"
posterCol5 = "#865856"
posterCol6 = "#B99796"
posterCol7 = "#F8E5DC"

viridisCol1 = "#79307D"
viridisCol2 = "#417C8C"
viridsCol3 = "#E57262"

colors = c( "#417C8C","#417C8C","#417C8C","#417C8C","#E57262","#79307D","#79307D","#79307D")

#colors = c( "#79307D","#79307D","#417C8C","#417C8C","#417C8C","#417C8C","#E57262","#79307D","#79307D","#79307D")

colors = c( "#8488AC","#BA6646","#9FB498")

fortyfive_pal = c("#a0cb6b","#8368cb","#c86c69","#cdd3e5","#dab594","#d692d1",
               "#7495c3","#9fdeca","#e2e8c3","#d8a4af","#71bed2","#bca9dd",
               "#8bb598","#e5cbd4","#6d7ecd","#e4d5ca","#a8dfa5","#a0bada",
               "#cbca6f","#c6926f","#cce7e6","#81a4b2","#ca69c3","#76bd75",
               "#d37dae","#abb684","#6ecda4","#a88bbb","#8cbbb6","#d18191",
               "#c1d8c2","#dac4e6","#e0adce","#a96dca","#e2bdb6","#aacedd",
               "#9992d9","#e0d0ab","#a2abdf","#b88b9e","#b4a382","#dba294",
               "#c9dd9b","#c8b67c","#b9847b")


poster_colors = c("#dc9298","#a51c2f","#3d3d3d","#f3e3e5","#417c8c","#e57262","#915493")


poster_colors = c("#dc9298","#a51c2f","#3d3d3d","#f3e3e5","#417c8c","#e57262","#915493")

poster_colors = c("#ce2e2d","#c00300","#88b8c5","#dd9298","#ab74a8","#fc3b3b","#a61d30","#eededf","#dc6c64","#467d91")



```



#SET THEMES
```{r}
library(ggplot2)

# Reset global theme and font
theme_set(theme_grey(base_family = "", base_size = 11))

# Remove theme_update() overrides
theme_update(
  axis.text.x        = element_text(),
  axis.text.y        = element_text(),
  axis.ticks.x       = element_line(),
  axis.ticks.y       = element_line(),
  axis.ticks.length.x = unit(0.15, "lines"),
  panel.grid         = element_line(),
  plot.background    = element_rect()
)

# Reset geom defaults that may have inherited "Chivo"
update_geom_defaults("text",  list(family = ""))
update_geom_defaults("label", list(family = ""))

showtext::showtext_auto(FALSE)

```

#FUNCTIONS
```{r}

pairwise.adonis <- function(x,factors, sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  co = combn(unique(factors),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(co[1,elem],co[2,elem]),] ~ factors[factors %in% c(co[1,elem],co[2,elem])] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}


#updated adonis2 function

pairwise.adonis2 <- function(x, factors, sim.method = 'bray', p.adjust.m = 'bonferroni') {
  library(vegan)
  
  co <- combn(unique(factors), 2) # all pair combinations
  pairs <- c()
  F.Model <- c()
  R2 <- c()
  p.value <- c()
  
  for (elem in 1:ncol(co)) {
    # Subset the data & factor
    sel <- factors %in% c(co[1, elem], co[2, elem])
    
    # Skip if any group has < 2 samples
    if (min(table(factors[sel])) < 2) {
      next
    }
    
    ad <- adonis2(
      x[sel, ] ~ factors[sel],
      method = sim.method
    )
    
    pairs     <- c(pairs, paste(co[1, elem], 'vs', co[2, elem]))
    F.Model   <- c(F.Model, ad$F[1])
    R2        <- c(R2, ad$R2[1])
    p.value   <- c(p.value, ad$`Pr(>F)`[1])
  }
  
  p.adjusted <- p.adjust(p.value, method = p.adjust.m)
  pairw.res  <- data.frame(pairs, F.Model, R2, p.value, p.adjusted)
  
  return(pairw.res)
}


```

#Load CSVs
```{r}
#load df
card_df <- read.table(
  "Input Data/combined_card_rpkm.tsv",
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE
) %>% 
  
  mutate(ARO.Accession = sub(".*(ARO:[^|]+)\\|.*", "\\1", gene_id))

aro_index <- read.table(
  "Input Data/aro_index.tsv",
  header = TRUE,
  sep = "\t",
  quote = "",
  fill = TRUE
)



card_df %>%
  group_by(SampleID) %>%
  summarise(total_RPKM = sum(RPKM))

```


#Load in Metadata
```{r}

metadata_df = read.csv("Input Data/metadata_Gi_AD.csv",header = T) %>% 
  
  #change name of Location
  mutate(Location =  sub("^[0-9]+\\.", "", Location)) %>% 
  
  #change month column
  mutate(Month = sub("^[0-9]+_", "", Month)) %>% 

  #change date column
  mutate(Date = sub("^[0-9]+_", "", Date)) %>% 
  mutate(Date = mdy(paste0(Date, " 2024"))) %>% 
  
  dplyr::select(-LocationMonth,-LocationType) %>% 
  
  dplyr::rename("SampleID"="Sample_id")

```


#Merge data
```{r}
# Merge datasets based on the 'Sample_ID' column
card_merged_df <- full_join(card_df, metadata_df, by = 'SampleID') %>% 
  distinct() %>% 
  ungroup() 
    
```


#Get ARGs in panel
```{r}
#load df with ARGs in RPIP panel
rpip_amr_df <- read_xlsx(
  "Input Data/RPIP_Panel_Target_List (1).xlsx",
  sheet = "Bacterial AMR",
  skip = 0
) 

rpip_alleles = rpip_amr_df %>% dplyr::select(Allele) %>% distinct(Allele) %>% pull(Allele)

rpip_gene_families = rpip_amr_df %>% distinct(`Gene Family`) %>% pull(`Gene Family`)

```

#Get proper ontologies
```{r}


# Merge datasets based on the 'Sample_ID' column
card_merged_clean_df <- card_merged_df %>% 
  distinct() %>% 
  ungroup() %>%
  mutate(Accession = str_extract(gene_id, "(?<=\\|)[^|]+(?=\\|ARO)"),
         `ARO Accession` = str_extract(gene_id, "ARO:\\d+"),
         Gene = str_extract(gene_id, "[^|]+$"))

library(jsonlite)
library(dplyr)
library(purrr)
library(tidyr)


#library(dplyr)
library(readr)

# Load
aro_df  <- read_tsv("/Users/amd689/Documents/Healy Lab/Biofilm Project/Biofilm YNHH Project R/Input Data/card-data/aro_index.tsv")             # ARO terms
aro_df2 <- read_tsv("/Users/amd689/Documents/Healy Lab/Biofilm Project/Biofilm YNHH Project R/Input Data/card-data/aro_categories.tsv")        # categories
aro_df3 <- read_tsv("/Users/amd689/Documents/Healy Lab/Biofilm Project/Biofilm YNHH Project R/Input Data/card-data/aro_categories_index.tsv")  # links

# Load the JSON file
aro <- fromJSON("/Users/amd689/Documents/Healy Lab/Biofilm Project/Biofilm YNHH Project R/Input Data/card-ontology/aro.json", flatten = TRUE) %>% 
  dplyr::rename( `ARO Accession` = accession)

aro_short_names = aro %>% 
  dplyr::select(`ARO Accession`,`CARD Short Name`)


#add short name
# card_merged_clean_df = card_merged_clean_df %>% 
#   full_join(aro_short_names,by = "accession")


# Load
card_to_megares_df  <- read_csv("/Users/amd689/Documents/Healy Lab/Biofilm Project/Biofilm YNHH Project R/Input Data/megares_to_external_header_mappings_v3.00.csv")  %>% dplyr::select(-UpdatedHeader) %>% dplyr::rename("header" = "MEGARes_header")
megares_annotations <- read_csv("/Users/amd689/Documents/Healy Lab/Biofilm Project/Biofilm YNHH Project R/Input Data/megares_annotations_v3.00.csv")        # categories

megares_to_card_annotation = megares_annotations %>%
  inner_join(card_to_megares_df, by =  "header") %>% filter(Database == "CARD") %>%
  mutate(`ARO Accession` = sub(".*(ARO:[0-9]+).*", "\\1", Source_header))

#final clean df to work with
card_merged_w_ont_df = card_merged_clean_df %>% 
  full_join(aro_df, by = "ARO Accession") %>% 
  left_join(megares_to_card_annotation, by = "ARO Accession") %>% 
  mutate(
    card_edited_drug_class = case_when(
      # Core antibiotic classes
      str_detect(tolower(`Drug Class`), "penicillin|carbapenem|cephalosporin|betalactam") ~ "Betalactams",
      str_detect(tolower(`Drug Class`), "sulfonamide") ~ "Sulfonamides",
      str_detect(tolower(`Drug Class`), "diaminopyrimidine") ~ "Trimethoprim",
      str_detect(tolower(`Drug Class`), "tetracycline") ~ "Tetracyclines",
      str_detect(tolower(`Drug Class`), "aminocoumarin") ~ "Aminocoumarin",
      str_detect(tolower(`Drug Class`), "aminoglycoside") ~ "Aminoglycoside",
      str_detect(tolower(`Drug Class`), "elfamycin") ~ "Elfamycin",
      str_detect(tolower(`Drug Class`), "phosphonic acid") ~ "Fosfomycin",
      str_detect(tolower(`Drug Class`), "glycopeptide") ~ "Glycopeptides",
      str_detect(tolower(`Drug Class`), "macrolide|lincosamide|streptogramin") ~ "MLS",
      str_detect(tolower(`Drug Class`), "peptide") ~ "Peptides",
      str_detect(tolower(`Drug Class`), "phenicol") ~ "Phenicol",
      str_detect(tolower(`Drug Class`), "fluoroquinolone|quinolone") ~ "Quinolones",
      str_detect(tolower(`Drug Class`), "rifamycin") ~ "Rifamycines",
      
      # Multi-drug / combination
      str_detect(`Drug Class`, ";") ~ "Multi-Drug Resistant",   # any cell with semicolon
      
      # Less common or rare antibiotics
      str_detect(tolower(`Drug Class`), "nucleoside") ~ "Nucleosides",
      str_detect(tolower(`Drug Class`), "nitroimidazole") ~ "Nitroimidazoles",
      str_detect(tolower(`Drug Class`), "mupirocin") ~ "Mupirocin-like",
      str_detect(tolower(`Drug Class`), "pleuromutilin") ~ "Pleuromutilins",
      str_detect(tolower(`Drug Class`), "bicyclomycin") ~ "Bicyclomycin-like",
      str_detect(tolower(`Drug Class`), "nitrofuran") ~ "Nitrofurans",
      str_detect(tolower(`Drug Class`), "fusidane") ~ "Fusidane-like",
      str_detect(tolower(`Drug Class`), "diarylquinoline") ~ "Diarylquinolines",
      str_detect(tolower(`Drug Class`), "isoniazid") ~ "Isoniazid-like",
      str_detect(tolower(`Drug Class`), "cycloserine") ~ "Cycloserine-like",
      str_detect(tolower(`Drug Class`), "pyrazine") ~ "Pyrazine-like",
      str_detect(tolower(`Drug Class`), "thioamide") ~ "Thioamide-like",
      str_detect(tolower(`Drug Class`), "salicylic acid") ~ "Salicylic Acid-like",
      str_detect(tolower(`Drug Class`), "tetracenomycin") ~ "Tetracenomycin-like",
      
      # Disinfectants / antiseptics
      str_detect(tolower(`Drug Class`), "disinfect|antiseptic") ~ "Disinfecting/Antiseptics",
      str_detect(tolower(`Drug Class`), "antibacterial free fatty acid") ~ "Fatty Acid Antibiotics",
      
      TRUE ~ `Drug Class`  # keep original if it doesn't match
    )
  )


drug_class_list = card_merged_w_ont_df %>%  dplyr::select(card_edited_drug_class ) %>% distinct() %>% pull()
  
categorized_drug_class = c("Betalactams", "Sulfonamides", "Trimethoprim", "Tetracyclines", "Aminocoumarin", "Aminoglycoside", "Elfamycin", "Fosfomycin", "Glycopeptides", "MLS", "Peptides", "Phenicol", "Quinolones", "Rifamycines")  
    
```
#NMDS
##NMDS on Resistome (CARD derived relative abundances)-ONLY SHARED DATES WITH BIOFILM
```{r}
#specify dates
dates = c("2024-09-09", "2024-09-11", "2024-09-12", "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09")

analysis_df <- card_merged_df %>%
  filter(
    Type %in% c("Drain", "Endcap", "H_WW")
  ) %>%
  mutate(
    Date   = as.Date(Date),
    Type   = factor(Type, levels = c("H_WW","Endcap","Drain" ))
  ) %>%
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>%
  drop_na(timepoint)


biofilm_dates <- analysis_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
  summarise(
    biofilm_ref_date = Date,
    .groups = "drop"
  ) %>% distinct()

filtered_df <- analysis_df %>%
  mutate(Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  ) %>% 
  
  dplyr::select(-best_date,-date_diff_days,-biofilm_ref_date)

metadata = filtered_df %>%
  
  #filter dates
  dplyr::select(SampleID,Location,Type,Month,Date) %>% distinct() %>% 
  
  
  column_to_rownames(var = "SampleID")

#filter for one outcome variable (using aro_accession because there are multiple aro_accessions per gene_name)
data <- filtered_df %>% 
  
  #  #filter dates
  # filter(
  #   (Type == "H_WW" & Date %in% dates) |
  #   (Type != "H_WW")
  # ) %>% 
  
  dplyr::select(SampleID,gene_id, RPKM) %>% distinct() %>% 
  
  group_by(SampleID,gene_id) %>% 
  mutate(RPKM = mean(RPKM,na.rm=T)) %>% 
  ungroup() %>% 
  
  distinct() 


# create dataframe with rows as samples names and columns as ARGs
counts <- data %>% tidyr::pivot_wider(names_from = gene_id, values_from = RPKM) %>% 
  
  #replace NAs with zeros
  mutate_all( ~replace(., lengths(.)==0, 0))

#make row names sample_id
counts <- counts %>% column_to_rownames(var = "SampleID")

#make sure it's a numeric dataframe
counts = as.data.frame(counts) %>% mutate(across(everything(), ~ as.numeric(as.character(.))))


#replace NAs with 0
counts[is.na(counts)] <- 0

#remove samples 
#counts <- counts[!(rownames(counts) %in% c("S209", "S210","S1","S2","S3")), ]


dist2 <- vegdist(counts)
argMDS <- isoMDS(dist2, trace = F)
scores <- as.data.frame(argMDS)
scores <- scores[,1:2]


# Assuming 'scores' and 'metadata' are your dataframes
common_rows <- intersect(rownames(counts), rownames(metadata))

# Filter 'scores' dataframe to keep only rows present in 'metadata'
filtered_counts <- counts[rownames(counts) %in% common_rows, ]

# Filter 'scores' dataframe to keep only rows present in 'metadata'
filtered_metadata <- metadata[rownames(metadata) %in% common_rows, ]

# now plot the NMDS plot with envfit variables
vare.mds <- metaMDS(filtered_counts, trace = FALSE, distance = "jaccard")

ef <- envfit(vare.mds, filtered_metadata, permu = 999,na.rm = TRUE)

#plote results
plot(vare.mds, display = "sites")


#extract envfit arrows data and filter out non significant one 
#en_coord_cont = as.data.frame(scores(ef, "vectors")) * ordiArrowMul(ef) 
#en_coord_cont$pval <- ef[["vectors"]][["pvals"]]
#en_coord_cont <- filter(en_coord_cont,pval<=0.05)


#extract nmds scores and associate to metadata 
nmds_scores <- as.data.frame(scores(vare.mds)$sites)

# Assuming your_data is your dataframe
nmds_scores$SampleID <- rownames(nmds_scores)
rownames(nmds_scores) <- NULL # now get rid of row names

full_metadata = metadata %>% rownames_to_column(var = "SampleID") %>% 
  
  #create new variable combined with site and sample type
  mutate("site_type"=paste0(Location,"_",Type))

nmdsMerged = left_join(nmds_scores,full_metadata, by = "SampleID") 


#specify Locations we want to incude in plot
locations = c("FIONA","SHREK","OSCAR","MARIO","LUIGI")

nmds_plotg <- nmdsMerged %>% 
  
  filter(Location %in% locations) %>% 
  
  # now replace Site names 
  mutate(Location = case_when(Location == "MARIO" ~ "Site 1",
                              Location == "FIONA" ~ "Site 2",
                              Location == "LUIGI" ~ "Site 3",
                              Location == "SHREK" ~ "Site 4",
                              Location == "OSCAR" ~ "Site 5",
                              .default = Location)) %>% 
  
  #filter out municipal ww
  filter(Type != "Mu_WW") %>% 
  
  #rename some variables
  mutate(Type = case_when(Type == "Mu_WW"~ "Municipal WW",
                          Type == "H_WW" ~ "Wastewater",
                          .default = Type)) %>% 
  
  ggplot(aes(x = NMDS1, y = NMDS2, color = Type, shape = Location, group = Type)) +
  geom_point() +
  stat_ellipse(linewidth = 0.75) +
  
  #annotate("text",x=NMDS.mean$MDS1,y=NMDS.mean$MDS2,label=NMDS.mean$group)+
  
  #fill colors
  scale_color_manual(values = c(resistome_colors)) +
  theme_bw() +
  labs(
    title = "NMDS for the Resistome",
    x = "NMDS1",
    y = "NMDS2",
    color = "",
    size = ""
  ) + 
  #geom_vline(xintercept = 0,linetype = "dashed")+ geom_hline(yintercept = 0,linetype = "dashed")+
  #geom_segment(aes(x = 0, y = 0, xend =NMDS1, yend = NMDS2), data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  #geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", fontface = "bold", label = row.names(en_coord_cont)) +

  
  theme(
    axis.text.y = element_blank(),
    legend.text = element_text(size = 13,color = "black"))

# print plot
nmds_plotg


#ggsave 
ggsave(plot = nmds_plotg, path = "Biofilm Project Figures", "CARD Based NMDS of Resistome-ONLY SHARED DATES WITH BIOFILM.png", dpi = 300, height = 6, width = 7, units = "in")
    
    


#prettier NMDS plot with shapes and colors
nmds_plotg <- nmdsMerged %>% 
  
  #filter out municipal ww
  filter(Type != "Mu_WW") %>% 
  
  #rename some variables
  mutate(Type = case_when(Type == "Mu_WW"~ "Municipal WW",
                          Type == "Endcap" ~ "Branch Biofilm",
                          Type == "Drain" ~ "Sink Biofilm",
                          Type == "H_WW" ~ "Wastewater",
                          .default = Type)) %>% 
  
  # now replace Site names 
  mutate(Location = case_when(Location == "MARIO" ~ "Site 1",
                              Location == "FIONA" ~ "Site 2",
                              Location == "LUIGI" ~ "Site 3",
                              Location == "SHREK" ~ "Site 4",
                              Location == "OSCAR" ~ "Site 5",
                              .default = Location)) %>% 
  
  #filter locations if needed
  #filter(Location %in% locations) %>% 
  ggplot(aes(x = NMDS1, y = NMDS2, color = Type, shape = Location, group = Type)) +
  geom_point() +
  stat_ellipse(linewidth = 0.75) +
  
  #change colors
  #scale_color_manual(values = c(
    #"4" = "lightsalmon", "12" = "grey",
   #                             )) +
  #scale_shape_manual(values = c(1, 2
                                #, 3, 4, 5, 6, 16, 17, 18, 19, 20, 21,22, 29
                                #)) +
  theme_bw() +
  labs(
    title = "NMDS Plot",
    x = "NMDS1",
    y = "NMDS2",
    color = "site_type",
    size = "Collection Date"
  ) + 
  
  #manualy change colors
  scale_color_manual(values = c("Sink Biofilm" = "#79307D", "Branch Biofilm" = "#417C8C", 
                                "Municipal WW" = resistome_colors[[3]], 
                                "Wastewater" = "#E57262"
                                )) +
  #scale_shape_manual(values = c("Drain" = 1, "Endcap" = 2, "Municipal WW" = 3,  "Hospital WW" = 4 )) +
  #geom_vline(xintercept = 0,linetype = "dashed")+ geom_hline(yintercept = 0,linetype = "dashed")+
  #geom_segment(aes(x = 0, y = 0, xend =NMDS1, yend = NMDS2), data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  #geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", fontface = "bold", label = row.names(en_coord_cont)) +
  labs(title = "",
       color = "",
       shape = "") +
  theme(axis.text.y = element_text(angle = 45, hjust = 1)) + 
  theme(axis.text.y = element_blank(),
        axis.title = element_text(size= 20,face = "bold",color = "black"),
        legend.text = element_text(size = 15,face = "bold", color = "black"))

# print plot
nmds_plotg

# Save the combined plot as a PNG file
ggsave(plot = nmds_plotg, path = "Biofilm Project Figures", "CARD Resistome derived NMDS (species level) by Site and Location-ONLY SHARED DATES WITH BIOFILM.png", dpi = 300, height = 5, width = 7.5, units = "in")

ggsave(plot = nmds_plotg, path = "Biofilm Project Figures", "CARD Resistome derived NMDS (species level) by Site and Location-ONLY SHARED DATES WITH BIOFILM.svg", height = 6, width = 7)



```


#PERMANOVA
##Permanova on resistome (CARD)
```{r}
#specify Locations we want to incude in plot
locations = c("FIONA","SHREK","OSCAR","MARIO","LUIGI")

# Find the right row indices
rows <- which(filtered_metadata[, 1] %in% locations)

# Subset the matrix based on the selected rows
subset_metadata <- filtered_metadata[rows, ]

# Get row names that are present in both data frames
common_rows <- intersect(rownames(subset_metadata), rownames(filtered_counts))

subset_adonis = filtered_counts[common_rows,]

#run pairwise adonis using function
pairwise_adonis_resistome = pairwise.adonis2(subset_adonis, subset_metadata$Location, sim.method = 'bray', p.adjust.m = 'bonferroni')

adonis2(subset_adonis ~ Location + Type , subset_metadata, permutations = 9999, by = "margin")

#adonis2 by Location
#adonis2(subset_adonis ~ Location, subset_metadata, permutations = 9999)


#run pairwise adonis using function
pairwise_adonis_microbiome = pairwise.adonis2(subset_adonis, subset_metadata$Type, sim.method = 'bray', p.adjust.m = 'bonferroni')

#adonis.pair(vegdist(subset_adonis),subset_metadata$subarea)

```

#Diversity
##Diversity Metrics for Resistome (CARD)-ONLY SHARED DATES WITH BIOFILM
```{r}
#specify dates
dates = c("2024-09-09", "2024-09-11", "2024-09-12", "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09")

analysis_df <- card_merged_df %>%
  filter(
    Type %in% c("Drain", "Endcap", "H_WW")
  ) %>%
  mutate(
    Date   = as.Date(Date),
    Type   = factor(Type, levels = c("H_WW","Endcap","Drain" ))
  ) %>%
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>%
  drop_na(timepoint)


biofilm_dates <- analysis_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
  mutate(
    biofilm_ref_date = Date
  ) %>% 
    ungroup() %>% 
  dplyr::distinct(Location,timepoint,biofilm_ref_date)

filtered_df <- analysis_df %>%
  mutate(Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  ) %>% 
  
  dplyr::select(-best_date,-date_diff_days,-biofilm_ref_date)

metadata = metadata_df %>% 
  column_to_rownames(var = "SampleID")


library(data.table)

# Instead of dplyr mutate_all
transposed <- filtered_df %>% 
  
  # #filter dates
  # filter(
  #   (Type == "H_WW" & Date %in% dates) |
  #   (Type != "H_WW")
  # ) %>% 
  
  dplyr::select(SampleID, gene_id, RPKM) %>% # average RPKM per sample-gene pair
  distinct() %>% # pivot long → wide 
  group_by(SampleID, gene_id) %>% # group in case of multiple entries
  summarise(rpkm = mean(RPKM, na.rm = TRUE), .groups = "drop") %>% # values in cells = RPKM
  pivot_wider(names_from = SampleID, values_from = rpkm, values_fill = 0) # replace missing with zero

# Use data.table for conversion
transposed_dt <- as.data.table(transposed) # turn into a data.table

# Convert to matrix with genes as rows, samples as cols
gene_matrix <- as.matrix(transposed_dt[, -"gene_id"])
rownames(gene_matrix) <- transposed_dt$gene_id

# Transpose so rows = samples, cols = genes
counts <- t(gene_matrix)

# Make sure it's integer
storage.mode(counts) <- "integer"

# Now metadata rownames should match
all(rownames(counts) %in% rownames(metadata))

# Assuming 'scores' and 'metadata' are your dataframes
common_rows <- intersect(rownames(counts), rownames(metadata))

# Filter 'scores' dataframe to keep only rows present in 'metadata'
filtered_counts <- counts[rownames(counts) %in% common_rows, ]


# Filter 'scores' dataframe to keep only rows present in 'metadata'
filtered_metadata <- metadata[rownames(metadata) %in% common_rows, ]

# Start clean from filtered_counts
wide <- as.data.frame(filtered_counts)

# # Force everything to numeric explicitly
# wide[] <- lapply(wide, function(x) as.numeric(as.character(x)))
# 
# # Round and coerce to integer
# wide[] <- lapply(wide, function(x) as.integer(round(x, 0)))

# Convert to matrix
wide <- as.matrix(wide)

# Replace NAs with 0
wide[is.na(wide)] <- 0

# Drop any rows that are all zeros
wide <- wide[rowSums(wide) > 0, ]

# # Double-check: must return "integer"
# str(wide)

#S <- specnumber(wide)
# raremax <- max(rowSums(wide))

#Srare <- as.data.frame(vegan::rrarefy(wide, raremax))

#shannon <- as.data.frame(vegan::diversity(as.matrix(Srare)))
#shannon <- as.data.frame(vegan::diversity(Srare))
# simpson <- as.data.frame(vegan::diversity(Srare, "simpson"))
# richness <- as.data.frame(specnumber(Srare)) 
# evenness <- as.data.frame(shannon/log(richness))

#shannon  <- vegan::diversity(wide, index = "shannon")
shannon = as.data.frame(vegan::diversity(wide,"shannon"))
simpson  <- as.data.frame(vegan::diversity(wide, index = "simpson")) 
richness <- as.data.frame(specnumber(wide))          # counts non-zero features
evenness <- as.data.frame(shannon / log(richness))

#spread <- as.data.frame(t(Srare))

#chao1 <- as.data.frame(apply(spread, 2, chao1))

diversity <- bind_cols(shannon, simpson, richness, evenness) %>% rownames_to_column() %>% distinct()

names(diversity)[1:5] <- c("SampleID", "Shannon", "Simpson", "Richness", "Evenness")

#specify Locations we want to incude in plot
locations = c("FIONA","SHREK","OSCAR","MARIO","LUIGI")

diversity_cleaned = diversity %>%   
  
  
  left_join(metadata_df, by = "SampleID") %>% 
  
  #filter subareas
  filter(Location %in% locations) %>% 
  
  #select needed variables
  dplyr::select(SampleID, Shannon,Location, Type) %>% 
    
    drop_na(Shannon) %>% 
  
  #clean
  distinct() 

#diversity$Process <- factor(diversity$Process, levels = c("PE", "AS", "SE", "FE"))

# Calculate the median values for each group
medians <- aggregate(Shannon ~ Location, diversity_cleaned, median)

# Reorder the levels of the 'Group' factor based on the median values
diversity_cleaned$Location <- factor(diversity_cleaned$Location, levels = medians$Location[order(medians$Shannon)])

plot = ggbetweenstats(
  data = diversity_cleaned,
  x    = Location, 
  y    = Shannon, 
  outlier.tagging = T,
  type = "nonparametric", 
  p.adjust.method = "bonferroni",
  conf.level = 0.95,
  #pairwise.display = "none",
  pairwise.display = "significant",
  mean.plotting = TRUE,
  mean.ci = FALSE,
  centrality.plotting = FALSE,
  #point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.4, size = 3, stroke = 0),
  #median.point.args = list(size = 0.2, color = "pink"),
  mean.label.args = list(size = 0.5)
) + 
  ylab("Shannon Diversity for Resistome")+
  xlab("")+
  labs(color = "")+
  guides(color = "none")+
  theme_classic()+
  theme(axis.text.x= element_text(size = 9,color= "black"))
  #theme(axis.text.x= element_text(angle = 45,size = 7,vjust = 0.9,hjust = 0.9,color= "black"))

print(plot)

# Save the combined plot as a PNG file
ggsave(plot = plot, path = "Biofilm Project Figures", "CARD Based Resistome Diversity by Location-ONLY FOR SHARED DATES WITH BIOFILM.png", dpi = 300, height = 6, width = 7, units = "in")




tidied_diversity = diversity_cleaned %>%
  
  # now replace Site names 
  mutate(Location = case_when(Location == "MARIO" ~ "Site 1",
                              Location == "FIONA" ~ "Site 2",
                              Location == "LUIGI" ~ "Site 3",
                              Location == "SHREK" ~ "Site 4",
                              Location == "OSCAR" ~ "Site 5",
                              .default = Location)) %>% 
  
  filter(Type != "Mu_WW") %>% 
  
  mutate(Type = case_when(Type =="H_WW"~ "Wastewater",Type == "Mu_WW" ~ "Municipal_WW",
                        Type == "Endcap"~ "Sewer Biofilm",
                        Type == "Drain" ~ "Sink Biofilm",
                          .default = Type)) 

#new ggstats plot wiht diversity by sample type
# Calculate the median values for each group
medians <- aggregate(Shannon ~ Type, tidied_diversity, median)

# Reorder the levels of the 'Group' factor based on the median values
tidied_diversity$Type <- factor(tidied_diversity$Type, levels = medians$Type[order(medians$Shannon)])

plot = ggbetweenstats(
  data = tidied_diversity,
  x    = Type, 
  y    = Shannon, 
  outlier.tagging = T,
  type = "nonparametric", 
  p.adjust.method = "bonferroni",
  conf.level = 0.95,
  #pairwise.display = "none",
  mean.plotting = TRUE,
  mean.ci = FALSE,
  #point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.4, size = 3, stroke = 0),
  #median.point.args = list(size = 0.2, color = "pink"),
  mean.label.args = list(size = 0.5)
) + 
  ylab("Shannon Diversity for Resistome")+
  xlab("")+
  labs(color = "")+
  guides(color = "none")+
  theme_classic()+
  theme(axis.text.x= element_text(size = 9,face= "bold",color= "black"))
  #theme(axis.text.x= element_text(angle = 45,size = 7,vjust = 0.9,hjust = 0.9,color= "black"))

print(plot)

ggsave(plot = plot, path = "Biofilm Project Figures", "CARD Based Resistome Diversity by Sample Type-ONLY FOR SHARED DATES WITH BIOFILM.png", dpi = 300, height = 6, width = 7, units = "in")

```





#Venn Diagrams
## Venn diagram for sites shared-ONLY SHARED DATES WITH BIOFILM-updated
```{r}

#specify dates
dates = c("2024-09-09", "2024-09-11", "2024-09-12", "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09")

# ---- Prepare genus-level matrix ----
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")



analysis_df <- card_merged_df %>%
  filter(
    Location %in% locations,
    Type %in% c("Drain", "Endcap", "H_WW"),
     (Type == "H_WW" & Date %in% dates) | (Type != "H_WW")
  ) %>%
  mutate(
    rel_ab = as.numeric(rel_abundance),
    Date   = as.Date(Date),
    Type   = factor(Type, levels = c("Drain", "Endcap", "H_WW"))
  ) %>%
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>%
  drop_na(timepoint)


biofilm_dates <- analysis_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
   mutate(
    biofilm_ref_date = Date
  ) %>% 
    ungroup() %>% 
  dplyr::distinct(Location,timepoint,biofilm_ref_date)

filtered_df <- analysis_df %>%
  mutate(Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  ) %>% 
  
  dplyr::select(-best_date,-date_diff_days,-biofilm_ref_date)


library(dplyr)
library(ggvenn)

# Colors
earthy_green  <- "#9FB498"
earthy_purple <- "#8488AC"
earthy_orange <- "#E57262"

viridisCol1 = "#79307D"
viridisCol2 = "#417C8C"
viridsCol3 = "#E57262"

# Subareas to include
locations <- c("FIONA","SHREK","OSCAR","MARIO","LUIGI")

# Species-level rows only
venn_df <- filtered_df 

# Loop through each location
for (loc in locations) {

  df <- venn_df %>%
    ungroup() %>%
    distinct(Type, Location, Date, gene_id, rel_ab) %>%
    filter(Location == loc) %>%
    filter(rel_ab != 0) %>%
    drop_na(Location)

  # Build Venn input
  ggvenn_df <- list(
    Sink_Biofilm   = df %>% filter(Type == "Drain") %>% pull(gene_id),
    Sewer_Biofilm = df %>% filter(Type == "Endcap") %>% pull(gene_id),
    Hospital_WW  = df %>% filter(Type == "H_WW") %>% pull(gene_id)
  )

  # Plot with a title for clarity
  print(
    ggvenn(
      ggvenn_df,
      fill_color = c(viridisCol1, viridisCol2, viridsCol3),
      stroke_size = 0.5,
      set_name_size = 4
    ) +
    ggtitle(paste("Species overlap at", loc))
  )
}


#now average overlap across Sites

library(dplyr)
library(purrr)
library(ggvenn)

# Function to compute overlaps at one location
get_overlap_counts <- function(df) {
  sets <- list(
    Sink_Biofilm   = df %>% filter(Type == "Drain") %>% pull(gene_id),
    Sewer_Biofilm = df %>% filter(Type == "Endcap") %>% pull(gene_id),
    Hospital_WW  = df %>% filter(Type == "H_WW") %>% pull(gene_id)
  )
  
  # All unique species
  all_species <- unique(unlist(sets))
  
  # Build matrix: species × sets
  mat <- sapply(sets, function(s) all_species %in% s)
  
  # Count combinations (like Venn regions)
  counts <- apply(mat, 1, function(x) paste(which(x), collapse = "_"))
  tibble(region = counts) %>%
    dplyr::count(region, gene_id = "n")
}

# ---- STEP 1: Compute overlaps per location ----
overlap_by_loc <- venn_df %>%
  ungroup() %>%
  distinct(Type, Location, Date, gene_id, rel_ab) %>%
  filter(rel_ab != 0, Location %in% locations) %>%
  split(.$Location) %>%
  map(~ get_overlap_counts(.x))

# ---- STEP 2: Average counts across locations ----
avg_overlap <- bind_rows(
  lapply(names(overlap_by_loc), function(loc) {
    overlap_by_loc[[loc]] %>%
      mutate(Location = loc)
  })
) %>%
  group_by(region) %>%
  summarise(avg_count = mean(n), .groups = "drop")

print(avg_overlap)

library(eulerr)

fit <- euler(c(
  "Sink_Biofilm"             = avg_overlap$avg_count[avg_overlap$region == "1"],
  "Sewer_Biofilm"           = avg_overlap$avg_count[avg_overlap$region == "2"],
  "Hospital_WW"            = avg_overlap$avg_count[avg_overlap$region == "3"],
  "Sink_Biofilm&Sewer_Biofilm"= avg_overlap$avg_count[avg_overlap$region == "1_2"],
  "Sink_Biofilm&Hospital_WW" = avg_overlap$avg_count[avg_overlap$region == "1_3"],
  "Sewer_Biofilm&Hospital_WW"= avg_overlap$avg_count[avg_overlap$region == "2_3"],
  "Sink_Biofilm&Sewer_Biofilm&Hospital_WW"= avg_overlap$avg_count[avg_overlap$region == "1_2_3"]
))

plot(fit, fills = c(viridisCol1, viridisCol2, viridsCol3), edges = TRUE)

#counts
plot(
  fit,
  fills = c(viridisCol1, viridisCol2, viridsCol3),
  edges = TRUE,
  labels = list(font = 2, cex = 0.8),   # bold, smaller text
  quantities = list(type = "counts", font = 2, cex = 0.8) # numbers inside
)



#percentages
plot(
  fit,
  fills = c(viridisCol1, viridisCol2, viridsCol3),
  edges = TRUE,
  labels = list(font = 2, cex = 0.8),   # bold, smaller text
  quantities = list(type = "percent", font = 2, cex = 0.8) # numbers inside
)



library(gridGraphics)
library(grid)

# Draw the base eulerr plot onto the grid device
grid.echo(function() {
  plot(
    fit,
    fills = c(viridisCol1, viridisCol2, viridsCol3),
    edges = TRUE,
    labels = list(font = 2, cex = 0.8),
    quantities = list(type = "counts", font = 2, cex = 0.8)
  )
})

# Capture whatever is currently on the device
g <- grid.grab()

# Save with ggsave
ggsave("Biofilm Project Figures/Kraken CT 0.15 Based Euler diagram Average Overlap -Counts -ONLY SHARED DATES WITH BIOFILM-RPIP only.svg", plot = g, width = 6, height = 6)

# Save with ggsave
ggsave("Biofilm Project Figures/Kraken CT 0.15 Based Euler diagram Average Overlap -Counts -ONLY SHARED DATES WITH BIOFILM-RPIP only.png", plot = g, width = 6, height = 6,dpi = 300)

# # Draw the base eulerr plot onto the grid device
# grid.echo(function() {
#   plot(
#     fit,
#     fills = c(viridisCol1, viridisCol2, viridsCol3),
#     edges = TRUE,
#     labels = list(font = 2, cex = 0.8),
#     quantities = list(type = "percent", font = 2, cex = 0.8)
#   )
# })
# 
# # Capture whatever is currently on the device
# g <- grid.grab()
# 
# # Save with ggsave
# ggsave("Biofilm Project Figures/Kraken CT 0.15 Based Euler diagram Average Overlap -Percentages -ONLY SHARED DATES WITH BIOFILM-RPIP only.svg", plot = g, width = 6, height = 6)
# # Save with ggsave
# ggsave("Biofilm Project Figures/Kraken CT 0.15 Based Euler diagram Average Overlap -Percentages -ONLY SHARED DATES WITH BIOFILM-RPIP only.png", plot = g, width = 6, height = 6,dpi = 300)


#overlap for total pathogens across all Sites------------------------


# Species-level rows only
venn_df <- filtered_df 



  df <- venn_df %>%
    ungroup() %>%
    distinct(Type, SampleID, gene_id, rel_ab) %>%
    filter(rel_ab != 0) %>% 
    mutate(presence = 1) %>% dplyr::select(-rel_ab) %>% distinct()
  
  df <- filtered_df %>%
  distinct(Type, SampleID, gene_id, rel_ab) %>%
  filter(rel_ab != 0) %>%
  dplyr::select(Type, gene_id) %>%
  distinct()
  
  sets_total <- list(
  Sink_Biofilm   = df %>% filter(Type == "Drain") %>% pull(gene_id),
  Sewer_Biofilm = df %>% filter(Type == "Endcap") %>% pull(gene_id),
  Hospital_WW   = df %>% filter(Type == "H_WW") %>% pull(gene_id)
)

all_species <- unique(unlist(sets_total))

presence_mat <- sapply(sets_total, function(s) all_species %in% s)

region_df <- tibble(
  region = apply(presence_mat, 1, function(x)
    paste(which(x), collapse = "_"))
) %>%
  dplyr::count(region, gene_id = "n")

  
print(region_df)


library(eulerr)

fit_total <- euler(c(
  "Sink_Biofilm"   = region_df$n[region_df$region == "1"],
  "Sewer_Biofilm" = region_df$n[region_df$region == "2"],
  "Hospital_WW"   = region_df$n[region_df$region == "3"],

  "Sink_Biofilm&Sewer_Biofilm" =
    region_df$n[region_df$region == "1_2"],

  "Sink_Biofilm&Hospital_WW" =
    region_df$n[region_df$region == "1_3"],

  "Sewer_Biofilm&Hospital_WW" =
    region_df$n[region_df$region == "2_3"],

  "Sink_Biofilm&Sewer_Biofilm&Hospital_WW" =
    region_df$n[region_df$region == "1_2_3"]
))

plot(
  fit_total,
  fills = c(viridisCol1, viridisCol2, viridsCol3),
  edges = TRUE,
  labels = list(font = 2, cex = 0.9),
  quantities = list(type = "counts", font = 2, cex = 0.9)
)

# Capture whatever is currently on the device
g <- grid.grab()

# Save with ggsave
ggsave("Biofilm Project Figures/Kraken CT 0.15 Based Euler diagram Total -Counts -ONLY SHARED DATES WITH BIOFILM-RPIP only.svg", plot = g, width = 6, height = 6)

# Save with ggsave
ggsave("Biofilm Project Figures/Kraken CT 0.15 Based Euler diagram Total -Counts -ONLY SHARED DATES WITH BIOFILM-RPIP only.png", plot = g, width = 6, height = 6,dpi = 300)

#percentages
plot(
  fit_total,
  fills = c(viridisCol1, viridisCol2, viridsCol3),
  edges = TRUE,
  labels = list(font = 2, cex = 0.9),
  quantities = list(type = "percent", font = 2, cex = 0.9)
)

# Capture whatever is currently on the device
g <- grid.grab()

# Save with ggsave
ggsave("Biofilm Project Figures/Kraken CT 0.15 Based Euler diagram Total -Percents -ONLY SHARED DATES WITH BIOFILM-RPIP only.svg", plot = g, width = 6, height = 6)

# Save with ggsave
ggsave("Biofilm Project Figures/Kraken CT 0.15 Based Euler diagram Total -Percents -ONLY SHARED DATES WITH BIOFILM-RPIP only-percents.png", plot = g, width = 6, height = 6,dpi = 300)

  
  # Build Venn input
  ggvenn_df <- list(
    Sink_Biofilm   = df %>% filter(Type == "Drain") %>% pull(gene_id),
    Sewer_Biofilm = df %>% filter(Type == "Endcap") %>% pull(gene_id),
    Hospital_WW  = df %>% filter(Type == "H_WW") %>% pull(gene_id)
  )

  # Plot with a title for clarity
  print(
    ggvenn(
      ggvenn_df,
      fill_color = c(viridisCol1, viridisCol2, viridsCol3),
      stroke_size = 0.5,
      set_name_size = 4
    ) +
    ggtitle("Species overlap")
  )

```


#UpSet Plots

##upset plots without facets  (ggplot and UpSet version) (was this gene ever seen in this environment?)-ONLY SHARED DATES WITH BIOFILM
```{r}

#colors--------

resistome_colors = c("#8488AC", "#9FB498", "#BA6646", "#565656", "#865856", "#B99796", "#F8E5DC")

posterCol1=  "#8488AC"
posterCol2 = "#9FB498"
posterCol3 = "#BA6646"
posterCol4 = "#565656"
posterCol5 = "#865856"
posterCol6 = "#B99796"
posterCol7 = "#F8E5DC"

viridisCol1 = "#79307D"
viridisCol2 = "#417C8C"
viridsCol3 = "#E57262"

colors = c( "#417C8C","#417C8C","#417C8C","#417C8C","#E57262","#79307D","#79307D","#79307D")

#colors = c( "#79307D","#79307D","#417C8C","#417C8C","#417C8C","#417C8C","#E57262","#79307D","#79307D","#79307D")

colors = c( "#8488AC","#BA6646","#9FB498")


poster_colors = c("#dc9298","#a51c2f","#3d3d3d","#f3e3e5","#417c8c","#e57262","#915493")


poster_colors = c("#dc9298","#a51c2f","#3d3d3d","#f3e3e5","#417c8c","#e57262","#915493")

poster_colors = c("#ce2e2d","#c00300","#88b8c5","#dd9298","#ab74a8","#fc3b3b","#a61d30","#eededf","#dc6c64","#467d91")


# Specify locations to include------
#specify dates
dates = c("2024-09-09", "2024-09-11", "2024-09-12", "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09")

locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame


# Specify AMR genes------
amr_genes <- c("blaTEM", "blaSHV", "blaCTX-M", "blaOXA", "blaKPC", "blaNDM", "blaVIM", "blaIMP",
               "mecA", "vanA", "vanB", "ermB", "tetM", "tetO", "aac(3)-I", "aph(3')-IIIa",
               "strA", "strB","sul1", "sul2", "dfrA", "dfrB", "qnrA", "qnrB", "qnrS")





# ---------- Extract and clean ----------
analysis_df <- card_merged_w_ont_df %>%
  filter(
    Type %in% c("Drain", "Endcap", "H_WW"),
    Location %in% locations
  ) %>%
  mutate(
    Date   = as.Date(Date),
    rel_ab = as.numeric(rel_abundance),
    Type   = factor(Type, levels = c("H_WW","Endcap","Drain" ))
  ) %>%
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>%
  drop_na(timepoint)


biofilm_dates <- analysis_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
  summarise(
    biofilm_ref_date = Date,
    .groups = "drop"
  ) %>% distinct()

filtered_df <- analysis_df %>%
  mutate(rel_ab = as.numeric(rel_abundance),
         Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  ) %>% 
  
  dplyr::select(-best_date,-date_diff_days,-biofilm_ref_date)


df <- filtered_df %>% 
  
  #clean ontology df
  dplyr::select(Gene,SampleID, Date, Location, Type,rel_ab,`Drug Class`, `Resistance Mechanism`, class, mechanism,card_edited_drug_class, group) %>% 
  
  
  #filter args in rpip list
  #filter(Species %in% c(amr_genes)) %>% 
  
  # #filter dates
  # filter(
  #   (Type == "H_WW" & Date %in% dates) |
  #   (Type != "H_WW")
  # ) %>% 
  
  mutate(
    Month = format(as.Date(Date), "%m-%Y")
  ) 



# Plot with gene as category for upset------
plot_df <- df %>%

  dplyr::select(Gene,SampleID, Date, Location, Type,rel_ab) %>% 
  
  
  #make presence/absence
  mutate(rel_ab = ifelse(rel_ab>0,1,0)) %>%
  mutate(
    Type  = case_when(Type =="H_WW"   ~ "Wastewater",
                      Type =="Endcap" ~ "Sewer Biofilm",
                   Type =="Drain"  ~ "Sink Biofilm",
                   .default = Type)
  )
  
  


wide_df <- plot_df %>%
  group_by(Gene, Type) %>%
  
  #filter for at least one detect
  summarise(detect = as.integer(sum(rel_ab) > 0), .groups = "drop") %>%
  pivot_wider(
    names_from = Type,
    values_from = detect,
    values_fill = 0
  )



library(ggupset)
library(forcats)

n_genes = 45


wide_long <- wide_df %>%
  pivot_longer(cols = c("Sink Biofilm","Sewer Biofilm","Wastewater"),
               names_to = "Type",
               values_to = "detect") %>%
  filter(detect == 1) %>%
  group_by(Gene) %>%
  summarise(sets = list(unique(Type)), .groups = "drop") %>%
  ungroup() 


# Calculate top genes globally
genes_freq <- wide_long %>%
  dplyr::count(Gene, sort = TRUE)

keep_genes <- genes_freq %>%
  slice_head(n = n_genes) %>%
  pull(Gene)

# Relabel the rest as "Other"
wide_long_cut <- wide_long %>%
  mutate(
    Genes_lumped = ifelse(Gene %in% keep_genes, Gene, "Other"),
    Genes_lumped = factor(Genes_lumped)
  )

# Verify the levels
length(levels(wide_long_cut$Genes_lumped))
levels(wide_long_cut$Genes_lumped)

genes_levels = levels(wide_long_cut$Genes_lumped)


# #viridis palette
# 
# library(viridisLite)
# 
# genes_levels <- levels(wide_long$Genes_lumped)
# 
# palette <- c(
#   viridis(length(genes_levels)-1, option = "D"),
#   "grey70"
# )
# names(palette) <- c(setdiff(genes_levels, "Other"), "Other")


# Assign colors: first N get a palette, "Other" gets grey
palette <- c(fortyfive_pal, "grey70")
names(palette) <- c(setdiff(genes_levels, "Other"), "Other")


# Plot 

plot = ggplot(wide_long_cut, aes(x = sets, fill = Genes_lumped)) +
  geom_bar() +
  scale_x_upset() +
  labs(
    fill = "Genes",
    x = "",
    y = "Number of Genes\nShared Between Groups"
  ) +
  scale_fill_manual(values = palette) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_rect(colour = "black"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_markdown(size = 7.5),  # <- enable markdown
    axis.title.y = element_text(margin = margin(r = 5)) # tighten y label
  ) +
  theme_combmatrix(
    combmatrix.panel.point.color.fill = "black",
    combmatrix.label.make_space = TRUE,   # <- stop reserving extra space
    combmatrix.label.text = element_text(color = "black", size = 8.5),
    #combmatrix.label.text = element_blank(),
    combmatrix.panel.line.color = "black",
    combmatrix.panel.striped_background = FALSE,
    combmatrix.panel.point.size = 2
  )

print(plot)

ggsave(path = "Biofilm Project Figures",file="Unfaceted CARD UpSet plot of ARGS detected-ONLY SHARED DATES WITH BIOFILM.png", plot=plot, width=7.5, height=4,dpi = 500)
ggsave(path = "Biofilm Project Figures",file="Unfaceted CARD UpSet plot of ARGS detected-ONLY SHARED DATES WITH BIOFILM.svg", plot=plot, width=7.5, height=4)




#plot at class level ---------------

  

plot_df <- df %>%

  dplyr::select(Gene,SampleID, Date, Location, Type,rel_ab,`Drug Class`) %>% 
  
  drop_na(Type) %>% 
  
  
  #make presence/absence
  mutate(rel_ab = ifelse(rel_ab>0,1,0)) %>%
  mutate(
    Type  = case_when(Type =="H_WW"   ~ "Wastewater",
                      Type =="Endcap" ~ "Sewer Biofilm",
                   Type =="Drain"  ~ "Sink Biofilm",
                   .default = Type)
  )
  
  


wide_df <- plot_df %>%
  group_by(Gene, `Drug Class` ,Type) %>%
  
  #filter for at least one detect
  summarise(detect = as.integer(sum(rel_ab) > 0), .groups = "drop") %>%
  pivot_wider(
    names_from = Type,
    values_from = detect,
    values_fill = 0
  )



library(ggupset)
library(forcats)

n_classes = 45


wide_long <- wide_df %>%
  pivot_longer(cols = c("Sink Biofilm","Sewer Biofilm","Wastewater"),
               names_to = "Type",
               values_to = "detect") %>%
  filter(detect == 1) %>%
  group_by(Gene, `Drug Class`) %>%
  summarise(sets = list(unique(Type)), .groups = "drop") %>%
  ungroup() 


# Calculate top genes globally
classes_freq <- wide_long %>%
  dplyr::count(`Drug Class`, sort = TRUE)

keep_classes <- classes_freq %>%
  slice_head(n = n_classes) %>%
  pull(`Drug Class`)

# Relabel the rest as "Other"
wide_long_cut <- wide_long %>%
  mutate(
    classes_lumped = ifelse(`Drug Class` %in% keep_classes, `Drug Class`, "Other"),
    classes_lumped = factor(classes_lumped)
  )

# Verify the levels
length(levels(wide_long_cut$classes_lumped))


class_levels = levels(wide_long_cut$classes_lumped)




# #viridis palette
# 
# library(viridisLite)
# 
# genes_levels <- levels(wide_long$Genes_lumped)
# 
# palette <- c(
#   viridis(length(genes_levels)-1, option = "D"),
#   "grey70"
# )
# names(palette) <- c(setdiff(genes_levels, "Other"), "Other")


# Assign colors: first N get a palette, "Other" gets grey
palette <- c(fortyfive_pal, "grey70")
names(palette) <- c(setdiff(class_levels, "Other"), "Other")


# Plot 

plot = ggplot(wide_long_cut, aes(x = sets, fill = classes_lumped)) +
  geom_bar() +
  scale_x_upset() +

  labs(
    fill = "CARD Drug Class",
    x = "",
    y = "Number of Genes\nShared Between Groups"
  ) +
  scale_fill_manual(values = palette) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_rect(colour = "black"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_markdown(size = 7.5),  # <- enable markdown
    axis.title.y = element_text(margin = margin(r = 5)) # tighten y label
  ) +
  theme_combmatrix(
    combmatrix.panel.point.color.fill = "black",
    combmatrix.label.make_space = TRUE,   # <- stop reserving extra space
    combmatrix.label.text = element_text(color = "black", size = 9),
    #combmatrix.label.text = element_blank(),
    combmatrix.panel.line.color = "black",
    combmatrix.panel.striped_background = FALSE,
    combmatrix.panel.point.size = 2
  )

print(plot)

ggsave(path = "Biofilm Project Figures",file="Unfaceted CARD UpSet plot of ARGS detected by CARD Drug Class-ONLY SHARED DATES WITH BIOFILM.png", plot=plot, width=7.5, height=8,dpi = 500)
#ggsave(path = "Biofilm Project Figures",file="Unfaceted CARD UpSet plot of ARGS detected by CARD Drug Class-ONLY SHARED DATES WITH BIOFILM.svg", plot=plot, width=7.5, height=4)



#plot at arg mechanism level ---------------


plot_df <- df %>%

  dplyr::select(Gene,SampleID, Date, Location, Type,rel_ab,`Resistance Mechanism`) %>% 
  
  drop_na(Type) %>% 
  
    mutate(`Resistance Mechanism` = ifelse(grepl(";", `Resistance Mechanism`), "Multi", `Resistance Mechanism`)) %>%
  
  
  #make presence/absence
  mutate(rel_ab = ifelse(rel_ab>0,1,0)) %>%
  mutate(
    Type  = case_when(Type =="H_WW"   ~ "Wastewater",
                      Type =="Endcap" ~ "Sewer Biofilm",
                   Type =="Drain"  ~ "Sink Biofilm",
                   .default = Type)
  )
  
  


wide_df <- plot_df %>%
  group_by(Gene, `Resistance Mechanism` ,Type) %>%
  
  #filter for at least one detect
  summarise(detect = as.integer(sum(rel_ab) > 0), .groups = "drop") %>%
  pivot_wider(
    names_from = Type,
    values_from = detect,
    values_fill = 0
  )



library(ggupset)
library(forcats)

n_classes = 45

wide_long <- wide_df %>%
  pivot_longer(cols = c("Sink Biofilm","Sewer Biofilm","Wastewater"),
               names_to = "Type",
               values_to = "detect") %>%
  filter(detect == 1) %>% 
  group_by(Gene, `Resistance Mechanism`) %>%
  mutate(sets = list(unique(Type))) %>% 
  ungroup()

wide_long_clean <- wide_df %>%
  pivot_longer(cols = c("Sink Biofilm","Sewer Biofilm","Wastewater"),
               names_to = "Type",
               values_to = "detect") %>%
  filter(detect == 1) %>% 
  group_by(Gene, `Resistance Mechanism`) %>%
  mutate(sets = list(unique(Type))) %>% 
  ungroup() %>% 
  
  dplyr::select(Gene, sets, `Resistance Mechanism`) %>%
  distinct() 


# Calculate top genes globally
classes_freq <- wide_long_clean %>%
  dplyr::count(`Resistance Mechanism`, sort = TRUE)

keep_classes <- classes_freq %>%
  slice_head(n = n_classes) %>%
  pull(`Resistance Mechanism`)

# Relabel the rest as "Other"
wide_long_cut <- wide_long_clean %>%
  mutate(
    classes_lumped = ifelse(`Resistance Mechanism` %in% keep_classes,`Resistance Mechanism`, "Other"),
    classes_lumped = factor(classes_lumped)
  )

# Verify the levels
length(levels(wide_long_cut$classes_lumped))


class_levels = levels(wide_long_cut$classes_lumped)




# #viridis palette
# 
# library(viridisLite)
# 
# genes_levels <- levels(wide_long$Genes_lumped)
# 
# palette <- c(
#   viridis(length(genes_levels)-1, option = "D"),
#   "grey70"
# )
# names(palette) <- c(setdiff(genes_levels, "Other"), "Other")


# Assign colors: first N get a palette, "Other" gets grey
palette <- c(fortyfive_pal, "grey70")
names(palette) <- c(setdiff(class_levels, "Other"), "Other")


# Plot 

plot = ggplot(wide_long_cut, aes(x = sets, fill = classes_lumped)) +
  geom_bar() +
  scale_x_upset() +

  labs(
    fill = "CARD Res Mechanism",
    x = "",
    y = "Number of Genes\nShared Between Groups"
  ) +
  scale_fill_manual(values = palette) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_rect(colour = "black"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_markdown(size = 7.5),  # <- enable markdown
    axis.title.y = element_text(margin = margin(r = 5)) # tighten y label
  ) +
  theme_combmatrix(
    combmatrix.panel.point.color.fill = "black",
    combmatrix.label.make_space = TRUE,   # <- stop reserving extra space
    combmatrix.label.text = element_text(color = "black", size = 9),
    #combmatrix.label.text = element_blank(),
    combmatrix.panel.line.color = "black",
    combmatrix.panel.striped_background = FALSE,
    combmatrix.panel.point.size = 2
  )

print(plot)

ggsave(path = "Biofilm Project Figures",file="Unfaceted CARD UpSet plot of ARGS detected by CARD Res Mechanism-ONLY SHARED DATES WITH BIOFILM.png", plot=plot, width=9, height=4,dpi = 500)
#ggsave(path = "Biofilm Project Figures",file="Unfaceted CARD UpSet plot of ARGS detected by CARD Res Mechanism-ONLY SHARED DATES WITH BIOFILM.svg", plot=plot, width=7.5, height=4)



#plot at MEGARes Drug Class level ---------------



plot_df <- df %>%
  dplyr::select(Gene, SampleID, Date, Location, Type, rel_ab, `Drug Class`, class) %>%
  drop_na(Type) %>%
  
  mutate(
    class = if_else(
      is.na(class) | class == "",
      # Take text before first colon OR semicolon
      str_trim(str_extract(`Drug Class`, "^[^:;]+")),
      class
    )
  ) %>%
  
  dplyr::select(-`Drug Class`) %>%
  mutate(
    rel_ab = ifelse(rel_ab > 0, 1, 0),
    Type = case_when(
      Type == "H_WW" ~ "Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain" ~ "Sink Biofilm",
      .default = Type
    )
  )


  


wide_df <- plot_df %>%
  group_by(Gene,  class ,Type) %>%
  
  #filter for at least one detect
  summarise(detect = as.integer(sum(rel_ab) > 0), .groups = "drop") %>%

  pivot_wider(
    names_from = Type,
    values_from = detect,
    values_fill = 0
  )



library(ggupset)
library(forcats)

n_classes = 45

wide_long <- wide_df %>%
  pivot_longer(cols = c("Sink Biofilm","Sewer Biofilm","Wastewater"),
               names_to = "Type",
               values_to = "detect") %>%
  filter(detect == 1) %>% 
  group_by(Gene, class) %>%
  mutate(sets = list(unique(Type))) %>% 
  ungroup()

wide_long_clean <- wide_df %>%
  pivot_longer(cols = c("Sink Biofilm","Sewer Biofilm","Wastewater"),
               names_to = "Type",
               values_to = "detect") %>%
  filter(detect == 1) %>% 
  group_by(Gene, class) %>%
  mutate(sets = list(unique(Type))) %>% 
  ungroup() %>% 
  
  dplyr::select(Gene, sets, class) %>%
  distinct() 


# Calculate top genes globally
classes_freq <- wide_long_clean %>%
  dplyr::count(class, sort = TRUE)

keep_classes <- classes_freq %>%
  slice_head(n = n_classes) %>%
  pull(class)

# Relabel the rest as "Other"
wide_long_cut <- wide_long_clean %>%
  mutate(
    classes_lumped = ifelse(class %in% keep_classes,class, "Other"),
    classes_lumped = factor(classes_lumped)
  )

# Verify the levels
length(levels(wide_long_cut$classes_lumped))


class_levels = levels(wide_long_cut$classes_lumped)




# #viridis palette
# 
# library(viridisLite)
# 
# genes_levels <- levels(wide_long$Genes_lumped)
# 
# palette <- c(
#   viridis(length(genes_levels)-1, option = "D"),
#   "grey70"
# )
# names(palette) <- c(setdiff(genes_levels, "Other"), "Other")


# Assign colors: first N get a palette, "Other" gets grey
palette <- c(fortyfive_pal, "grey70")
names(palette) <- c(setdiff(class_levels, "Other"), "Other")


# Plot 

plot = ggplot(wide_long_cut, aes(x = sets, fill = classes_lumped)) +
  geom_bar() +
  scale_x_upset() +

  labs(
    fill = "MEGARes Drug Class",
    x = "",
    y = "Number of Genes\nShared Between Groups"
  ) +
  scale_fill_manual(values = palette) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_rect(colour = "black"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_markdown(size = 7.5),  # <- enable markdown
    axis.title.y = element_text(margin = margin(r = 5)) # tighten y label
  ) +
  theme_combmatrix(
    combmatrix.panel.point.color.fill = "black",
    combmatrix.label.make_space = TRUE,   # <- stop reserving extra space
    combmatrix.label.text = element_text(color = "black", size = 8),
    #combmatrix.label.text = element_blank(),
    combmatrix.panel.line.color = "black",
    combmatrix.panel.striped_background = FALSE,
    combmatrix.panel.point.size = 2
  )

print(plot)

ggsave(path = "Biofilm Project Figures",file="Unfaceted CARD UpSet plot of ARGS detected by MEGARes Drug Class plus CARD-ONLY SHARED DATES WITH BIOFILM.png", plot=plot, width=8.5, height=4,dpi = 500)
ggsave(path = "Biofilm Project Figures",file="Unfaceted CARD UpSet plot of ARGS detected by MEGARes Drug Class plus CARD-ONLY SHARED DATES WITH BIOFILM.svg", plot=plot, width=8.5, height=4)





#plot at Drug Class level drop CARD annotations that don't have MEGARES matches ---------------




plot_df <- df %>%
  dplyr::select(Gene, SampleID, Date, Location, Type, rel_ab, `Drug Class`, class) %>%
  drop_na(Type) %>%
  
  drop_na(class) %>% 
  
  dplyr::select(-`Drug Class`) %>%
  
  mutate(
    rel_ab = ifelse(rel_ab > 0, 1, 0),
    Type = case_when(
      Type == "H_WW" ~ "Hospital Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain" ~ "Sink Biofilm",
      .default = Type
    )
  )


  


wide_df <- plot_df %>%
  group_by(Gene,  class ,Type) %>%
  
  #filter for at least one detect
  summarise(detect = as.integer(sum(rel_ab) > 0), .groups = "drop") %>%
  pivot_wider(
    names_from = Type,
    values_from = detect,
    values_fill = 0
  )



library(ggupset)
library(forcats)

n_classes = 45

wide_long <- wide_df %>%
  pivot_longer(cols = c("Sink Biofilm","Sewer Biofilm","Hospital Wastewater"),
               names_to = "Type",
               values_to = "detect") %>%
  filter(detect == 1) %>% 
  group_by(Gene, class) %>%
  mutate(sets = list(unique(Type))) %>% 
  ungroup()

wide_long_clean <- wide_df %>%
  pivot_longer(cols = c("Sink Biofilm","Sewer Biofilm","Hospital Wastewater"),
               names_to = "Type",
               values_to = "detect") %>%
  filter(detect == 1) %>% 
  group_by(Gene, class) %>%
  mutate(sets = list(unique(Type))) %>% 
  ungroup() %>% 
  
  dplyr::select(Gene, sets, class) %>%
  distinct() 


# Calculate top genes globally
classes_freq <- wide_long_clean %>%
  dplyr::count(class, sort = TRUE)

keep_classes <- classes_freq %>%
  slice_head(n = n_classes) %>%
  pull(class)

# Relabel the rest as "Other"
wide_long_cut <- wide_long_clean %>%
  mutate(
    classes_lumped = ifelse(class %in% keep_classes,class, "Other"),
    classes_lumped = factor(classes_lumped)
  )

# Verify the levels
length(levels(wide_long_cut$classes_lumped))


class_levels = levels(wide_long_cut$classes_lumped)




# #viridis palette
# 
# library(viridisLite)
# 
# genes_levels <- levels(wide_long$Genes_lumped)
# 
# palette <- c(
#   viridis(length(genes_levels)-1, option = "D"),
#   "grey70"
# )
# names(palette) <- c(setdiff(genes_levels, "Other"), "Other")


# Assign colors: first N get a palette, "Other" gets grey
palette <- c(fortyfive_pal, "grey70")
names(palette) <- c(setdiff(class_levels, "Other"), "Other")


# Plot 

plot = ggplot(wide_long_cut, aes(x = sets, fill = classes_lumped)) +
  geom_bar() +
  scale_x_upset() +

  labs(
    fill = "MEGARes Drug Class",
    x = "",
    y = "Number of ARGs\nShared Between Groups"
  ) +
  scale_fill_manual(values = palette) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_rect(colour = "black"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_markdown(size = 7.5),  # <- enable markdown
    axis.title.y = element_text(margin = margin(r = 5)) # tighten y label
  ) +
  theme_combmatrix(
    combmatrix.panel.point.color.fill = "black",
    combmatrix.label.make_space = TRUE,   # <- stop reserving extra space
    combmatrix.label.text = element_text(color = "black", size = 8),
    #combmatrix.label.text = element_blank(),
    combmatrix.panel.line.color = "black",
    combmatrix.panel.striped_background = FALSE,
    combmatrix.panel.point.size = 2
  )

print(plot)

ggsave(path = "Biofilm Project Figures",file="Unfaceted CARD UpSet plot of ARGS detected by MEGARes Drug Class-ONLY SHARED DATES WITH BIOFILM.png", plot=plot, width=7.5, height=4,dpi = 500)
ggsave(path = "Biofilm Project Figures",file="Unfaceted CARD UpSet plot of ARGS detected by MEGARes Drug Class-ONLY SHARED DATES WITH BIOFILM.svg", plot=plot, width=7.5, height=4)




#plot at Drug Class level drop CARD annotations that don't have MEGARES matches ---------------




plot_df <- df %>%
  dplyr::select(Gene, SampleID, Date, Location, Type, rel_ab, card_edited_drug_class) %>%
  drop_na(Type) %>%
  
  
  mutate(
    rel_ab = ifelse(rel_ab > 0, 1, 0),
    Type = case_when(
      Type == "H_WW" ~ "Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain" ~ "Sink Biofilm",
      .default = Type
    )
  )


  


wide_df <- plot_df %>%
  group_by(Gene,card_edited_drug_class ,Type) %>%
  
  #filter for at least one detect
  summarise(detect = as.integer(sum(rel_ab) > 0), .groups = "drop") %>%
  pivot_wider(
    names_from = Type,
    values_from = detect,
    values_fill = 0
  )



library(ggupset)
library(forcats)

n_classes = 45



wide_long_clean <- wide_df %>%
  pivot_longer(cols = c("Sink Biofilm","Sewer Biofilm","Wastewater"),
               names_to = "Type",
               values_to = "detect") %>%
  filter(detect == 1) %>% 
  group_by(Gene, card_edited_drug_class) %>%
  mutate(sets = list(unique(Type))) %>% 
  ungroup() %>% 
  
  dplyr::select(Gene, sets, card_edited_drug_class) %>%
  distinct() 


# Calculate top genes globally
classes_freq <- wide_long_clean %>%
  dplyr::count(card_edited_drug_class, sort = TRUE)

keep_classes <- classes_freq %>%
  slice_head(n = n_classes) %>%
  pull(card_edited_drug_class)

categorized_drug_class = c("Betalactams", "Sulfonamides", "Trimethoprim", "Tetracyclines", "Aminocoumarin", "Aminoglycoside", "Elfamycin", "Fosfomycin", "Glycopeptides", "MLS", "Peptides", "Phenicol", "Quinolones", "Rifamycines")  

keep_classes  = categorized_drug_class

# Relabel the rest as "Other"
wide_long_cut <- wide_long_clean %>%
  mutate(
    classes_lumped = ifelse(card_edited_drug_class %in% keep_classes,card_edited_drug_class, "Other"),
    classes_lumped = factor(classes_lumped)
  )

# Verify the levels
length(levels(wide_long_cut$classes_lumped))


class_levels = levels(wide_long_cut$classes_lumped)




# #viridis palette
# 
# library(viridisLite)
# 
# genes_levels <- levels(wide_long$Genes_lumped)
# 
# palette <- c(
#   viridis(length(genes_levels)-1, option = "D"),
#   "grey70"
# )
# names(palette) <- c(setdiff(genes_levels, "Other"), "Other")


# Assign colors: first N get a palette, "Other" gets grey
palette <- c(fortyfive_pal, "grey70")
names(palette) <- c(setdiff(class_levels, "Other"), "Other")


# Plot 

plot = ggplot(wide_long_cut, aes(x = sets, fill = classes_lumped)) +
  geom_bar() +
  scale_x_upset() +

  labs(
    fill = "Drug Class",
    x = "",
    y = "Number of ARGs\nShared Between Groups"
  ) +
  scale_fill_manual(values = palette) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_rect(colour = "black"),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_markdown(size = 7.5),  # <- enable markdown
    axis.title.y = element_text(margin = margin(r = 5)) # tighten y label
  ) +
  theme_combmatrix(
    combmatrix.panel.point.color.fill = "black",
    combmatrix.label.make_space = TRUE,   # <- stop reserving extra space
    combmatrix.label.text = element_text(color = "black", size = 8),
    #combmatrix.label.text = element_blank(),
    combmatrix.panel.line.color = "black",
    combmatrix.panel.striped_background = FALSE,
    combmatrix.panel.point.size = 2
  )

print(plot)

ggsave(path = "Biofilm Project Figures",file="Unfaceted CARD UpSet plot of ARGS detected by Pruden Lab Drug Class-ONLY SHARED DATES WITH BIOFILM.png", plot=plot, width=4.5, height=4,dpi = 500)
ggsave(path = "Biofilm Project Figures",file="Unfaceted CARD UpSet plot of ARGS detected by Pruden Lab Drug Class-ONLY SHARED DATES WITH BIOFILM.svg", plot=plot, width=4.5, height=4)

#---------------Get counts of ARGs in each category for Pruden lab categories---------
# Create a summary table of the intersections
intersection_counts <- wide_long_clean %>%
  mutate(set_string = sapply(sets, function(x) paste(sort(x), collapse = " & "))) %>%
  group_by(set_string) %>%
  summarise(n_genes = n(), .groups = "drop") %>%
  arrange(desc(n_genes))

# View the results
print(intersection_counts)

```

#UpSet Plot- Genes Shared by Relative Abundance or Presence Absence
##upset plot set themes
```{r}
font_add_google("Reem Kufi", "Reem Kufi")   # <-- REQUIRED
showtext_auto()                             # <-- activates showtext


## ggplot theme
theme_set(theme_minimal(base_family = "Reem Kufi", base_size = 12))

theme_update(
  plot.title = element_text(size = 27,
                            face = "bold",
                            hjust = .5,
                            margin = margin(10, 0, 30, 0)),
  plot.caption = element_text(size = 9,
                              color = "grey40",
                              hjust = .5,
                              margin = margin(20, 0, 5, 0)),
  axis.text.y = element_blank(),
  axis.title = element_blank(),
  plot.background = element_rect(fill = "grey88", color = NA),
  panel.background = element_rect(fill = NA, color = NA),
  panel.grid = element_blank(),
  panel.spacing.y = unit(0, "lines"),
  strip.text.y = element_blank(),
  legend.position = "bottom",
  legend.text = element_text(size = 9, color = "grey40"),
  legend.box.margin = margin(t = 30), 
  legend.background = element_rect(color = "grey40", 
                                   size = .3, 
                                   fill = "grey95"),
  legend.key.height = unit(.25, "lines"),
  legend.key.width = unit(2.5, "lines"),
  plot.margin = margin(rep(20, 4))
)

#second-------------------

font_add_google("Chivo", "Chivo")
font_add_google("Passion One", "Passion One")

## ggplot theme
theme_set(theme_minimal(base_family = "Chivo"))

theme_update(axis.text.x = element_text(size = 11, color = "grey20"),
             axis.text.y = element_text(size = 13, color = "black", face = "bold"),
             axis.ticks.x = element_line(color = "grey45"),
             axis.ticks.y = element_blank(),
             axis.ticks.length.x = unit(.4, "lines"),
             panel.grid = element_blank(),
             plot.background = element_rect(fill = "grey60", color = "grey60"))
```

##Upset Plot - ARGs (genes of concern) by Presence Absence- All Sites in One-ONLY SHARED DATES WITH BIOFILM-filtered for what's enriched (listed approach)
```{r eval=FALSE, include=FALSE}
#theme set----------
#font_add_google("Chivo", "Chivo")
#font_add_google("Passion One", "Passion One")

## ggplot theme
theme_set(theme_minimal(base_family = "Chivo"))

theme_update(axis.text.x = element_text(size = 11, color = "grey20"),
             axis.text.y = element_text(size = 13, color = "black"),
             axis.ticks.x = element_line(color = "grey45"),
             axis.ticks.y = element_blank(),
             axis.ticks.length.x = unit(.4, "lines"),
             panel.grid = element_blank(),
             plot.background = element_rect(fill = "grey60", color = "grey60"))

#font_add_google("Reem Kufi", "Reem Kufi")   # <-- REQUIRED
showtext_auto()                             # <-- activates showtext


## ggplot theme
theme_set(theme_minimal(base_family = "Reem Kufi", base_size = 12))

theme_update(
  plot.title = element_text(size = 27,
                            face = "bold",
                            hjust = .5,
                            margin = margin(10, 0, 30, 0)),
  plot.caption = element_text(size = 9,
                              color = "grey40",
                              hjust = .5,
                              margin = margin(20, 0, 5, 0)),
  axis.text.y = element_blank(),
  axis.title = element_blank(),
  plot.background = element_rect(fill = "grey88", color = NA),
  panel.background = element_rect(fill = NA, color = NA),
  panel.grid = element_blank(),
  panel.spacing.y = unit(0, "lines"),
  strip.text.y = element_blank(),
  legend.position = "bottom",
  legend.text = element_text(size = 9, color = "grey40"),
  legend.box.margin = margin(t = 30), 
  legend.background = element_rect(color = "grey40", 
                                   size = .3, 
                                   fill = "grey95"),
  legend.key.height = unit(.25, "lines"),
  legend.key.width = unit(2.5, "lines"),
  plot.margin = margin(rep(20, 4))
)

#R prep---------------------------------------
## packages
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggtext)
library(showtext)
library(sysfonts)
library(ggplot2)
library(ggstatsplot)
library(dplyr)
library(tidyverse)
library(readxl)
library(vegan)
library(OTUtable)
library(vegan)
library(MASS)
library(ggtext)
library(ggupset)
library(forcats)

#font_add_google("Chivo", "Chivo")
#font_add_google("Passion One", "Passion One")

## ggplot theme
theme_set(theme_minimal(base_family = "Chivo"))

theme_update(axis.text.x = element_text(size = 11, color = "grey10"),
             axis.text.y = element_text(size = 13, color = "black"),
             axis.ticks.x = element_line(color = "white"),
             axis.ticks.y = element_blank(),
             axis.ticks.length.x = unit(.4, "lines"),
             title = element_text(size = 55),
             panel.grid = element_blank(),
             plot.background = element_rect(fill = "white", color = "white"))

#R data -------------------------
dates = c("2024-09-09", "2024-09-11", "2024-09-12", "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09")

# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame

#genes

#blactx-m, mecA, tetA, tetW,ermB, vanA, qnrS, blatem

#filter for contains CTX-M-1 in ARO Name
#filter for tet(A) as ARO Name
# filter for contains TEM- in ARO Name
# filter for mecA in ARO Name (colistin gene)
# filter for tet(W) as ARO Name
# filter for ErmB as ARO Name
#filter for Gen = vanA or ARO is exactly vanA
#filter for Gen = sul1 or ARO is exactly sul1
#filter for has qnrS in ARO Name
#filter for Gen = vanB or ARO is exactly vanB
#filter for Gen = sul2 or ARO is exactly sul2
#filter for Gen = sul3 or ARO is exactly sul3
#filter for contains NDM in ARO Name
#filter for contains KPC in ARO Name
#filter for contains SHV in ARO Name
#filter for contains MCR in ARO Name
#filter for contains cfrA in ARO Name
#filter for contains QnrA in ARO Name
#filter for contains tet(M) in ARO Name
#filter for qac genes

# from  Diamantina Lymperatou 2025 paper
#primary hosts of ARGs within HWS were found to be Escherichia coli and Klebsiella pneumoniae

#most common ARGs in wastewater:
#betalactams: blaCTX-M,bla-SHV, bla-TEM, bla-OXA-ESBLs
#Carbapenems: bla-KPC, bla-NDM, bla-VIM, bla-OXA-48, bla-IMP 
#Colistin: mcr-1 to mcr-10, pmrA, pmrB
#Tetracyclines: tetA, tetB, tetM, tetO, tetX, tetC, tetD, tetG
#Sulfonamides: sul1, sul2, sul3
#vancomycin: vanA, vanB, vanC, vanD, vanE, vanG, vanL, vanM, vanN

#most of concer:
#mcr-5.1

#MRSA related: mecA, FemB, and SCCmec elements

#from Xinyi Shuai 2024 paper
#highest priority: 
#Betalactams: blaOXA-181, blaOXA-212, bla-OXA-309, blaOXA-1, blaMOX-7, blaGES-1, bla OXA 13, blaVIM-2, blaIMP-1, blaBEL-2, blaOXA-236, blaFOX-7, blaPER-3, blaOXA-5, blaNPS-1, bla FOX-4, bla KPC-2, bla OXA-36, bla TEM-1, blaOXA-101, blaOXA-333, blaFOX-6, bla MOX-6, blaOXA-119, blaTEM-83, blaGES-5, blaOXA-17, blaOXA-58, blaOXA-24, blaNDM-1, blaOXA-4, blaCTX-M-61, bla OXA-211, bla OXA-10, bla CTX-M, bla VEB-3, blaIMP-25, blaOXA-2,blaCMY-111
#MLS: ermB, ermF, msrC, ereB, ereA, mefA
#multidrug: emrE, adeK, emrB, mexB, adeJ, EmrB-QacA, mdtE, TolC, mdtK, adeC, mexT, mexI
# Polymyxin: mcr-3, mcr-5
#quinolone: qnrS, qnrB

#high priority
#vancomycin: vanH, vanS, vanB, vanX
#aminoglycoside: APH(6)-Id, ant(3'')-lh-aac(6'-lld), ant(2'')-I, aadB, aph(3'')-I, aph(6)-I, aac(6')-30-aac(6')-lb', aac(3)-II, aadA
#fosfomycin: fosA
#fosmidomycin: rosB
#Glycylcyclines: tetX2

#highly important:
#Chloramphenicol: cmlA
#sulfonamide: sul1, sul2
#tetracyclines: tetA, tetQ, tetE, tetO, tetM, tetW, tet39
#trimethoprim: dfrA17, dfrA23,dfrA12
#bacitracin: bacA




#filter for contains OXA in ARO Name
library(stringr)
library(stringr)
library(dplyr)



analysis_df <- card_merged_w_ont_df %>%
  filter(
    Type %in% c("Drain", "Endcap", "H_WW")
  ) %>%
  mutate(
    Date   = as.Date(Date),
    rel_ab = as.numeric(rel_abundance),
    Type   = factor(Type, levels = c("H_WW","Endcap","Drain" ))
  ) %>%
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>%
  drop_na(timepoint)


biofilm_dates <- analysis_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
  mutate(
    biofilm_ref_date = Date
  ) %>% 
  ungroup() %>% 
  dplyr::distinct(Location,timepoint,biofilm_ref_date)

filtered_df <- analysis_df %>%
  mutate(rel_ab = as.numeric(rel_abundance),
         Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  ) %>% 
  
  dplyr::select(-best_date,-date_diff_days,-biofilm_ref_date)


# 1. Clean your target list (ensure it's a simple character vector)
target_genes <- tolower(as.character(rpip_alleles))

library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(tibble)



df <- filtered_df %>%

  mutate(
    TypeLabel = case_when(
      Type == "H_WW"   ~ "Wastewater",
      Type == "Drain"  ~ "Sink Biofilm",
      Type == "Endcap" ~ "Sewer Biofilm",
      TRUE ~ as.character(Type)
    ),
    Month = format(Date, "%m-%Y")
  ) %>%
  distinct(Gene, SampleID, Date, Location, Type, card_edited_drug_class, rel_ab, TypeLabel, Month, `Model Name`,`ARO Name`) %>% 

  # 2. Filter using the 'Model Name' column 
  # (This matches the nomenclature in your list best)
  #filter(tolower(`Model Name`) %in% target_genes) %>% 



#filter args in candidate list
filter(
  # contains CTX-M-1
  str_detect(`ARO Name`, regex("CTX-M-", ignore_case = TRUE)) |
  # exactly tet(A)
  str_detect(`ARO Name`, regex("tet\\(A\\)", ignore_case = TRUE)) |
  # exactly tet(M)
  str_detect(`ARO Name`, regex("tet\\(M\\)", ignore_case = TRUE)) |
  # contains TEM- (for example TEM-1, TEM-2, etc)
  str_detect(`ARO Name`, regex("TEM-", ignore_case = TRUE)) |
  # contains NDM-
  str_detect(`ARO Name`, regex("NDM-", ignore_case = TRUE)) |
  # contains KPC-
  str_detect(`ARO Name`, regex("KPC-", ignore_case = TRUE)) |
  # contains SHV-
  str_detect(`ARO Name`, regex("SHV-", ignore_case = TRUE)) |
  # contains OXA-
  str_detect(`ARO Name`, regex("OXA-", ignore_case = TRUE)) |
  # contains MCR-
  str_detect(`ARO Name`, regex("MCR-", ignore_case = TRUE)) |
  # exactly mecA
  str_to_lower(`ARO Name`) == "meca" |
  # exactly tet(W)
  str_detect(`ARO Name`, regex("tet\\(W\\)", ignore_case = TRUE)) |
  # exactly tet(M)
  str_detect(`ARO Name`, regex("tet\\(M\\)", ignore_case = TRUE)) |
  # exactly ErmB (or case-insensitive)
  str_to_lower(`ARO Name`) == "ermb" |
  # exactly QnrA (or case-insensitive)
  str_to_lower(`ARO Name`) == "qnra" |
  # Gene == vanA OR ARO Name exactly vanA
  Gene == "vanA" |
  str_to_lower(`ARO Name`) == "vana" |  # if your ARO Name column is just "vanA"
  # Gene == cfrA OR ARO Name exactly cfrA
  Gene == "cfrA" |
  str_to_lower(`ARO Name`) == "cfra" |  # if your ARO Name column is just "vanA"
    # Gene == vanB OR ARO Name exactly vanB
  Gene == "vanA" |
  str_to_lower(`ARO Name`) == "vanb" |  # if your ARO Name column is just "vanB"
  # Gene == sul1 OR ARO Name exactly sul1
  Gene == "sul1" |
  str_to_lower(`ARO Name`) == "sul1" |
  # Gene == sul2 OR ARO Name exactly sul2
  Gene == "sul2" |
  str_to_lower(`ARO Name`) == "sul2" |
  # Gene == sul3 OR ARO Name exactly sul3
  Gene == "sul3" |
  str_to_lower(`ARO Name`) == "sul3" |
  # contains qnrS anywhere
  str_detect(`ARO Name`, regex("qnrS", ignore_case = TRUE)) |

    Gene == "MCR-5.1" |
  #filter for qac genes
  str_detect(`ARO Name`, regex("qac", ignore_case = TRUE))
) %>%

  mutate(Gene = case_when(str_detect(Gene, regex("CTX-M-1\\b", ignore_case = TRUE)) ~ "bla-CTX-M-1",
                        str_detect(Gene, regex("CTX-M-", ignore_case = TRUE)) ~ "bla-CTX-M",
                        str_detect(Gene, regex("tet\\(A\\)", ignore_case = TRUE))~ "tetA",
                        str_detect(Gene, regex("tet\\(M\\)", ignore_case = TRUE))~ "tetM",
                        str_detect(Gene, regex("TEM-", ignore_case = TRUE)) ~ "blaTEM",
                        str_detect(Gene, regex("NDM-", ignore_case = TRUE)) ~ "blaNDM",
                        str_detect(Gene, regex("KPC-", ignore_case = TRUE)) ~ "blaKPC",
                        str_detect(Gene, regex("SHV-", ignore_case = TRUE)) ~ "blaSHV",
                        str_detect(Gene, regex("OXA-", ignore_case = TRUE)) ~ "blaOXA",
                        str_detect(Gene, regex("MCR-5.1", ignore_case = TRUE)) ~ "mcr-5.1",
                        str_detect(Gene, regex("MCR-4.3", ignore_case = TRUE)) ~ "mcr-4.3",
                        str_detect(Gene, regex("MCR-1", ignore_case = TRUE)) ~ "mcr-1",
                        str_detect(Gene, regex("MCR-2", ignore_case = TRUE)) ~ "mcr-2",
                        str_detect(Gene, regex("MCR-3", ignore_case = TRUE)) ~ "mcr-3",
                        str_detect(Gene, regex("MCR-6", ignore_case = TRUE)) ~ "mcr-6",
                        str_detect(Gene, regex("MCR-7", ignore_case = TRUE)) ~ "mcr-7",
                        str_detect(Gene, regex("MCR-8", ignore_case = TRUE)) ~ "mcr-8",
                        str_detect(Gene, regex("MCR-9", ignore_case = TRUE)) ~ "mcr-9",
                        str_detect(Gene, regex("MCR-", ignore_case = TRUE)) ~ "mcr",
                        str_detect(Gene, regex("tet\\(W\\)", ignore_case = TRUE)) ~ "tetW",
                        str_detect(Gene, regex("meca", ignore_case = TRUE)) ~ "mecA",
                        str_detect(Gene, regex("ermb", ignore_case = TRUE)) ~ "ermB",
                        str_detect(Gene, regex("vana", ignore_case = TRUE)) ~ "vanA",
                        str_detect(Gene, regex("vanb", ignore_case = TRUE)) ~ "vanB",
                        str_detect(Gene, regex("cfra", ignore_case = TRUE)) ~ "cfrA",
                        str_detect(Gene, regex("qnra", ignore_case = TRUE)) ~ "QnrA",
                        str_detect(Gene, regex("sul1", ignore_case = TRUE)) ~ "sul1",
                        str_detect(Gene, regex("sul2", ignore_case = TRUE)) ~ "sul2",
                        str_detect(Gene, regex("sul3", ignore_case = TRUE)) ~ "sul3",
                        str_detect(Gene, regex("qnrs", ignore_case = TRUE)) ~ "qnrS",
                        str_detect(Gene, regex("qac", ignore_case = TRUE)) ~ "qac",
                        .default = Gene)) %>%



  
  #clean ontology df
  dplyr::select(Gene,SampleID, Date, Location, Type,rel_ab,`Model Name`,
                #`Drug Class`, `Resistance Mechanism`, class, mechanism, group,
                card_edited_drug_class) %>% 
  
  
  # #filter dates
  # filter(
  #   (Type == "H_WW" & Date %in% dates) |
  #   (Type != "H_WW")
  # ) %>% 
  
  mutate(
    TypeLabel = case_when(Type == "H_WW" ~ "Wastewater", 
                          Type == "Drain" ~ "Sink Biofilm", 
                          Type == "Endcap" ~ "Sewer Biofilm", 
                          TRUE ~ as.character(Type)),
    Month = format(as.Date(Date), "%m-%Y")
  ) %>% 
  
  #distinct genes per sample id, date, sample type, etc.
  distinct(Gene,SampleID,Date,Location,Type,card_edited_drug_class,rel_ab,TypeLabel,Month,`Model Name`)
  
  # 3. Check if you missed anything
print(paste("Original genes in list:", length(unique(target_genes))))
print(paste("Genes found in dataframe:", length(unique(df$`Model Name`))))
         
#Create duplicate rows so mcr 5.1 is included inthe broader mcr count
# Define your overlaps here
overlaps <- tribble(
  ~Specific,   ~Broad,
  "mcr-5.1",   "mcr",
  "mcr-4.3",   "mcr",
  "mcr-1",   "mcr",
  "mcr-2",   "mcr",
  "mcr-3",   "mcr",
  "mcr-4",   "mcr",
  "mcr-5",   "mcr",
  "mcr-6",   "mcr",
  "mcr-7",   "mcr",
  "mcr-8",   "mcr",
  "mcr-9",   "mcr",
  "mcr-10",   "mcr",
  "bla-CTX-M-1",   "bla-CTX-M"
)

# Create the duplicated rows based on the overlap table
extra_rows <- df %>%
  inner_join(overlaps, by = c("Gene" = "Specific")) %>%
  mutate(Gene = Broad) %>%  # Re-label the specific gene as the Broad name
  dplyr::select(-Broad)     # Remove the helper column

# Combine back with your original data
df_expanded <- bind_rows(df, extra_rows)
  
  
# Count detects per gene × Type
df_detects <- df_expanded %>%
  mutate(detect = rel_ab > 0) %>% 
  
  #distinct genes per sample id, date, sample type, etc.
  distinct(Gene,SampleID,Date,Location,Type,card_edited_drug_class,detect,TypeLabel,Month) %>% 
  
  group_by(Gene, Type) %>%
  summarise(
    n_detects_type = sum(detect, na.rm = TRUE),   # how many samples detected?
    .groups = "drop"
  ) %>%
  ungroup()

# Total detects per species (for ordering)
detect_order <- df_detects %>%
  group_by(Gene) %>%
  summarise(n_detects = sum(n_detects_type)
         #, .groups = "drop"
         ) %>%
  arrange(desc(n_detects))

# compute detection + average rel_ab
df_ranks_type <- df_detects %>%
  left_join(detect_order, by = "Gene") %>%
  
mutate(ID = rank(-n_detects_type, ties.method = "first")) %>%
  
  mutate(
    Gene = factor(
      Gene,
      levels = unique(Gene[order(n_detects, -n_detects_type)])
    )
  ) %>% 
  
  mutate(
    TypeLabel = case_when(
      Type == "H_WW" ~ "Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain" ~ "Sink Biofilm",
      .default = Type
    ),
    TypeLabel = factor(TypeLabel, 
                       levels = c("Wastewater", 
                                  "Sewer Biofilm", 
                                  "Sink Biofilm")),
    type_num = as.numeric(TypeLabel)
  ) %>% 
  
  mutate(
    Gene_label = fct_reorder(Gene, n_detects, .desc = FALSE)
  )
  


#Plot-------------
max_detect <- max(df_ranks_type$n_detects)

# candidate step
raw_step <- max_detect / 4

# function to snap to nice values (2, 5, 10, etc.)
nice_step <- function(x) {
  base <- 10^floor(log10(x))
  if (x/base <= 2) return(2 * base)
  if (x/base <= 5) return(5 * base)
  return(10 * base)
}

break_step <- 5


bars <- 
  df_ranks_type %>% 
  ggplot(aes(Gene_label, -n_detects_type)) +
    geom_col(aes(fill = TypeLabel),
             #color = "white",
             size = .5,
             width = 1.02) +
    # geom_curve(aes(x = 51.2, xend = 47, 
    #                y = -148, yend = -166),
    #            curvature = -.4) +
    # annotate("text", x = 47, y = -300, 
    #          label = "Each rectangle represents\none song included in the\nBBC ranking, its length\n the total points and the\ncolor indicates the rank",
    #          family = "Chivo",
    #          size = 2,
    #          lineheight = .9) +
    # annotate("text", x = 21.5, y = -120, 
    #          label = 'The Top Artists featured in the BBC´s\n"Greatest Hip-Hop Songs of All Time"',
    #          family = "Chivo",
    #          fontface = "bold",
    #          size = 4,
    #          lineheight = .9) +
    # annotate("text", x = 17, y = -120,
    #          label = 'In Autumn 2019, 108 hip-hop and music experts ranked their 5 favorites out of\n311 nominated songs in an online survey by the BBC. The graphic shows points\nscored in total and per song for the top ranked artists and broken down by era.',
    #          family = "Chivo",
    #          fontface = "bold",
    #          color = "grey30",
    #          size = 2.5,
    #          lineheight = .9) +
    coord_flip(clip = "off") +
    scale_x_discrete(position = "top") +
    scale_y_continuous(
      expand = c(.02, .02),
      limits = c(
        -max_detect - max_detect/4,
        0
      ),
      breaks = seq(
        -ceiling(max_detect / break_step) * break_step,
        0,
        by = break_step
      ),
      labels = function(x) abs(as.integer(x)),
      position = "right"
    )+
    # nord::scale_fill_nord(palette = "halifax_harbor",discrete = T,reverse = F
    #                       #guide = F
    #                       ) +
    # scale_fill_viridis_d(option = "D")+
   scale_fill_manual(values = c("#6d7ecd", "#9fdeca", "#c86c69"))+
    #scale_fill_manual(values = fortyfive_pal)+
    theme(
          axis.text.y.right = element_blank(),
          axis.ticks.x = element_line(color = "black"),
          axis.ticks.length = unit(0.5, "pt"),
          #axis.text.y.right = element_markdown(hjust = .5,size = 10),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 16),
          legend.position = "left",
          legend.text = element_text(size = 17),
          plot.margin = margin(5, 0, 5, 5)) +
    labs(x = NULL, 
         fill = "",
         y = "Number of ARG detections")

dots <- df_ranks_type %>% 

  ggplot(aes(Gene_label, type_num, group = Gene_label)) +
    geom_point(aes(Gene_label, 1), color = "white", size = 2,show.legend = FALSE) +
    geom_point(aes(Gene_label, 2), color = "white", size = 2,show.legend = FALSE,) +
    geom_point(aes(Gene_label, 3), color = "white", size = 2,show.legend = FALSE) +
  
    geom_segment(aes(x = Gene_label, xend = Gene_label, 
                     y = 1, yend = 3), 
                 color = "white",
                 size = .06) +
    geom_line(color = "black",
              show.legend = FALSE,
              size = .03) +
    geom_point(aes(fill = TypeLabel, size = n_detects_type),
               show.legend = FALSE,
               shape = 21,
               color = "black",
               stroke = 0.3) +
  # geom_point(aes(fill = TypeLabel, size = n_detects_type), 
  #          shape = 21, 
  #          color = "black", 
  #          stroke = 1) +

    #geom_curve(aes(x = 47, xend = 51, 
    # geom_curve(aes(x = 45, xend = 49, 
    #                y = 6.1, yend = 4.3),
    #            curvature = .4) +
    #annotate("text", x = 45.1, y = 6.1, 
    # annotate("text", x = 43.1, y = 6.1, 
    #          label = "The dot size indicates\nthe number of songs,\nthe dot color the best\nrank in each era",
    #          family = "Chivo",
    #          size = 3.8, 
    #          lineheight = .9) +
    coord_flip() +
    scale_y_continuous(limits = c(.5, 7.3),
                       breaks = 1:3,
                       labels = c("Wastewater", 
                                  "Sewer Biofilm", 
                                  "Sink Biofilm"), 
                       position = "right") +
  scale_size(range = c(0.9, 2.2), guide = FALSE)+
  #scale_size(range = c(2, 5.5), guide = F) +
  scale_fill_manual(values = c("#6d7ecd", "#9fdeca", "#c86c69"))+
  # scale_fill_viridis_d(option = "D")+
    # nord::scale_fill_nord(palette = "halifax_harbor", discrete = T, reverse = F, guide = F) +
  #scale_fill_manual(values = fortyfive_pal)+
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 11, 
                                     #face = "bold", 
                                     hjust = .1, vjust = 0, angle = 45),
          #axis.text.y = element_blank(),
          axis.text.y = element_text(size = 12),
          plot.margin = margin(5, 5, 5, 0),
          plot.caption = element_text(face = "bold", color = "grey30", 
                                      size = 7, margin = margin(t = 15))) +
  guides(shape = "none", color = "none", size = "none")+
    labs(x = NULL, y = NULL,
         caption = "")



plot = bars + dots + plot_layout(widths = c(1, .5))

#plot_spacer() + bars + dots + plot_layout(widths = c(1, .35), heights = c(.1, 1))
# Tell showtext exactly what resolution you'll be using later
showtext_opts(dpi = 300)


ggsave(path = "Biofilm Project Figures",file="CARD-Upset Plot Detections per Top ARGs-only relevant-ONLY SHARED DATES WITH BIOFILM.png", plot=plot, width=7, height= 5.5,dpi = 300)
ggsave(path = "Biofilm Project Figures",file="CARD-Upset Plot Detections per Top ARGs-only relevant-ONLY SHARED DATES WITH BIOFILM.svg", plot=plot, width=7, height=5.5)

print(plot)

```

##Upset Plot - Drug Classes by Presence Absence- All Sites in One-ONLY SHARED DATES WITH BIOFILM
```{r eval=FALSE, include=FALSE}
#theme set----------
font_add_google("Chivo", "Chivo")
font_add_google("Passion One", "Passion One")

## ggplot theme
theme_set(theme_minimal(base_family = "Chivo"))

theme_update(axis.text.x = element_text(size = 11, color = "grey20"),
             axis.text.y = element_text(size = 13, color = "black", face = "bold"),
             axis.ticks.x = element_line(color = "grey45"),
             axis.ticks.y = element_blank(),
             axis.ticks.length.x = unit(.4, "lines"),
             panel.grid = element_blank(),
             plot.background = element_rect(fill = "grey60", color = "grey60"))

font_add_google("Reem Kufi", "Reem Kufi")   # <-- REQUIRED
showtext_auto()                             # <-- activates showtext


## ggplot theme
theme_set(theme_minimal(base_family = "Reem Kufi", base_size = 12))

theme_update(
  plot.title = element_text(size = 27,
                            face = "bold",
                            hjust = .5,
                            margin = margin(10, 0, 30, 0)),
  plot.caption = element_text(size = 9,
                              color = "grey40",
                              hjust = .5,
                              margin = margin(20, 0, 5, 0)),
  axis.text.y = element_blank(),
  axis.title = element_blank(),
  plot.background = element_rect(fill = "grey88", color = NA),
  panel.background = element_rect(fill = NA, color = NA),
  panel.grid = element_blank(),
  panel.spacing.y = unit(0, "lines"),
  strip.text.y = element_blank(),
  legend.position = "bottom",
  legend.text = element_text(size = 9, color = "grey40"),
  legend.box.margin = margin(t = 30), 
  legend.background = element_rect(color = "grey40", 
                                   size = .3, 
                                   fill = "grey95"),
  legend.key.height = unit(.25, "lines"),
  legend.key.width = unit(2.5, "lines"),
  plot.margin = margin(rep(20, 4))
)

#R prep---------------------------------------
## packages
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggtext)
library(showtext)
library(sysfonts)
library(ggplot2)
library(ggstatsplot)
library(dplyr)
library(tidyverse)
library(readxl)
library(vegan)
library(OTUtable)
library(vegan)
library(MASS)
library(ggtext)
library(ggupset)
library(forcats)

font_add_google("Chivo", "Chivo")
font_add_google("Passion One", "Passion One")

## ggplot theme
theme_set(theme_minimal(base_family = "Chivo"))

theme_update(axis.text.x = element_text(size = 11, color = "grey10"),
             axis.text.y = element_text(size = 13, color = "black", face = "bold"),
             axis.ticks.x = element_line(color = "white"),
             axis.ticks.y = element_blank(),
             axis.ticks.length.x = unit(.4, "lines"),
             title = element_text(size = 55),
             panel.grid = element_blank(),
             plot.background = element_rect(fill = "white", color = "white"))

#R data -------------------------
dates = c("2024-09-09", "2024-09-11", "2024-09-12", "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09")

# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame



#final_list = c( "bla-CTX-M", "tetA", "vanA","sul1")

analysis_df <- card_merged_w_ont_df %>%
  filter(
    Type %in% c("Drain", "Endcap", "H_WW"),
    Location %in% locations
  ) %>%
  mutate(
    Date   = as.Date(Date),
    rel_ab = as.numeric(rel_abundance),
    Type   = factor(Type, levels = c("H_WW","Endcap","Drain" ))
  ) %>%
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>%
  drop_na(timepoint)


biofilm_dates <- analysis_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
  mutate(
    biofilm_ref_date = Date
  ) %>% 
    ungroup() %>% 
  dplyr::distinct(Location,timepoint,biofilm_ref_date)

filtered_df <- analysis_df %>%
  mutate(rel_ab = as.numeric(rel_abundance),
         Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  ) %>% 
  
  dplyr::select(-best_date,-date_diff_days,-biofilm_ref_date)


#dataframe
df <- filtered_df %>% 
  
  
  #clean ontology df
  dplyr::select(card_edited_drug_class,SampleID, Date, Location, Type,rel_ab
                #`Drug Class`, `Resistance Mechanism`, class, mechanism, group,
                ) %>% 

  
  # #filter dates
  # filter(
  #   (Type == "H_WW" & Date %in% dates) |
  #   (Type != "H_WW")
  # ) %>% 
  
  mutate(
    TypeLabel = case_when(Type == "H_WW" ~ "Wastewater", 
                          Type == "Drain" ~ "Sink Biofilm", 
                          Type == "Endcap" ~ "Sewer Biofilm", 
                          TRUE ~ as.character(Type)),
    Month = format(as.Date(Date), "%m-%Y")
  ) %>% 

  
  #distinct genes per sample id, date, sample type, etc.
  distinct(card_edited_drug_class,SampleID,Date,Location,Type,rel_ab,TypeLabel,Month)
         
  
  
# Count detects per gene × Type
df_detects <- df %>%
  mutate(detect = rel_ab > 0) %>% 
  
  #distinct genes per sample id, date, sample type, etc.
  distinct(card_edited_drug_class,SampleID,Date,Location,Type,detect,TypeLabel,Month) %>% 
  
  group_by(card_edited_drug_class, Type) %>%
  summarise(
    n_detects_type = sum(detect, na.rm = TRUE),   # how many samples detected?
    .groups = "drop"
  ) %>%
  ungroup()

# Total detects per species (for ordering)
detect_order <- df_detects %>%
  group_by(card_edited_drug_class) %>%
  summarise(n_detects = sum(n_detects_type)
         #, .groups = "drop"
         ) %>%
  arrange(desc(n_detects))

# compute detection + average rel_ab
df_ranks_type <- df_detects %>%
  left_join(detect_order, by = "card_edited_drug_class") %>%
  
mutate(ID = rank(-n_detects_type, ties.method = "first")) %>%
  
  mutate(
    card_edited_drug_class = factor(
      card_edited_drug_class,
      levels = unique(card_edited_drug_class[order(n_detects, -n_detects_type)])
    )
  ) %>% 
  
  mutate(
    TypeLabel = case_when(
      Type == "H_WW" ~ "Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain" ~ "Sink Biofilm",
      .default = Type
    ),
    TypeLabel = factor(TypeLabel, 
                       levels = c("Wastewater", 
                                  "Sewer Biofilm", 
                                  "Sink Biofilm")),
    type_num = as.numeric(TypeLabel)
  ) %>% 
  
  mutate(
    class_label = fct_reorder(card_edited_drug_class, n_detects, .desc = FALSE)
  )
  


#Plot-------------
max_detect <- max(df_ranks_type$n_detects)

# candidate step
raw_step <- max_detect / 4

# function to snap to nice values (2, 5, 10, etc.)
nice_step <- function(x) {
  base <- 10^floor(log10(x))
  if (x/base <= 2) return(2 * base)
  if (x/base <= 5) return(5 * base)
  return(10 * base)
}

break_step <- 10


bars <- 
  df_ranks_type %>% 
  ggplot(aes(class_label, -n_detects_type)) +
    geom_col(aes(fill = TypeLabel),
             #color = "white",
             size = .5,
             width = 1.02) +
    # geom_curve(aes(x = 51.2, xend = 47, 
    #                y = -148, yend = -166),
    #            curvature = -.4) +
    # annotate("text", x = 47, y = -300, 
    #          label = "Each rectangle represents\none song included in the\nBBC ranking, its length\n the total points and the\ncolor indicates the rank",
    #          family = "Chivo",
    #          size = 2,
    #          lineheight = .9) +
    # annotate("text", x = 21.5, y = -120, 
    #          label = 'The Top Artists featured in the BBC´s\n"Greatest Hip-Hop Songs of All Time"',
    #          family = "Chivo",
    #          fontface = "bold",
    #          size = 4,
    #          lineheight = .9) +
    # annotate("text", x = 17, y = -120,
    #          label = 'In Autumn 2019, 108 hip-hop and music experts ranked their 5 favorites out of\n311 nominated songs in an online survey by the BBC. The graphic shows points\nscored in total and per song for the top ranked artists and broken down by era.',
    #          family = "Chivo",
    #          fontface = "bold",
    #          color = "grey30",
    #          size = 2.5,
    #          lineheight = .9) +
    coord_flip(clip = "off") +
    scale_x_discrete(position = "top") +
    scale_y_continuous(
      expand = c(.02, .02),
      limits = c(
        -max_detect - max_detect/4,
        0
      ),
      breaks = seq(
        -ceiling(max_detect / break_step) * break_step,
        0,
        by = break_step
      ),
      labels = function(x) abs(as.integer(x)),
      position = "right"
    )+
    # nord::scale_fill_nord(palette = "halifax_harbor",discrete = T,reverse = F
    #                       #guide = F
    #                       ) +
    #scale_fill_viridis_d(option = "D")+
   scale_fill_manual(values = c("#6d7ecd", "#9fdeca", "#c86c69"))+
    #scale_fill_manual(values = fortyfive_pal)+
    theme(
          axis.text.y.right = element_blank(),
          #axis.text.y.right = element_markdown(hjust = .5,size = 10),
          axis.ticks.x = element_line(color = "black"),
          axis.ticks.length = unit(0.5, "pt"),
          axis.text.x = element_text(size = 44),
          axis.title = element_text(size = 62,face = "bold"),
          legend.position = "left",
          legend.text = element_text(size = 58),
          plot.margin = margin(5, 0, 5, 5)) +
    labs(x = NULL, 
         fill = "",
         y = "Number of ARG detections")

dots <- df_ranks_type %>% 

  ggplot(aes(class_label, type_num, group = class_label)) +
    geom_point(aes(class_label, 1), color = "white", size = 2,show.legend = FALSE) +
    geom_point(aes(class_label, 2), color = "white", size = 2,show.legend = FALSE,) +
    geom_point(aes(class_label, 3), color = "white", size = 2,show.legend = FALSE) +
  
    geom_segment(aes(x = class_label, xend = class_label, 
                     y = 1, yend = 3), 
                 color = "white",
                 size = .06) +
    geom_line(color = "black",
              show.legend = FALSE,
              size = .03) +
    geom_point(aes(fill = TypeLabel, size = n_detects_type),
               show.legend = FALSE,
               shape = 21,
               color = "black",
               stroke = 0.3) +
  # geom_point(aes(fill = TypeLabel, size = n_detects_type), 
  #          shape = 21, 
  #          color = "black", 
  #          stroke = 1) +

    #geom_curve(aes(x = 47, xend = 51, 
    # geom_curve(aes(x = 45, xend = 49, 
    #                y = 6.1, yend = 4.3),
    #            curvature = .4) +
    #annotate("text", x = 45.1, y = 6.1, 
    # annotate("text", x = 43.1, y = 6.1, 
    #          label = "The dot size indicates\nthe number of songs,\nthe dot color the best\nrank in each era",
    #          family = "Chivo",
    #          size = 3.8, 
    #          lineheight = .9) +
    coord_flip() +
    scale_y_continuous(limits = c(.5, 7.3),
                       breaks = 1:3,
                       labels = c("Wastewater", 
                                  "Sewer Biofilm", 
                                  "Sink Biofilm"), 
                       position = "right") +
  scale_size(range = c(0.9, 2.2), guide = FALSE)+
  #scale_size(range = c(2, 5.5), guide = F) +
  # scale_fill_viridis_d(option = "D")+
  scale_fill_manual(values = c("#6d7ecd", "#9fdeca", "#c86c69"))+
    #nord::scale_fill_nord(palette = "halifax_harbor", discrete = T, reverse = F, guide = F) +
  # scale_fill_manual(values = fortyfive_pal)+
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 43, face = "bold", 
                                     hjust = .1, vjust = 0, angle = 35),
          #axis.text.y = element_blank(),
          axis.text.y = element_text(size = 45),
          plot.margin = margin(5, 5, 5, 0),
          plot.caption = element_text(face = "bold", color = "grey30", 
                                      size = 32, margin = margin(t = 15))) +
  guides(shape = "none", color = "none", size = "none")+
    labs(x = NULL, y = NULL,
         caption = "")



plot = bars + dots + plot_layout(widths = c(1, .35))

#plot_spacer() + bars + dots + plot_layout(widths = c(1, .35), heights = c(.1, 1))



ggsave(path = "Biofilm Project Figures",file="CARD-Upset Plot Detections per drug classes-ONLY SHARED DATES WITH BIOFILM.png", plot=plot, width=13, height= 8,dpi = 900)
ggsave(path = "Biofilm Project Figures",file="CARD-UpSet plot of drug classes-ONLY SHARED DATES WITH BIOFILM.svg", plot=plot, width=13, height=8)

print(plot)

```


#Stacked Bars
##stacked bars  sink biofilm v sewer biofilm v ww (facet nested by Location and Type) only top 25 species averaged in sewer branch biofilm across all locations -- x axis is unique sample type and date combination as row numbers 
```{r}


# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame


library(dplyr)
library(ggplot2)
library(forcats)
library(viridis)
library(ggh4x)

# ---------- STEP 1: Extract and clean ----------


df <- card_merged_w_ont_df %>% 
  
  #clean ontology df
  dplyr::select(Gene,SampleID, Date, Location, Type,rel_abundance,`Drug Class`, `Resistance Mechanism`, class, mechanism, group) %>% 
  
  mutate(rel_ab = as.numeric(rel_abundance)) %>% 
  filter(Location %in% locations) %>% 
  
  #filter args in rpip list
  #filter(Species %in% c(amr_genes)) %>% 
  
  #filter for sample type
  filter(Type %in% c("Drain", "Endcap", "H_WW")) %>% 
  
  #
  mutate(
    Type = factor(Type, levels = c("Drain", "Endcap", "H_WW")),
    TypeLabel = case_when(Type == "H_WW" ~ "Hospital Wastewater", 
                          TRUE ~ as.character(Type)),
    Month = format(as.Date(Date), "%m-%Y")
  ) 


# ---------- STEP 2: Collapse into top 25 species ----------
# get top 25 globally by total rel_ab
top25_args <- df %>%
  filter(Type == "Endcap") %>% 
  group_by(Gene) %>%
  summarise(total_ab = mean(rel_ab, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_ab)) %>%
  slice_head(n = 11) %>%
  pull(Gene)

# assign species or "Other"
plot_df <- df %>%
  dplyr::select(Gene, SampleID, Date, Location, Type, rel_ab) %>%
  mutate(
    Type = case_when(
      Type == "H_WW"   ~ "Hospital Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain"  ~ "Sink Biofilm",
      .default = Type
    ),
    class = factor(ifelse(Gene %in% top25_args, Gene, "Other")),
    Type  = factor(Type, levels = c("Sink Biofilm", "Sewer Biofilm", "Hospital Wastewater"))
  )

# ---------- STEP 3: Plot ----------

library(RColorBrewer)
ref = "Other"
length(levels(plot_df$class))
myColors <- brewer.pal(length(levels(plot_df$class)),"Set3")
names(myColors) <- levels(plot_df$class)
myColors[names(myColors)==ref] <- "grey"


library(dplyr)
library(ggplot2)
library(ggh4x)
library(ggtext)   # for element_markdown

# 1) Build a per-facet ordering of samples (one row per unique sample)
sample_order <- plot_df %>%
  distinct(Location, Type, SampleID, Date) %>%
  group_by(Location, Type) %>%
  arrange(Date, SampleID, .by_group = TRUE) %>%
  mutate(
    Event = row_number(),                                # 1..N per facet (no gaps)
    EventLabel = format(as.Date(Date), "%Y-%m-%d")       # optional: a label you can show
  ) %>%
  ungroup()

# 2) Join that back so every row has its facet-local Event index
plot_df2 <- plot_df %>%
  left_join(sample_order, by = c("Location","Type","SampleID","Date"))




# 3) Plot using Event (factor) on x → no gaps within each facet
ggplot(plot_df2, aes(x = factor(Event), y = rel_ab, fill = class)) +
  geom_col(position = "stack") +
  facet_nested(Location ~ Type, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = myColors
                    ) +
  labs(x = NULL, y = "Fraction of Annotated Reads", fill = "ARGs") +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_blank(),
    panel.spacing    = unit(0.05, "lines"),
    legend.position  = "bottom",
    legend.key.size  = unit(0.05, "cm"),
    legend.text      = ggtext::element_markdown(size = 8)   # allow markdown (italics) in legend
  )







```

##stacked bars  sink biofilm v sewer biofilm v ww (facet nested by Location and Type) genes across all locations -- x axis is unique sample type and date combination as row numbers -- this time out of 100- total number of pathogens as the denominator 
```{r}


# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame


# AMR genes
amr_genes <- c("blaTEM", "blaSHV", "blaCTX-M", "blaOXA", "blaKPC", "blaNDM", "blaVIM", "blaIMP",
               "mecA", "vanA", "vanB", "ermB", "tetM", "tetO", "aac(3)-I", "aph(3')-IIIa",
               "strA", "strB","sul1", "sul2", "dfrA", "dfrB", "qnrA", "qnrB", "qnrS")


library(dplyr)
library(ggplot2)
library(forcats)
library(viridis)
library(ggh4x)

# ---------- STEP 1: Extract and clean ----------

df <- card_merged_w_ont_df %>% 
  
  #clean ontology df
  dplyr::select(Gene,SampleID, Date, Location, Type,rel_abundance,`Drug Class`, `Resistance Mechanism`, class, mechanism,card_edited_drug_class, group) %>% 
  
  mutate(rel_ab = as.numeric(rel_abundance)) %>% 
  filter(Location %in% locations) %>% 
  
  #filter args in rpip list
  #filter(Species %in% c(amr_genes)) %>% 
  
  #filter for sample type
  filter(Type %in% c("Drain", "Endcap", "H_WW")) %>% 
  
  #
  mutate(
    Type = factor(Type, levels = c("Drain", "Endcap", "H_WW")),
    TypeLabel = case_when(Type == "H_WW" ~ "Hospital Wastewater", 
                          TRUE ~ as.character(Type)),
    Month = format(as.Date(Date), "%m-%Y")
  ) %>% 
  
  # make denominator of relative abundance total counts for rpip targets
  group_by(SampleID, Location, Type) %>%
  mutate(rel_ab = rel_ab / sum(rel_ab, na.rm = TRUE)) %>%
  ungroup()

# ---------- STEP 2: Collapse into top 25 species ----------
# get top ___ globally by total rel_ab
top25_args <- df %>%
  filter(Type == "Endcap") %>% 
  group_by(Gene) %>%
  summarise(total_ab = mean(rel_ab, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_ab)) %>%
  slice_head(n = 10) %>%
  pull(Gene)

# assign species or "Other"
plot_df <- df %>%
  dplyr::select(Gene, SampleID, Date, Location, Type, rel_ab) %>%
  mutate(
    Type = case_when(
      Type == "H_WW"   ~ "Hospital Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain"  ~ "Sink Biofilm",
      .default = Type
    ),
    class = factor(ifelse(Gene %in% top25_args, Gene, "Other")),
    Type  = factor(Type, levels = c("Sink Biofilm", "Sewer Biofilm", "Hospital Wastewater"))
  )

# ---------- STEP 3: Plot ----------

#colors
# ref = "Other"
# length(levels(plot_df$class))
# myColors <- brewer.pal(length(levels(plot_df$class)),"Set3")
# names(myColors) <- levels(plot_df$class)
# myColors[names(myColors)==ref] <- "grey"

#viridis colors
# ref = "Other"
# myColors <- viridis(length(levels(plot_df$class)),option = "A")
# names(myColors) <- levels(plot_df$class)
# myColors[names(myColors)==ref] <- "grey"


poster_colors = c("#dc9298","#a51c2f","#3d3d3d","#f3e3e5","#417c8c","#e57262","#915493")

poster_colors = c("#ce2e2d","#c00300","#88b8c5","#dd9298","#ab74a8","#fc3b3b","#a61d30","#eededf","#dc6c64","#467d91","#915493","#3d3d3d")

poster_colors = c("#ce2e2d","#88b8c5","#dd9298","#efc1c0","#fc3b3b","#a61d30","#eededf","#dc6c64","#915493","#3d3d3d")

# Assign colors: first 15 get a palette, "Other" gets grey
# ref = "Other"
# myColors <- poster_colors[1:length(levels(plot_df$class))]
# names(myColors) <- levels(plot_df$class)
# myColors[names(myColors)==ref] <- "grey"


# Assign colors: 
ref = "Other"
myColors <- fortyfive_pal[1:length(levels(plot_df$class))]
names(myColors) <- levels(plot_df$class)
myColors[names(myColors)==ref] <- "grey"



##actual plotting-----------
library(dplyr)
library(ggplot2)
library(ggh4x)
library(ggtext)   # for element_markdown

# 1) Build a per-facet ordering of samples (one row per unique sample)
sample_order <- plot_df %>%
  distinct(Location, Type, SampleID, Date) %>%
  group_by(Location, Type) %>%
  arrange(Date, SampleID, .by_group = TRUE) %>%
  mutate(
    Event = row_number(),                                # 1..N per facet (no gaps)
    EventLabel = format(as.Date(Date), "%Y-%m-%d")       # optional: a label you can show
  ) %>%
  ungroup()

# 2) Join that back so every row has its facet-local Event index
plot_df2 <- plot_df %>%
  left_join(sample_order, by = c("Location","Type","SampleID","Date"))



# Plot using Event (factor) on x → no gaps within each facet


library(scales)  # for percent_format()

plot_df2 <- plot_df2 %>%
  mutate(Type = case_when(Type =="Sink Biofilm" ~ "Sink\nBiofilm",
                       Type =="Sewer Biofilm" ~ "Sewer\nBiofilm",
                       Type == "Hospital Wastewater" ~ "Hospital\nWastewater",
                       .default = Type))


plot = ggplot(plot_df2, aes(x = factor(Event), y = rel_ab, fill = class)) +
  geom_col(position = "stack") +
  facet_nested(Location ~ Type, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = myColors) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75),   # tick marks at 25, 50, 75
    labels = c("0","25", "50", "75"),   # plain numbers, no %
    limits = c(0, 1)                 # optional: keep to 100%
  ) +
  labs(x = NULL, y = "Percent of Annotated Reads", fill = "ARGs") +
theme_minimal(base_size = 11) +
theme(
  panel.background = element_rect(fill = "white", colour = "black"),
  strip.background = element_rect(fill = "white", colour = "black"),
  strip.text = element_text(
    face = "bold", 
    hjust = 0.5,          # center horizontally
    vjust = 0.5,          # center vertically
    lineheight = 0.9      # tighter line spacing if wrapped
  ),
  axis.text.x        = element_blank(),
  panel.spacing      = unit(0.05, "lines"),
  legend.position    = "bottom",
  legend.key.size    = unit(0.5, "cm"),
  legend.text        = ggtext::element_markdown(size = 8),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  #panel.grid.major.y = element_line(color = "grey80"),
  panel.grid.minor.y = element_blank(),
  panel.ontop        = FALSE   # <- ensures geoms (bars) are drawn over grid lines
)


print(plot)

  ggsave(path = "Biofilm Project Figures",file="CARD stacked bar plot of ARGs detected.png", plot=plot, width=9, height=5,dpi = 400)
  ggsave(path = "Biofilm Project Figures",file="CARD stacked bar plot of ARGs detected.svg", plot=plot, width=9, height=5)
  
  
##now vertical-----------
  
  plot = ggplot(plot_df2, aes(y = factor(Event), x = rel_ab, fill = class)) +
  geom_col(position = "stack") +
  facet_nested(Type ~ Location, scales = "free_y", space = "free_y",switch = "y") +
  scale_fill_manual(values = myColors,
                    guide = guide_legend(
    nrow = 4,              # number of rows in the legend
    byrow = TRUE,          # fill across rows instead of down columns
    title.position = "top" # keep the title above
  )
  ) +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.50, 0.75),   # tick marks at 25, 50, 75
    labels = c("0","25", "50", "75"),  # plain numbers, no %
    limits = c(0, 1)
  ) +
  labs(y = NULL, x = "Percent of Annotated Reads", fill = "ARGs") +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(
      face = "bold", 
      hjust = 0.5,
      vjust = 0.5,
      lineheight = 0.9
    ),
    axis.text.y        = element_blank(),  # hide if crowded
    panel.spacing      = unit(0.05, "lines"),
    legend.position    = "bottom",
    legend.key.size    = unit(0.5, "cm"),
    legend.text        = ggtext::element_markdown(size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.ontop        = FALSE
  )

  
  print(plot)

  ggsave(path = "Biofilm Project Figures",file="CARD stacked bar plot of ARGs detected vertical.png", plot=plot, width=6.5, height=10,dpi = 500)
  ggsave(path = "Biofilm Project Figures",file="CARD stacked bar plot of ARGs detected vertical.svg", plot=plot, width=5, height=9)


##now by ARG class-----------
  
  # get top 25 globally by total rel_ab
top25_classes <- df %>%
  dplyr::select(Gene, SampleID, Date, Location, Type, rel_ab, class) %>%
  drop_na(Type) %>%
  
  drop_na(class) %>% 
  
  filter(Type == "Endcap") %>% 
  group_by(class) %>%
  summarise(total_ab = mean(rel_ab, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_ab)) %>%
  slice_head(n = 10) %>%
  pull(class)

# assign species or "Other"
plot_df <- df %>%
  dplyr::select(Gene, SampleID, Date, Location, Type, rel_ab, class) %>%
  drop_na(Type) %>%
  
  drop_na(class) %>% 
  
  #dplyr::select(Gene, SampleID, Date, Location, Type, rel_ab) %>%
  mutate(
    Type = case_when(
      Type == "H_WW"   ~ "Hospital Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain"  ~ "Sink Biofilm",
      .default = Type
    ),
    class = factor(ifelse(class %in% top25_classes, class, "Other")),
    Type  = factor(Type, levels = c("Sink Biofilm", "Sewer Biofilm", "Hospital Wastewater"))
  )


#colors
# ref = "Other"
# length(levels(plot_df$class))
# myColors <- brewer.pal(length(levels(plot_df$class)),"Set3")
# names(myColors) <- levels(plot_df$class)
# myColors[names(myColors)==ref] <- "grey"

#viridis colors
# ref = "Other"
# myColors <- viridis(length(levels(plot_df$class)),option = "A")
# names(myColors) <- levels(plot_df$class)
# myColors[names(myColors)==ref] <- "grey"


poster_colors = c("#dc9298","#a51c2f","#3d3d3d","#f3e3e5","#417c8c","#e57262","#915493")

poster_colors = c("#ce2e2d","#c00300","#88b8c5","#dd9298","#ab74a8","#fc3b3b","#a61d30","#eededf","#dc6c64","#467d91","#915493","#3d3d3d")

poster_colors = c("#ce2e2d","#88b8c5","#dd9298","#efc1c0","#fc3b3b","#a61d30","#eededf","#dc6c64","#915493","#3d3d3d")

# Assign colors: first 15 get a palette, "Other" gets grey
ref = "Other"
myColors <- fortyfive_pal[1:length(levels(plot_df$class))]
names(myColors) <- levels(plot_df$class)
myColors[names(myColors)==ref] <- "grey"



library(dplyr)
library(ggplot2)
library(ggh4x)
library(ggtext)   # for element_markdown

# 1) Build a per-facet ordering of samples (one row per unique sample)
sample_order <- plot_df %>%
  distinct(Location, Type, SampleID, Date) %>%
  group_by(Location, Type) %>%
  arrange(Date, SampleID, .by_group = TRUE) %>%
  mutate(
    Event = row_number(),                                # 1..N per facet (no gaps)
    EventLabel = format(as.Date(Date), "%Y-%m-%d")       # optional: a label you can show
  ) %>%
  ungroup()

# 2) Join that back so every row has its facet-local Event index
plot_df2 <- plot_df %>%
  left_join(sample_order, by = c("Location","Type","SampleID","Date"))



# Plot using Event (factor) on x → no gaps within each facet


library(scales)  # for percent_format()

plot_df2 <- plot_df2 %>%
  mutate(Type = case_when(Type =="Sink Biofilm" ~ "Sink\nBiofilm",
                       Type =="Sewer Biofilm" ~ "Sewer\nBiofilm",
                       Type == "Hospital Wastewater" ~ "Hospital\nWastewater",
                       .default = Type))


plot = ggplot(plot_df2, aes(x = factor(Event), y = rel_ab, fill = class)) +
  geom_col(position = "stack") +
  facet_nested(Location ~ Type, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = myColors) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75),   # tick marks at 25, 50, 75
    labels = c("0","25", "50", "75"),   # plain numbers, no %
    limits = c(0, 1)                 # optional: keep to 100%
  ) +
  labs(x = NULL, y = "Percent of Annotated Reads", fill = "") +
theme_minimal(base_size = 11) +
theme(
  panel.background = element_rect(fill = "white", colour = "black"),
  strip.background = element_rect(fill = "white", colour = "black"),
  strip.text = element_text(
    face = "bold", 
    hjust = 0.5,          # center horizontally
    vjust = 0.5,          # center vertically
    lineheight = 0.9      # tighter line spacing if wrapped
  ),
  axis.text.x        = element_blank(),
  panel.spacing      = unit(0.05, "lines"),
  legend.position    = "bottom",
  legend.key.size    = unit(0.5, "cm"),
  legend.text        = ggtext::element_markdown(size = 8),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  #panel.grid.major.y = element_line(color = "grey80"),
  panel.grid.minor.y = element_blank(),
  panel.ontop        = FALSE   # <- ensures geoms (bars) are drawn over grid lines
)


print(plot)

  ggsave(path = "Biofilm Project Figures",file="CARD stacked bar plot of ARGs detected by Drug Class.png", plot=plot, width=9, height=5,dpi = 400)
  ggsave(path = "Biofilm Project Figures",file="CARD stacked bar plot of ARGs detected by Drug Class.svg", plot=plot, width=9, height=5)
  

##now vertical-----------
  
  plot = ggplot(plot_df2, aes(y = factor(Event), x = rel_ab, fill = class)) +
  geom_col(position = "stack") +
  facet_nested(Type ~ Location, scales = "free_y", space = "free_y",switch = "y") +
  scale_fill_manual(values = myColors,
                    guide = guide_legend(
    nrow = 4,              # number of rows in the legend
    byrow = TRUE,          # fill across rows instead of down columns
    title.position = "top" # keep the title above
  )
  ) +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.50, 0.75),   # tick marks at 25, 50, 75
    labels = c("0","25", "50", "75"),  # plain numbers, no %
    limits = c(0, 1)
  ) +
  labs(y = NULL, x = "Percent of Annotated Reads", fill = "") +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(
      face = "bold", 
      hjust = 0.5,
      vjust = 0.5,
      lineheight = 0.9
    ),
    axis.text.y        = element_blank(),  # hide if crowded
    panel.spacing      = unit(0.05, "lines"),
    legend.position    = "bottom",
    legend.key.size    = unit(0.5, "cm"),
    legend.text        = ggtext::element_markdown(size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.ontop        = FALSE
  )

  
  print(plot)

  ggsave(path = "Biofilm Project Figures",file="CARD stacked bar plot of ARGs detected by Drug Class vertical.png", plot=plot, width=6.5, height=10,dpi = 500)
  ggsave(path = "Biofilm Project Figures",file="CARD stacked bar plot of ARGs detected by Drug Class vertical.svg", plot=plot, width=5, height=9)
  
  
##now by CARD ARG Pruden Lab Drug class-----------
  
categorized_drug_class = c("Betalactams", "Sulfonamides", "Trimethoprim", "Tetracyclines", "Aminocoumarin", "Aminoglycoside", "Elfamycin", "Fosfomycin", "Glycopeptides", "MLS", "Peptides", "Phenicol", "Quinolones", "Rifamycines")  

keep_classes  = categorized_drug_class


# assign species or "Other"
plot_df <- df %>%
  dplyr::select(Gene, SampleID, Date, Location, Type, rel_ab, card_edited_drug_class) %>%
  drop_na(Type) %>%
  
  
  #dplyr::select(Gene, SampleID, Date, Location, Type, rel_ab) %>%
  mutate(
    Type = case_when(
      Type == "H_WW"   ~ "Hospital Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain"  ~ "Sink Biofilm",
      .default = Type
    ),
    class = factor(ifelse(card_edited_drug_class %in% keep_classes, card_edited_drug_class, "Other")),
    Type  = factor(Type, levels = c("Sink Biofilm", "Sewer Biofilm", "Hospital Wastewater"))
  )


#colors
# ref = "Other"
# length(levels(plot_df$class))
# myColors <- brewer.pal(length(levels(plot_df$class)),"Set3")
# names(myColors) <- levels(plot_df$class)
# myColors[names(myColors)==ref] <- "grey"

#viridis colors
# ref = "Other"
# myColors <- viridis(length(levels(plot_df$class)),option = "A")
# names(myColors) <- levels(plot_df$class)
# myColors[names(myColors)==ref] <- "grey"


poster_colors = c("#dc9298","#a51c2f","#3d3d3d","#f3e3e5","#417c8c","#e57262","#915493")

poster_colors = c("#ce2e2d","#c00300","#88b8c5","#dd9298","#ab74a8","#fc3b3b","#a61d30","#eededf","#dc6c64","#467d91","#915493","#3d3d3d")

poster_colors = c("#ce2e2d","#88b8c5","#dd9298","#efc1c0","#fc3b3b","#a61d30","#eededf","#dc6c64","#915493","#3d3d3d")

# Assign colors: first 15 get a palette, "Other" gets grey
ref = "Other"
myColors <- fortyfive_pal[1:length(levels(plot_df$class))]
names(myColors) <- levels(plot_df$class)
myColors[names(myColors)==ref] <- "grey"



library(dplyr)
library(ggplot2)
library(ggh4x)
library(ggtext)   # for element_markdown

# 1) Build a per-facet ordering of samples (one row per unique sample)
sample_order <- plot_df %>%
  distinct(Location, Type, SampleID, Date) %>%
  group_by(Location, Type) %>%
  arrange(Date, SampleID, .by_group = TRUE) %>%
  mutate(
    Event = row_number(),                                # 1..N per facet (no gaps)
    EventLabel = format(as.Date(Date), "%Y-%m-%d")       # optional: a label you can show
  ) %>%
  ungroup()

# 2) Join that back so every row has its facet-local Event index
plot_df2 <- plot_df %>%
  left_join(sample_order, by = c("Location","Type","SampleID","Date"))



# Plot using Event (factor) on x → no gaps within each facet


library(scales)  # for percent_format()

plot_df2 <- plot_df2 %>%
  mutate(Type = case_when(Type =="Sink Biofilm" ~ "Sink\nBiofilm",
                       Type =="Sewer Biofilm" ~ "Sewer\nBiofilm",
                       Type == "Hospital Wastewater" ~ "Hospital\nWastewater",
                       .default = Type))


plot = ggplot(plot_df2, aes(x = factor(Event), y = rel_ab, fill = class)) +
  geom_col(position = "stack") +
  facet_nested(Location ~ Type, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = myColors) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75),   # tick marks at 25, 50, 75
    labels = c("0","25", "50", "75"),   # plain numbers, no %
    limits = c(0, 1)                 # optional: keep to 100%
  ) +
  labs(x = NULL, y = "Percent of Annotated Reads", fill = "") +
theme_minimal(base_size = 11) +
theme(
  panel.background = element_rect(fill = "white", colour = "black"),
  strip.background = element_rect(fill = "white", colour = "black"),
  strip.text = element_text(
    face = "bold", 
    hjust = 0.5,          # center horizontally
    vjust = 0.5,          # center vertically
    lineheight = 0.9      # tighter line spacing if wrapped
  ),
  axis.text.x        = element_blank(),
  panel.spacing      = unit(0.05, "lines"),
  legend.position    = "bottom",
  legend.key.size    = unit(0.5, "cm"),
  legend.text        = ggtext::element_markdown(size = 8),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  #panel.grid.major.y = element_line(color = "grey80"),
  panel.grid.minor.y = element_blank(),
  panel.ontop        = FALSE   # <- ensures geoms (bars) are drawn over grid lines
)


print(plot)

  ggsave(path = "Biofilm Project Figures",file="CARD stacked bar plot of ARGs detected by Pruden Lab Drug Class.png", plot=plot, width=9, height=5,dpi = 400)
  ggsave(path = "Biofilm Project Figures",file="CARD stacked bar plot of ARGs detected by Pruden Lab Drug Class.svg", plot=plot, width=9, height=5)
  
#for slides
  
plot = ggplot(plot_df2, aes(x = factor(Event), y = rel_ab, fill = class)) +
  geom_col(position = "stack") +
  facet_nested(Location ~ Type, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = myColors) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75),   # tick marks at 25, 50, 75
    labels = c("0","25", "50", "75"),   # plain numbers, no %
    limits = c(0, 1)                 # optional: keep to 100%
  ) +
  labs(x = NULL, y = "Percent of Annotated Reads", fill = "", title = "") +
theme_minimal(base_size = 11) +
theme(
  panel.background = element_rect(fill = "white", colour = "black"),
  strip.background = element_rect(fill = "white", colour = "black"),
  strip.text = element_text(
    face = "bold", 
    hjust = 0.5,          # center horizontally
    vjust = 0.5,          # center vertically
    lineheight = 0.9      # tighter line spacing if wrapped
  ),
  axis.text.x        = element_blank(),
  panel.spacing      = unit(0.05, "lines"),
  legend.position    = "bottom",
  legend.key.size    = unit(0.1, "cm"),
  legend.spacing.x = unit(0.2, "cm"),   # space between items horizontally
  legend.spacing.y = unit(0.2, "cm"),    # space between rows if wrapped
  legend.text        = ggtext::element_markdown(size = 7.5),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  legend.key.width  = unit(1, "cm"),
  legend.key.height = unit(0.5, "cm"),
  #panel.grid.major.y = element_line(color = "grey80"),
  panel.grid.minor.y = element_blank(),
  panel.ontop        = FALSE   # <- ensures geoms (bars) are drawn over grid lines
)+
  guides(fill = guide_legend(ncol = 3))


print(plot)

  ggsave(path = "Biofilm Project Figures",file="CARD stacked bar plot of ARGs detected by Pruden Lab Drug Class-slides.png", plot=plot, width=7, height=8,dpi = 600)
  ggsave(path = "Biofilm Project Figures",file="CARD stacked bar plot of ARGs detected by Pruden Lab Drug Class-slides.svg", plot=plot, width=7, height=8)
  

##now vertical-----------
  
  plot = ggplot(plot_df2, aes(y = factor(Event), x = rel_ab, fill = class)) +
  geom_col(position = "stack") +
  facet_nested(Type ~ Location, scales = "free_y", space = "free_y",switch = "y") +
  scale_fill_manual(values = myColors,
                    guide = guide_legend(
    nrow = 4,              # number of rows in the legend
    byrow = TRUE,          # fill across rows instead of down columns
    title.position = "top" # keep the title above
  )
  ) +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.50, 0.75),   # tick marks at 25, 50, 75
    labels = c("0","25", "50", "75"),  # plain numbers, no %
    limits = c(0, 1)
  ) +
  labs(y = NULL, x = "Percent of Annotated Reads", fill = "") +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(
      face = "bold", 
      hjust = 0.5,
      vjust = 0.5,
      lineheight = 0.9
    ),
    axis.text.y        = element_blank(),  # hide if crowded
    panel.spacing      = unit(0.05, "lines"),
    legend.position    = "bottom",
    legend.key.size    = unit(0.5, "cm"),
    legend.text        = ggtext::element_markdown(size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.ontop        = FALSE
  )

  
  print(plot)

  ggsave(path = "Biofilm Project Figures",file="CARD stacked bar plot of ARGs detected by Drug Class vertical.png", plot=plot, width=6.5, height=10,dpi = 500)
  ggsave(path = "Biofilm Project Figures",file="CARD stacked bar plot of ARGs detected by Drug Class vertical.svg", plot=plot, width=5, height=9)




```





##stacked bars  sink biofilm v sewer biofilm v ww (facet nested by Location and Type) genes across all locations -- x axis is unique sample type and date combination as row numbers -- this time out of 100- total number of pathogens as the denominator-ONLY SHARED DATES FOR WW WITH BIOFILM
```{r}
#specify dates
dates <- as.Date(c(
  "2024-09-09", "2024-09-11", "2024-09-12",
  "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09"
))


# Specify locations to include
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame


# AMR genes
amr_genes <- c("blaTEM", "blaSHV", "blaCTX-M", "blaOXA", "blaKPC", "blaNDM", "blaVIM", "blaIMP",
               "mecA", "vanA", "vanB", "ermB", "tetM", "tetO", "aac(3)-I", "aph(3')-IIIa",
               "strA", "strB","sul1", "sul2", "dfrA", "dfrB", "qnrA", "qnrB", "qnrS")


library(dplyr)
library(ggplot2)
library(forcats)
library(viridis)
library(ggh4x)

# ---------- STEP 1: Extract and clean ----------

df <- card_merged_w_ont_df %>% 
  
  #clean ontology df
  dplyr::select(Gene,SampleID, Date, Location, Type,rel_abundance,`Drug Class`, `Resistance Mechanism`, class, mechanism,card_edited_drug_class, group) %>% 
  
  mutate(rel_ab = as.numeric(rel_abundance)) %>% 
  
  mutate(Date = as.Date(Date)) %>% 
  
  filter(Location %in% locations) %>% 
  
  #filter args in rpip list
  #filter(Species %in% c(amr_genes)) %>% 
  
  #filter for sample type
  filter(Type %in% c("Drain", "Endcap", "H_WW")) %>% 
  
  # now replace Site names 
  mutate(Location = case_when(Location == "MARIO" ~ "Site 1",
                              Location == "FIONA" ~ "Site 2",
                              Location == "LUIGI" ~ "Site 3",
                              Location == "SHREK" ~ "Site 4",
                              Location == "OSCAR" ~ "Site 5",
                              .default = Location)) %>% 
  
 
  #filter dates
  filter(
    (Type == "H_WW" & Date %in% dates) |
    (Type != "H_WW")
  ) %>% 
  
  #
  mutate(
    Type = factor(Type, levels = c("Drain", "Endcap", "H_WW")),
    TypeLabel = case_when(Type == "H_WW" ~ "Hospital Wastewater", 
                          TRUE ~ as.character(Type)),
    Month = format(as.Date(Date), "%m-%Y")
  ) %>% 
  
  # make denominator of relative abundance total counts for rpip targets
  group_by(SampleID, Location, Type) %>%
  mutate(rel_ab = rel_ab / sum(rel_ab, na.rm = TRUE)) %>%
  ungroup()

# ---------- STEP 2: Collapse into top 25 species ----------
# get top ___ globally by total rel_ab
top25_args <- df %>%
  filter(Type == "Endcap") %>% 
  group_by(Gene) %>%
  summarise(total_ab = mean(rel_ab, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_ab)) %>%
  slice_head(n = 10) %>%
  pull(Gene)

# assign species or "Other"
plot_df <- df %>%
  dplyr::select(Gene, SampleID, Date, Location, Type, rel_ab) %>%
  mutate(
    Type = case_when(
      Type == "H_WW"   ~ "Hospital Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain"  ~ "Sink Biofilm",
      .default = Type
    ),
    class = factor(ifelse(Gene %in% top25_args, Gene, "Other")),
    Type  = factor(Type, levels = c("Sink Biofilm", "Sewer Biofilm", "Hospital Wastewater"))
  )

# ---------- STEP 3: Plot ----------

  
  
##now by CARD ARG Pruden Lab Drug class-----------
  
categorized_drug_class = c("Betalactams", "Sulfonamides", "Trimethoprim", "Tetracyclines", "Aminocoumarin", "Aminoglycoside", "Elfamycin", "Fosfomycin", "Glycopeptides", "MLS", "Peptides", "Phenicol", "Quinolones", "Rifamycines")  

keep_classes  = categorized_drug_class


# assign species or "Other"
plot_df <- df %>%
  dplyr::select(Gene, SampleID, Date, Location, Type, rel_ab, card_edited_drug_class) %>%
  drop_na(Type) %>%
  
  
  #dplyr::select(Gene, SampleID, Date, Location, Type, rel_ab) %>%
  mutate(
    Type = case_when(
      Type == "H_WW"   ~ "Hospital Wastewater",
      Type == "Endcap" ~ "Sewer Biofilm",
      Type == "Drain"  ~ "Sink Biofilm",
      .default = Type
    ),
    class = factor(ifelse(card_edited_drug_class %in% keep_classes, card_edited_drug_class, "Other")),
    Type  = factor(Type, levels = c("Sink Biofilm", "Sewer Biofilm", "Hospital Wastewater"))
  )


#colors
# ref = "Other"
# length(levels(plot_df$class))
# myColors <- brewer.pal(length(levels(plot_df$class)),"Set3")
# names(myColors) <- levels(plot_df$class)
# myColors[names(myColors)==ref] <- "grey"

#viridis colors
# ref = "Other"
# myColors <- viridis(length(levels(plot_df$class)),option = "A")
# names(myColors) <- levels(plot_df$class)
# myColors[names(myColors)==ref] <- "grey"


poster_colors = c("#dc9298","#a51c2f","#3d3d3d","#f3e3e5","#417c8c","#e57262","#915493")

poster_colors = c("#ce2e2d","#c00300","#88b8c5","#dd9298","#ab74a8","#fc3b3b","#a61d30","#eededf","#dc6c64","#467d91","#915493","#3d3d3d")

poster_colors = c("#ce2e2d","#88b8c5","#dd9298","#efc1c0","#fc3b3b","#a61d30","#eededf","#dc6c64","#915493","#3d3d3d")

# Assign colors: first 15 get a palette, "Other" gets grey
ref = "Other"
myColors <- fortyfive_pal[1:length(levels(plot_df$class))]
names(myColors) <- levels(plot_df$class)
myColors[names(myColors)==ref] <- "grey"



library(dplyr)
library(ggplot2)
library(ggh4x)
library(ggtext)   # for element_markdown

# 1) Build a per-facet ordering of samples (one row per unique sample)
sample_order <- plot_df %>%
  distinct(Location, Type, SampleID, Date) %>%
  group_by(Location, Type) %>%
  arrange(Date, SampleID, .by_group = TRUE) %>%
  mutate(
    Event = row_number(),                                # 1..N per facet (no gaps)
    EventLabel = format(as.Date(Date), "%Y-%m-%d")       # optional: a label you can show
  ) %>%
  ungroup()

# 2) Join that back so every row has its facet-local Event index
plot_df2 <- plot_df %>%
  left_join(sample_order, by = c("Location","Type","SampleID","Date"))



# Plot using Event (factor) on x → no gaps within each facet


library(scales)  # for percent_format()

plot_df2 <- plot_df2 %>%
  mutate(Type = case_when(Type =="Sink Biofilm" ~ "Sink\nBiofilm",
                       Type =="Sewer Biofilm" ~ "Sewer\nBiofilm",
                       Type == "Hospital Wastewater" ~ "Hospital\nWastewater",
                       .default = Type))


plot = ggplot(plot_df2, aes(x = factor(Event), y = rel_ab, fill = class)) +
  geom_col(position = "stack") +
  facet_nested(Location ~ Type, scales = "free_x", space = "free_x",switch = "y") +
  scale_fill_manual(values = myColors) +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75),   # tick marks at 25, 50, 75
    labels = c("0","25", "50", "75"),   # plain numbers, no %
    limits = c(0, 1)                 # optional: keep to 100%
  ) +
  labs(x = NULL, y = "Percent of Annotated Reads", fill = "") +
theme_minimal(base_size = 11) +
theme(
  panel.background = element_rect(fill = "white", colour = "black"),
  strip.background = element_rect(fill = "white", colour = "black"),
  strip.text = element_text(
    face = "bold", 
    hjust = 0.5,          # center horizontally
    vjust = 0.5,          # center vertically
    lineheight = 0.9      # tighter line spacing if wrapped
  ),
  axis.text.x        = element_blank(),
  panel.spacing      = unit(0.05, "lines"),
  legend.position    = "bottom",
  legend.key.size    = unit(0.5, "cm"),
  legend.text        = ggtext::element_markdown(size = 4.5),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.major.y = element_blank(),
  #panel.grid.major.y = element_line(color = "grey80"),
  panel.grid.minor.y = element_blank(),
  panel.ontop        = FALSE   # <- ensures geoms (bars) are drawn over grid lines
)+
  guides(
    fill = guide_legend(ncol = 8)  # <- set number of columns here
  )


print(plot)

  ggsave(path = "Biofilm Project Figures",file="CARD stacked bar plot of ARGs detected by Pruden Lab Drug Class-ONLY SHARED DATES WITH BIOFILM.png", plot=plot, width=5, height=5,dpi = 400)
  ggsave(path = "Biofilm Project Figures",file="CARD stacked bar plot of ARGs detected by Pruden Lab Drug Class-ONLY SHARED DATES WITH BIOFILM.svg", plot=plot, width=5, height=5)
  

##now vertical-----------
  
  plot = ggplot(plot_df2, aes(y = factor(Event), x = rel_ab, fill = class)) +
  geom_col(position = "stack") +
  facet_nested(Type ~ Location, scales = "free_y", space = "free_y",switch = "y") +
  scale_fill_manual(values = myColors,
                    guide = guide_legend(
    nrow = 4,              # number of rows in the legend
    byrow = TRUE,          # fill across rows instead of down columns
    title.position = "top" # keep the title above
  )
  ) +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.50, 0.75),   # tick marks at 25, 50, 75
    labels = c("0","25", "50", "75"),  # plain numbers, no %
    limits = c(0, 1)
  ) +
  labs(y = NULL, x = "Percent of Annotated Reads", fill = "") +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", colour = "black"),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(
      face = "bold", 
      hjust = 0.5,
      vjust = 0.5,
      lineheight = 0.9
    ),
    axis.text.y        = element_blank(),  # hide if crowded
    panel.spacing      = unit(0.05, "lines"),
    legend.position    = "bottom",
    legend.key.size    = unit(0.5, "cm"),
    legend.text        = ggtext::element_markdown(size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.ontop        = FALSE
  )

  
  print(plot)

  ggsave(path = "Biofilm Project Figures",file="CARD stacked bar plot of ARGs detected by Drug Class vertical-ONLY SHARED DATES WITH BIOFILm.png", plot=plot, width=6.5, height=10,dpi = 500)
  ggsave(path = "Biofilm Project Figures",file="CARD stacked bar plot of ARGs detected by Drug Class vertical-ONLY SHARED DATES WITH BIOFILM.svg", plot=plot, width=5, height=9)




```





#MaAsLin
##Maaslin3-ONLY SHARED DATES WITH BIOFILM- ARGs
```{r}

#specify dates-------------
dates <- as.Date(c(
  "2024-09-09", "2024-09-11", "2024-09-12",
  "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09"
))

# Specify locations to include--------
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame



#libraries--------
library(ggplot2)
library(readxl)
library(here)
library(cowplot)
library(tidyverse)
library(patchwork)
library(readr)
library(phyloseq)
library(microshades)
library(ape)
library(vegan)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(maaslin3)
library(writexl)
library(ggrepel)

#maaslin3----------------
analysis_df <- card_merged_w_ont_df %>%
  filter(
    Type %in% c("Drain", "Endcap", "H_WW")
  ) %>%
  mutate(
    Date   = as.Date(Date),
    RPKM = as.numeric(RPKM),
    Type   = factor(Type, levels = c("H_WW","Endcap","Drain" ))
  ) %>%
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>%
  drop_na(timepoint)


biofilm_dates <- analysis_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
  mutate(
    biofilm_ref_date = Date
  ) %>% 
    ungroup() %>% 
  dplyr::distinct(Location,timepoint,biofilm_ref_date)

filtered_df <- analysis_df %>%
  mutate(
         Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  ) %>% 
  
  dplyr::select(-best_date,-date_diff_days,-biofilm_ref_date)

df <- card_merged_w_ont_df %>% 
  
  #clean ontology df
  dplyr::select(Gene,SampleID, Date, Location, Type,RPKM,`Drug Class`, `Resistance Mechanism`, class, mechanism,card_edited_drug_class, group) %>% 
  
  mutate(Date = as.Date(Date)) %>% 
  
  filter(Location %in% locations) %>% 
  
  #filter args in rpip list
  #filter(Species %in% c(amr_genes)) %>% 
  
  #filter for sample type
  filter(Type %in% c("Endcap", "H_WW")) %>% 
 
  #filter specific dates shared with biofilm
  filter(
    (Type == "H_WW" & Date %in% dates) |
    (Type != "H_WW")
  ) %>%
  
  #change label names
  mutate(
    Type = factor(Type, levels = c("Endcap", "H_WW")),
    TypeLabel = case_when(Type == "H_WW" ~ "Hospital Wastewater", 
                          TRUE ~ as.character(Type))
  ) 

df_sum <- df %>%
  group_by(SampleID, Gene) %>%
  summarise(
    RPKM = sum(RPKM, na.rm = TRUE),
    .groups = "drop"
  )

otu <- df_sum %>%
  dplyr::select(SampleID, Gene, RPKM) %>%
  tidyr::pivot_wider(
    names_from = Gene,
    values_from = RPKM,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("SampleID")

meta <- df %>%
  dplyr::select(
    SampleID,
    Date,
    Location,
    Type
  ) %>%
  distinct() %>%
  tibble::column_to_rownames("SampleID")

#MaAsLin3------------------------

fit_data <- maaslin3(
  input_data = otu,
  input_metadata = meta,
  output = "maaslin3_output",
  fixed_effects = c("Type"),
  random_effects = c("Location"),
  normalization = "NONE",
  transform = "LOG",
  standardize = TRUE,
  min_prevalence = 0.1,
  verbosity = "ERROR"
)



#plotting it-------------
plot_df <- fit_data$normalized %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "feature",
    values_to = "abundance"
  ) %>%
  left_join(
    fit_data$metadata %>%
      rownames_to_column("SampleID"),
    by = "SampleID"
  )


#volcano plot----------

#get mean abundance directly from MaAslin3 normalized data
# Significant features
sig_abund_joint <- subset(fit_data$fit_data_abundance$results, qval_joint < 0.05)
sig_prev_joint <- subset(fit_data$fit_data_prevalence$results, qval_joint < 0.05)



sig_feats <- union(sig_abund_joint$feature, sig_prev_joint$feature)

top_feat <- sig_abund_joint %>%
  filter(feature %in% sig_feats) %>%
  arrange(qval_joint)

# Compute mean abundance directly from MaAsLin3 normalized table
mean_abund <- colMeans(fit_data$normalized[, top_feat$feature, drop = FALSE], na.rm = TRUE)
top_feat$mean_abund <- mean_abund

# Optionally join extra info if you have a taxonomy table
# top_feat <- left_join(top_feat, top_feat_unique, by = "feature")

# Example: volcano for a subset of features (adjust as needed)
b_rg <- top_feat %>%
  # filter(name.x == "growthregrowth") %>% # if you have a variable for before/regrowth
  filter(!is.na(feature))  # just keep valid features

# Label a few genera if you have that info
# b_rg$label_it <- ifelse(b_rg$genus %in% c("Cupriavidus","Ralstonia","Burkholderia","Pseudomonas","Alcaligenes","Stenotrophomonas"), b_rg$genus, NA)

poster_colors = c("#dc9298","#a51c2f","#3d3d3d","#f3e3e5","#417c8c","#e57262","#915493")

poster_colors = c("#ce2e2d","#c00300","#88b8c5","#dd9298","#ab74a8","#fc3b3b","#a61d30","#eededf","#dc6c64","#467d91")


volc <- ggplot(b_rg, aes(
    x = coef,
    y = -log10(qval_joint),
    color = coef > 0,
    size = mean_abund
)) +
  geom_point() +
  # geom_text(aes(label = label_it), vjust = -1, size = 3, na.rm = TRUE, show.legend = FALSE, fontface = "italic") +
  scale_color_manual(values = c("#6d7ecd", "#c86c69"),labels = c("Biofilm","Wastewater")) +
  scale_size_continuous(range = c(1, 5)) +
    scale_y_continuous(
    breaks = c(0, 2.5, 5, 7.5, 10, 12.5, 15),
    limits = c(0, 15)
  ) +
  labs(
    title = "",
    x = "Effect Size (MaAsLin3 coefficient)",
    y = "-log10(q-value)",
    color = "",
    size = "Mean RPKM"
  ) +
  theme_minimal() +
  coord_cartesian(clip = "off")

volc

#get top ARGs-------

#label top args in new volcano plot----------
# Get top ARGs (e.g., top 5 by q-value)


top_args <- b_rg %>%
  arrange(qval_joint) %>%
  dplyr::slice(1:15) %>%        # top 15 most significant
  pull(feature)


top_args = b_rg %>% 
  #filter(-log10(qval_joint)>150) %>% 
  filter(-log10(qval_joint)>7.5 | (-log10(qval_joint)>4 & coef > 3|mean_abund>25000)) %>% 
  dplyr::select(feature) %>% distinct() %>% pull(); top_args


# Add a column for labels in the volcano data
b_rg$label_it <- ifelse(b_rg$feature %in% top_args, b_rg$feature, NA)

# Volcano plot with top ARG labels
volc <- ggplot(b_rg, aes(
    x = coef,
    y = -log10(qval_joint),
    color = coef > 0,
    size = mean_abund
)) +
  geom_point() +
    geom_text_repel(
    aes(label = label_it),
    size = 3,                 # smaller text
    max.overlaps = Inf,         # don’t drop labels
    box.padding = 0.3,          # space around text
    point.padding = 0.3,        # space from points
    min.segment.length = 0,     # always draw connecting line
    segment.size = 0.1,
    na.rm = TRUE,
    show.legend = FALSE,
    fontface = "italic"
  ) +
  #geom_text(aes(label = label_it), vjust = -1, size = 3, na.rm = TRUE, show.legend = FALSE, fontface = "italic") +
  scale_color_manual(values = c("#6d7ecd", "#c86c69"),labels = c("Biofilm","Wastewater")) +
  scale_size_continuous(range = c(1, 5)) +
    scale_y_continuous(
    breaks = c(0, 2.5, 5, 7.5, 10, 12.5, 15),
    limits = c(0, 15)
  ) +
  labs(
    title = "",
    x = "Effect Size (MaAsLin3 coefficient)",
    y = "-log10(q-value)",
    color = "",
    size = "Mean RPKM"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey88"),  # <--- major grid color
    panel.grid.minor = element_line(color = "grey88")  # <--- minor grid color
    #panel.border = element_rect(color = "", fill = NA)  # optional border
  ) +
  coord_cartesian(clip = "off")

volc

# Save final heatmap
ggsave(
  plot = volc,
  filename = "Biofilm Project Figures/CARD-MaAsLin3 Volcano Plot by ARG.png",
  width = 8,
  height = 4,
  units = "in",
  dpi = 300
)

# Save final heatmap
ggsave(
  plot = volc,
  filename = "Biofilm Project Figures/CARD-MaAsLin3 Volcano Plot by ARG.svg",
  width = 8,
  height = 4,
  units = "in",
  dpi = 500
)

# Get plot with all ARGs for visualizing---

top_args <- b_rg %>%
  arrange(qval_joint) %>%
  #slice(1:15) %>%        # top 15 most significant
  pull(feature)


# Add a column for labels in the volcano data
b_rg$label_it <- ifelse(b_rg$feature %in% top_args, b_rg$feature, NA)

# Volcano plot with top ARG labels
volc <- ggplot(b_rg, aes(
    x = coef,
    y = -log10(qval_joint),
    color = coef > 0,
    size = mean_abund
)) +
  geom_point() +
    geom_text_repel(
    aes(label = label_it),
    size = 3,                 # smaller text
    max.overlaps = Inf,         # don’t drop labels
    box.padding = 0.3,          # space around text
    point.padding = 0.3,        # space from points
    min.segment.length = 0,     # always draw connecting line
    segment.size = 0.1,
    na.rm = TRUE,
    show.legend = FALSE,
    fontface = "italic"
  ) +
  #geom_text(aes(label = label_it), vjust = -1, size = 3, na.rm = TRUE, show.legend = FALSE, fontface = "italic") +
  scale_color_manual(values = c("#6d7ecd", "#c86c69"),labels = c("Biofilm","Wastewater")) +
  scale_size_continuous(range = c(1, 5)) +
  scale_y_continuous(
    breaks = c(0, 2.5, 5, 7.5, 10, 12.5, 15),
    limits = c(0, 15)
  ) +
  labs(
    title = "",
    x = "Effect Size (MaAsLin3 coefficient)",
    y = "-log10(q-value)",
    color = "",
    size = "Mean Rel. Abundance"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey88"),  # <--- major grid color
    panel.grid.minor = element_line(color = "grey88")  # <--- minor grid color
    #panel.border = element_rect(color = "", fill = NA)  # optional border
  ) +
  coord_cartesian(clip = "off")

volc

# Save final heatmap
ggsave(
  plot = volc,
  filename = "Biofilm Project Figures/CARD-MaAsLin3 Volcano Plot by ARG All ARGs Labeled.png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)



#filtering for boxplots-----

sig_feats <- union(
  sig_abund_joint$feature,
  sig_prev_joint$feature
)

boxplot_df <- plot_df %>%
  filter(feature %in% sig_feats)

#join MaAsLin3 stats for box plots (effect sizes, q-values) -------
top_feat <- sig_abund_joint %>%
  filter(feature %in% sig_feats) %>%
  arrange(qval_joint)

boxplot_joined_df <- boxplot_df %>%
  left_join(
    top_feat %>% select(feature, coef, qval_joint),
    by = "feature"
  )


#loop through top args------------

for (arg_var in top_args) {


plot <- boxplot_joined_df %>% 
  filter(feature == arg_var) %>% 
  ggplot(., aes(x = Type, y = abundance, colour = Type)) +
  theme_minimal() +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(y = paste0(arg_var," Abundance"))

print(plot)

}

#create arg box plots using faceting approach

library(ggpubr)

plot <- boxplot_joined_df %>% 
  filter(feature %in% top_args) %>% 
  ggplot(aes(x = Type, y = abundance)) +
  
  geom_boxplot(
    outlier.shape = NA,
    color = "black",
    width = 0.6
  ) +
  
  geom_jitter(
    aes(color = Type),
    width = 0.15,
    size = 1.2,
    alpha = 0.6
  ) +
  
  stat_compare_means(
    method = "kruskal.test",
    label = "p.format",
    size = 3
  ) +
  
  facet_wrap(
    ~ feature,
    scales = "free_y",
    ncol = 4
  ) +
  
  labs(
    x = "",
    y = "Abundance"
  ) +
  
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 5),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "italic", size = 9),
    legend.position = "bottom"
  )

plot

# Save final heatmap
ggsave(
  plot = plot,
  filename = "Biofilm Project Figures/CARD-Box Plots Faceted of Top ARGs.png",
  width = 8,
  height = 12,
  units = "in",
  dpi = 300
)




#stats-------------
# 
# # Cupriavidus
# stat.test <- plot_df %>%
#   filter(genus == "Cupriavidus") %>% 
#   pairwise_wilcox_test(
#     abundance ~ growth,
#     p.adjust.method = "fdr")
# 
# c_plot <- cup_plot + stat_pvalue_manual(stat.test, label = "p.adj", y.position = c(0.55, 0.65, 0.75))


#MaAsLin3 table------

# mas_table <- test %>% 
#   select(name.x, qval_joint.x, mean_abund, genus) %>% 
#   filter(qval_joint.x > 0) %>% 
#   arrange(qval_joint.x) %>%  
#   slice(1:21)
# 
# write_xlsx(mas_table, "maaslin3_table.xlsx")



#------Extract All ARGs in Table Format by Mean Abundance, Coefficient, 
ontologies = card_merged_w_ont_df %>% 
  
  dplyr::select(Gene,`AMR Gene Family`,`Resistance Mechanism`,`card_edited_drug_class`,`AMR Gene Family`)

args = b_rg %>% 
  dplyr::select(feature, coef,N, mean_abund,qval_joint) %>% 
  
  mutate(enriched = ifelse(coef>0, "Wastewater", "Biofilm")) %>% 
  
  dplyr::rename("Gene" = "feature") %>% 
  
  left_join(ontologies,by = c("Gene")) %>% distinct()

#------Stacked Bar of What's Enriched in WW vs Biofilm----------------

stacked_bar =  args %>% 
  
  group_by(enriched) %>% 
  mutate(total = n()) %>% 
  ungroup() %>% 
  
  
  group_by(enriched,card_edited_drug_class) %>% 
  mutate(count = n()) %>% 
  dplyr::select(enriched,card_edited_drug_class,total,count) %>% distinct() %>% 
  summarise(percent = count/total * 100) %>% 
  ungroup() 
  

  
  
  
ggplot(stacked_bar, aes(x = factor(enriched), y = percent, fill = card_edited_drug_class)) +
  geom_col(position = "stack") 
  


```
 
 
## Maaslin3-ONLY SHARED DATES WITH BIOFILM- Drug Class
```{r}

#specify dates-------------
dates <- as.Date(c(
  "2024-09-09", "2024-09-11", "2024-09-12",
  "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09"
))

# Specify locations to include--------
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame



#libraries--------
library(ggplot2)
library(readxl)
library(here)
library(cowplot)
library(tidyverse)
library(patchwork)
library(readr)
library(phyloseq)
library(microshades)
library(ape)
library(vegan)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(maaslin3)
library(writexl)
library(ggrepel)

#maaslin3----------------
analysis_df <- card_merged_w_ont_df %>%
  filter(
    Type %in% c("Drain", "Endcap", "H_WW")
  ) %>%
  mutate(
    Date   = as.Date(Date),
    RPKM = as.numeric(RPKM),
    Type   = factor(Type, levels = c("H_WW","Endcap","Drain" ))
  ) %>%
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>%
  drop_na(timepoint)


biofilm_dates <- analysis_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
  mutate(
    biofilm_ref_date = Date
  ) %>% 
    ungroup() %>% 
  dplyr::distinct(Location,timepoint,biofilm_ref_date)

filtered_df <- analysis_df %>%
  mutate(
         Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  ) %>% 
  
  dplyr::select(-best_date,-date_diff_days,-biofilm_ref_date)

df <- card_merged_w_ont_df %>% 
  
  #clean ontology df
  dplyr::select(Gene,SampleID, Date, Location, Type,RPKM,`Drug Class`, `Resistance Mechanism`, class, mechanism,card_edited_drug_class, group) %>% 
  
  mutate(Date = as.Date(Date)) %>% 
  
  filter(Location %in% locations) %>% 
  
  #filter args in rpip list
  #filter(Species %in% c(amr_genes)) %>% 
  
  #filter for sample type
  filter(Type %in% c("Endcap", "H_WW")) %>% 
 
  #filter specific dates shared with biofilm
  filter(
    (Type == "H_WW" & Date %in% dates) |
    (Type != "H_WW")
  ) %>%
  
  #change label names
  mutate(
    Type = factor(Type, levels = c("Endcap", "H_WW")),
    TypeLabel = case_when(Type == "H_WW" ~ "Hospital Wastewater", 
                          TRUE ~ as.character(Type))
  ) 

df_sum <- df %>%
  group_by(SampleID, card_edited_drug_class) %>%
  summarise(
    RPKM = sum(RPKM, na.rm = TRUE),
    .groups = "drop"
  )

otu <- df_sum %>%
  dplyr::select(SampleID, card_edited_drug_class, RPKM) %>%
  tidyr::pivot_wider(
    names_from = card_edited_drug_class,
    values_from = RPKM,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("SampleID")

meta <- df %>%
  dplyr::select(
    SampleID,
    Date,
    Location,
    Type
  ) %>%
  distinct() %>%
  tibble::column_to_rownames("SampleID")

#MaAsLin3------------------------
fit_data <- maaslin3(
  input_data = otu,
  input_metadata = meta,
  output = "maaslin3_output",
  fixed_effects = c("Type"),
  random_effects = c("Location"),
  normalization = "NONE",
  transform = "LOG",
  standardize = TRUE,
  min_prevalence = 0.1,
  verbosity = "ERROR"
)



#plotting it-------------
plot_df <- fit_data$normalized %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "feature",
    values_to = "abundance"
  ) %>%
  left_join(
    fit_data$metadata %>%
      rownames_to_column("SampleID"),
    by = "SampleID"
  )


#volcano plot----------

#get mean abundance directly from MaAslin3 normalized data
# Significant features
sig_abund_joint <- subset(fit_data$fit_data_abundance$results, qval_joint < 0.05)
sig_prev_joint <- subset(fit_data$fit_data_prevalence$results, qval_joint < 0.05)



sig_feats <- union(sig_abund_joint$feature, sig_prev_joint$feature)

top_feat <- sig_abund_joint %>%
  filter(feature %in% sig_feats) %>%
  arrange(qval_joint)

# Compute mean abundance directly from MaAsLin3 normalized table
mean_abund <- colMeans(fit_data$normalized[, top_feat$feature, drop = FALSE], na.rm = TRUE)
top_feat$mean_abund <- mean_abund

# Optionally join extra info if you have a taxonomy table
# top_feat <- left_join(top_feat, top_feat_unique, by = "feature")

# Example: volcano for a subset of features (adjust as needed)
b_rg <- top_feat %>%
  # filter(name.x == "growthregrowth") %>% # if you have a variable for before/regrowth
  filter(!is.na(feature))  # just keep valid features

# Label a few genera if you have that info
# b_rg$label_it <- ifelse(b_rg$genus %in% c("Cupriavidus","Ralstonia","Burkholderia","Pseudomonas","Alcaligenes","Stenotrophomonas"), b_rg$genus, NA)

volc <- ggplot(b_rg, aes(
    x = coef,
    y = -log10(qval_joint),
    color = coef > 0,
    size = mean_abund
)) +
  geom_point() +
  # geom_text(aes(label = label_it), vjust = -1, size = 3, na.rm = TRUE, show.legend = FALSE, fontface = "italic") +
  scale_color_manual(values =c("#6d7ecd", "#c86c69")) +
  scale_size_continuous(range = c(1, 5)) +
  scale_y_continuous(
    breaks = c(0, 2.5, 5, 7.5),
    limits = c(0, 7.5)
  ) +
  labs(
    title = "",
    x = "Effect Size (MaAsLin3 coefficient)",
    y = "-log10(q-value)",
    color = "",
    size = "Mean RPKM"
  ) +
  theme_minimal() +
  coord_cartesian(clip = "off")

volc

#get top classes-------

#label top args in new volcano plot----------
# Get top ARGs (e.g., top 5 by q-value)

top_classes = b_rg %>% 
  #filter(-log10(qval_joint)>8 | (-log10(qval_joint)>5 & coef > 3) | mean_abund > 0.06) %>% 
  dplyr::select(feature) %>% distinct() %>% pull(); top_classes


# Add a column for labels in the volcano data
b_rg$label_it <- ifelse(b_rg$feature %in% top_classes, b_rg$feature, NA)

# Volcano plot with top ARG labels
volc <- ggplot(b_rg, aes(
    x = coef,
    y = -log10(qval_joint),
    color = coef > 0,
    size = mean_abund
)) +
  geom_point() +
    geom_text_repel(
    aes(label = label_it),
    size = 3,                 # smaller text
    max.overlaps = Inf,         # don’t drop labels
    box.padding = 0.3,          # space around text
    point.padding = 0.3,        # space from points
    min.segment.length = 0,     # always draw connecting line
    segment.size = 0.5,
    na.rm = TRUE,
    show.legend = FALSE,
    fontface = "italic"
  ) +
  #geom_text(aes(label = label_it), vjust = -1, size = 3, na.rm = TRUE, show.legend = FALSE, fontface = "italic") +
  scale_color_manual(values = c("#6d7ecd", "#c86c69"),labels = c("Biofilm","Wastewater")) +
  scale_size_continuous(range = c(1, 5)) +
  scale_y_continuous(
    breaks = c(0, 2.5, 5, 7.5),
    limits = c(0, 7.5)
  ) +
  labs(
    title = "",
    x = "Effect Size (MaAsLin3 coefficient)",
    y = "-log10(q-value)",
    color = "",
    size = "Mean RPKM"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey88"),  # <--- major grid color
    panel.grid.minor = element_line(color = "grey88")  # <--- minor grid color
    #panel.border = element_rect(color = "", fill = NA)  # optional border
  ) +
  coord_cartesian(clip = "off")

volc



# Save final heatmap
ggsave(
  plot = volc,
  filename = "Biofilm Project Figures/CARD-MaAsLin3 Volcano Plot by Drug Class.png",
  width = 8,
  height = 4,
  units = "in",
  dpi = 300
)

# Save final heatmap
ggsave(
  plot = volc,
  filename = "Biofilm Project Figures/CARD-MaAsLin3 Volcano Plot by Drug Class.svg",
  width = 8,
  height = 4,
  units = "in",
  dpi = 500
)




```
 
 
 #Time Series
##Maaslin3-ONLY SHARED DATES WITH BIOFILM- ARGs
```{r}

#specify dates-------------
dates <- as.Date(c(
  "2024-09-09", "2024-09-11", "2024-09-12",
  "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09"
))

# Specify locations to include--------
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame



#libraries--------
library(ggplot2)
library(readxl)
library(here)
library(cowplot)
library(tidyverse)
library(patchwork)
library(readr)
library(phyloseq)
library(microshades)
library(ape)
library(vegan)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(maaslin3)
library(writexl)
library(ggrepel)

#maaslin3----------------
analysis_df <- card_merged_w_ont_df %>%
  filter(
    Type %in% c("Drain", "Endcap", "H_WW")
  ) %>%
  mutate(
    Date   = as.Date(Date),
    RPKM = as.numeric(RPKM),
    Type   = factor(Type, levels = c("H_WW","Endcap","Drain" ))
  ) %>%
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>%
  drop_na(timepoint)


biofilm_dates <- analysis_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
  mutate(
    biofilm_ref_date = Date
  ) %>% 
    ungroup() %>% 
  dplyr::distinct(Location,timepoint,biofilm_ref_date)

filtered_df <- analysis_df %>%
  mutate(
         Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  ) %>% 
  
  dplyr::select(-best_date,-date_diff_days,-biofilm_ref_date)

df <- card_merged_w_ont_df %>% 
  
  #clean ontology df
  dplyr::select(Gene,SampleID, Date, Location, Type,RPKM,`Drug Class`, `Resistance Mechanism`, class, mechanism,card_edited_drug_class, group) %>% 
  
  mutate(Date = as.Date(Date)) %>% 
  
  filter(Location %in% locations) %>% 
  
  #filter args in rpip list
  #filter(Species %in% c(amr_genes)) %>% 
  
  #filter for sample type
  filter(Type %in% c("Endcap", "H_WW")) %>% 
 
  #filter specific dates shared with biofilm
  filter(
    (Type == "H_WW" & Date %in% dates) |
    (Type != "H_WW")
  ) %>%
  
  #change label names
  mutate(
    Type = factor(Type, levels = c("Endcap", "H_WW")),
    TypeLabel = case_when(Type == "H_WW" ~ "Hospital Wastewater", 
                          TRUE ~ as.character(Type))
  ) 

df_sum <- df %>%
  group_by(SampleID, Gene) %>%
  summarise(
    RPKM = sum(RPKM, na.rm = TRUE),
    .groups = "drop"
  )

otu <- df_sum %>%
  dplyr::select(SampleID, Gene, RPKM) %>%
  tidyr::pivot_wider(
    names_from = Gene,
    values_from = RPKM,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("SampleID")

meta <- df %>%
  dplyr::select(
    SampleID,
    Date,
    Location,
    Type
  ) %>%
  distinct() %>%
  tibble::column_to_rownames("SampleID")

#MaAsLin3------------------------

fit_data <- maaslin3(
  input_data = otu,
  input_metadata = meta,
  output = "maaslin3_output",
  fixed_effects = c("Type"),
  random_effects = c("Location"),
  normalization = "NONE",
  transform = "LOG",
  standardize = TRUE,
  min_prevalence = 0.1,
  verbosity = "ERROR"
)



#plotting it-------------
plot_df <- fit_data$normalized %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "feature",
    values_to = "abundance"
  ) %>%
  left_join(
    fit_data$metadata %>%
      rownames_to_column("SampleID"),
    by = "SampleID"
  )


#volcano plot----------

#get mean abundance directly from MaAslin3 normalized data
# Significant features
sig_abund_joint <- subset(fit_data$fit_data_abundance$results, qval_joint < 0.05)
sig_prev_joint <- subset(fit_data$fit_data_prevalence$results, qval_joint < 0.05)



sig_feats <- union(sig_abund_joint$feature, sig_prev_joint$feature)

top_feat <- sig_abund_joint %>%
  filter(feature %in% sig_feats) %>%
  arrange(qval_joint)

# Compute mean abundance directly from MaAsLin3 normalized table
mean_abund <- colMeans(fit_data$normalized[, top_feat$feature, drop = FALSE], na.rm = TRUE)
top_feat$mean_abund <- mean_abund

# Optionally join extra info if you have a taxonomy table
# top_feat <- left_join(top_feat, top_feat_unique, by = "feature")

# Example: volcano for a subset of features (adjust as needed)
b_rg <- top_feat %>%
  # filter(name.x == "growthregrowth") %>% # if you have a variable for before/regrowth
  filter(!is.na(feature))  # just keep valid features

# Label a few genera if you have that info
# b_rg$label_it <- ifelse(b_rg$genus %in% c("Cupriavidus","Ralstonia","Burkholderia","Pseudomonas","Alcaligenes","Stenotrophomonas"), b_rg$genus, NA)

poster_colors = c("#dc9298","#a51c2f","#3d3d3d","#f3e3e5","#417c8c","#e57262","#915493")

poster_colors = c("#ce2e2d","#c00300","#88b8c5","#dd9298","#ab74a8","#fc3b3b","#a61d30","#eededf","#dc6c64","#467d91")


volc <- ggplot(b_rg, aes(
    x = coef,
    y = -log10(qval_joint),
    color = coef > 0,
    size = mean_abund
)) +
  geom_point() +
  # geom_text(aes(label = label_it), vjust = -1, size = 3, na.rm = TRUE, show.legend = FALSE, fontface = "italic") +
  scale_color_manual(values = c("#6d7ecd", "#c86c69"),labels = c("Biofilm","Wastewater")) +
  scale_size_continuous(range = c(1, 5)) +
    scale_y_continuous(
    breaks = c(0, 2.5, 5, 7.5, 10, 12.5, 15),
    limits = c(0, 15)
  ) +
  labs(
    title = "",
    x = "Effect Size (MaAsLin3 coefficient)",
    y = "-log10(q-value)",
    color = "",
    size = "Mean RPKM"
  ) +
  theme_minimal() +
  coord_cartesian(clip = "off")

volc

#get top ARGs-------

#label top args in new volcano plot----------
# Get top ARGs (e.g., top 5 by q-value)


top_args <- b_rg %>%
  arrange(qval_joint) %>%
  dplyr::slice(1:15) %>%        # top 15 most significant
  pull(feature)


top_args = b_rg %>% 
  #filter(-log10(qval_joint)>150) %>% 
  filter(-log10(qval_joint)>7.5 | (-log10(qval_joint)>4 & coef > 3|mean_abund>25000)) %>% 
  dplyr::select(feature) %>% distinct() %>% pull(); top_args


# Add a column for labels in the volcano data
b_rg$label_it <- ifelse(b_rg$feature %in% top_args, b_rg$feature, NA)

# Volcano plot with top ARG labels
volc <- ggplot(b_rg, aes(
    x = coef,
    y = -log10(qval_joint),
    color = coef > 0,
    size = mean_abund
)) +
  geom_point() +
    geom_text_repel(
    aes(label = label_it),
    size = 3,                 # smaller text
    max.overlaps = Inf,         # don’t drop labels
    box.padding = 0.3,          # space around text
    point.padding = 0.3,        # space from points
    min.segment.length = 0,     # always draw connecting line
    segment.size = 0.1,
    na.rm = TRUE,
    show.legend = FALSE,
    fontface = "italic"
  ) +
  #geom_text(aes(label = label_it), vjust = -1, size = 3, na.rm = TRUE, show.legend = FALSE, fontface = "italic") +
  scale_color_manual(values = c("#6d7ecd", "#c86c69"),labels = c("Biofilm","Wastewater")) +
  scale_size_continuous(range = c(1, 5)) +
    scale_y_continuous(
    breaks = c(0, 2.5, 5, 7.5, 10, 12.5, 15),
    limits = c(0, 15)
  ) +
  labs(
    title = "",
    x = "Effect Size (MaAsLin3 coefficient)",
    y = "-log10(q-value)",
    color = "",
    size = "Mean RPKM"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey88"),  # <--- major grid color
    panel.grid.minor = element_line(color = "grey88")  # <--- minor grid color
    #panel.border = element_rect(color = "", fill = NA)  # optional border
  ) +
  coord_cartesian(clip = "off")

volc

# Save final heatmap
ggsave(
  plot = volc,
  filename = "Biofilm Project Figures/CARD-MaAsLin3 Volcano Plot by ARG.png",
  width = 8,
  height = 4,
  units = "in",
  dpi = 300
)

# Save final heatmap
ggsave(
  plot = volc,
  filename = "Biofilm Project Figures/CARD-MaAsLin3 Volcano Plot by ARG.svg",
  width = 8,
  height = 4,
  units = "in",
  dpi = 500
)

# Get plot with all ARGs for visualizing---

top_args <- b_rg %>%
  arrange(qval_joint) %>%
  #slice(1:15) %>%        # top 15 most significant
  pull(feature)


# Add a column for labels in the volcano data
b_rg$label_it <- ifelse(b_rg$feature %in% top_args, b_rg$feature, NA)

# Volcano plot with top ARG labels
volc <- ggplot(b_rg, aes(
    x = coef,
    y = -log10(qval_joint),
    color = coef > 0,
    size = mean_abund
)) +
  geom_point() +
    geom_text_repel(
    aes(label = label_it),
    size = 3,                 # smaller text
    max.overlaps = Inf,         # don’t drop labels
    box.padding = 0.3,          # space around text
    point.padding = 0.3,        # space from points
    min.segment.length = 0,     # always draw connecting line
    segment.size = 0.1,
    na.rm = TRUE,
    show.legend = FALSE,
    fontface = "italic"
  ) +
  #geom_text(aes(label = label_it), vjust = -1, size = 3, na.rm = TRUE, show.legend = FALSE, fontface = "italic") +
  scale_color_manual(values = c("#6d7ecd", "#c86c69"),labels = c("Biofilm","Wastewater")) +
  scale_size_continuous(range = c(1, 5)) +
  scale_y_continuous(
    breaks = c(0, 2.5, 5, 7.5, 10, 12.5, 15),
    limits = c(0, 15)
  ) +
  labs(
    title = "",
    x = "Effect Size (MaAsLin3 coefficient)",
    y = "-log10(q-value)",
    color = "",
    size = "Mean Rel. Abundance"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey88"),  # <--- major grid color
    panel.grid.minor = element_line(color = "grey88")  # <--- minor grid color
    #panel.border = element_rect(color = "", fill = NA)  # optional border
  ) +
  coord_cartesian(clip = "off")

volc

# Save final heatmap
ggsave(
  plot = volc,
  filename = "Biofilm Project Figures/CARD-MaAsLin3 Volcano Plot by ARG All ARGs Labeled.png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)



#filtering for boxplots-----

sig_feats <- union(
  sig_abund_joint$feature,
  sig_prev_joint$feature
)

boxplot_df <- plot_df %>%
  filter(feature %in% sig_feats)

#join MaAsLin3 stats for box plots (effect sizes, q-values) -------
top_feat <- sig_abund_joint %>%
  filter(feature %in% sig_feats) %>%
  arrange(qval_joint)

boxplot_joined_df <- boxplot_df %>%
  left_join(
    top_feat %>% select(feature, coef, qval_joint),
    by = "feature"
  )


#loop through top args------------

for (arg_var in top_args) {


plot <- boxplot_joined_df %>% 
  filter(feature == arg_var) %>% 
  ggplot(., aes(x = Type, y = abundance, colour = Type)) +
  theme_minimal() +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  labs(y = paste0(arg_var," Abundance"))

print(plot)

}

#create arg box plots using faceting approach

library(ggpubr)

plot <- boxplot_joined_df %>% 
  filter(feature %in% top_args) %>% 
  ggplot(aes(x = Type, y = abundance)) +
  
  geom_boxplot(
    outlier.shape = NA,
    color = "black",
    width = 0.6
  ) +
  
  geom_jitter(
    aes(color = Type),
    width = 0.15,
    size = 1.2,
    alpha = 0.6
  ) +
  
  stat_compare_means(
    method = "kruskal.test",
    label = "p.format",
    size = 3
  ) +
  
  facet_wrap(
    ~ feature,
    scales = "free_y",
    ncol = 4
  ) +
  
  labs(
    x = "",
    y = "Abundance"
  ) +
  
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 5),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "italic", size = 9),
    legend.position = "bottom"
  )

plot

# Save final heatmap
ggsave(
  plot = plot,
  filename = "Biofilm Project Figures/CARD-Box Plots Faceted of Top ARGs.png",
  width = 8,
  height = 12,
  units = "in",
  dpi = 300
)




#stats-------------
# 
# # Cupriavidus
# stat.test <- plot_df %>%
#   filter(genus == "Cupriavidus") %>% 
#   pairwise_wilcox_test(
#     abundance ~ growth,
#     p.adjust.method = "fdr")
# 
# c_plot <- cup_plot + stat_pvalue_manual(stat.test, label = "p.adj", y.position = c(0.55, 0.65, 0.75))


#MaAsLin3 table------

# mas_table <- test %>% 
#   select(name.x, qval_joint.x, mean_abund, genus) %>% 
#   filter(qval_joint.x > 0) %>% 
#   arrange(qval_joint.x) %>%  
#   slice(1:21)
# 
# write_xlsx(mas_table, "maaslin3_table.xlsx")



#------Extract All ARGs in Table Format by Mean Abundance, Coefficient, 
ontologies = card_merged_w_ont_df %>% 
  
  dplyr::select(Gene,`AMR Gene Family`,`Resistance Mechanism`,`card_edited_drug_class`,`AMR Gene Family`)

args = b_rg %>% 
  dplyr::select(feature, coef,N, mean_abund,qval_joint) %>% 
  
  mutate(enriched = ifelse(coef>0, "Wastewater", "Biofilm")) %>% 
  
  dplyr::rename("Gene" = "feature") %>% 
  
  left_join(ontologies,by = c("Gene")) %>% distinct()

#------Stacked Bar of What's Enriched in WW vs Biofilm----------------

stacked_bar =  args %>% 
  
  group_by(enriched) %>% 
  mutate(total = n()) %>% 
  ungroup() %>% 
  
  
  group_by(enriched,card_edited_drug_class) %>% 
  mutate(count = n()) %>% 
  dplyr::select(enriched,card_edited_drug_class,total,count) %>% distinct() %>% 
  summarise(percent = count/total * 100) %>% 
  ungroup() 
  

  
  
  
ggplot(stacked_bar, aes(x = factor(enriched), y = percent, fill = card_edited_drug_class)) +
  geom_col(position = "stack") 
  


```
 
 
## Maaslin3-ONLY SHARED DATES WITH BIOFILM- Drug Class
```{r}

#specify dates-------------
dates <- as.Date(c(
  "2024-09-09", "2024-09-11", "2024-09-12",
  "2024-10-30", "2024-10-31", "2024-12-04", "2024-12-09"
))

# Specify locations to include--------
locations <- c("FIONA", "SHREK", "OSCAR", "MARIO", "LUIGI")
# Assuming 'metaphlan_output' is your MetaPhlAn data frame



#libraries--------
library(ggplot2)
library(readxl)
library(here)
library(cowplot)
library(tidyverse)
library(patchwork)
library(readr)
library(phyloseq)
library(microshades)
library(ape)
library(vegan)
library(ggpubr)
library(rstatix)
library(RColorBrewer)
library(maaslin3)
library(writexl)
library(ggrepel)

#maaslin3----------------
analysis_df <- card_merged_w_ont_df %>%
  filter(
    Type %in% c("Drain", "Endcap", "H_WW")
  ) %>%
  mutate(
    Date   = as.Date(Date),
    RPKM = as.numeric(RPKM),
    Type   = factor(Type, levels = c("H_WW","Endcap","Drain" ))
  ) %>%
  left_join(
    sixteens_df,
    by = c("Type", "Location", "Date")
  ) %>%
  drop_na(timepoint)


biofilm_dates <- analysis_df %>%
  filter(Type %in% c("Endcap")) %>% 
  group_by(Location, timepoint) %>%
  mutate(
    biofilm_ref_date = Date
  ) %>% 
    ungroup() %>% 
  dplyr::distinct(Location,timepoint,biofilm_ref_date)

filtered_df <- analysis_df %>%
  mutate(
         Date = as.Date(Date)) %>%
  left_join(biofilm_dates, by = c("Location", "timepoint")) %>% 

  mutate(
    date_diff_days = abs(as.numeric(biofilm_ref_date-Date))
  ) %>% 

  #match temporal subsampling
  group_by(Location, timepoint, Type) %>%
  mutate(
    best_date = min(date_diff_days)
  ) %>%
  ungroup() %>% 

  filter(
    (best_date == date_diff_days)
  ) %>% 
  
  dplyr::select(-best_date,-date_diff_days,-biofilm_ref_date)

df <- card_merged_w_ont_df %>% 
  
  #clean ontology df
  dplyr::select(Gene,SampleID, Date, Location, Type,RPKM,`Drug Class`, `Resistance Mechanism`, class, mechanism,card_edited_drug_class, group) %>% 
  
  mutate(Date = as.Date(Date)) %>% 
  
  filter(Location %in% locations) %>% 
  
  #filter args in rpip list
  #filter(Species %in% c(amr_genes)) %>% 
  
  #filter for sample type
  filter(Type %in% c("Endcap", "H_WW")) %>% 
 
  #filter specific dates shared with biofilm
  filter(
    (Type == "H_WW" & Date %in% dates) |
    (Type != "H_WW")
  ) %>%
  
  #change label names
  mutate(
    Type = factor(Type, levels = c("Endcap", "H_WW")),
    TypeLabel = case_when(Type == "H_WW" ~ "Hospital Wastewater", 
                          TRUE ~ as.character(Type))
  ) 

df_sum <- df %>%
  group_by(SampleID, card_edited_drug_class) %>%
  summarise(
    RPKM = sum(RPKM, na.rm = TRUE),
    .groups = "drop"
  )

otu <- df_sum %>%
  dplyr::select(SampleID, card_edited_drug_class, RPKM) %>%
  tidyr::pivot_wider(
    names_from = card_edited_drug_class,
    values_from = RPKM,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("SampleID")

meta <- df %>%
  dplyr::select(
    SampleID,
    Date,
    Location,
    Type
  ) %>%
  distinct() %>%
  tibble::column_to_rownames("SampleID")

#MaAsLin3------------------------
fit_data <- maaslin3(
  input_data = otu,
  input_metadata = meta,
  output = "maaslin3_output",
  fixed_effects = c("Type"),
  random_effects = c("Location"),
  normalization = "NONE",
  transform = "LOG",
  standardize = TRUE,
  min_prevalence = 0.1,
  verbosity = "ERROR"
)



#plotting it-------------
plot_df <- fit_data$normalized %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(
    cols = -SampleID,
    names_to = "feature",
    values_to = "abundance"
  ) %>%
  left_join(
    fit_data$metadata %>%
      rownames_to_column("SampleID"),
    by = "SampleID"
  )


#volcano plot----------

#get mean abundance directly from MaAslin3 normalized data
# Significant features
sig_abund_joint <- subset(fit_data$fit_data_abundance$results, qval_joint < 0.05)
sig_prev_joint <- subset(fit_data$fit_data_prevalence$results, qval_joint < 0.05)



sig_feats <- union(sig_abund_joint$feature, sig_prev_joint$feature)

top_feat <- sig_abund_joint %>%
  filter(feature %in% sig_feats) %>%
  arrange(qval_joint)

# Compute mean abundance directly from MaAsLin3 normalized table
mean_abund <- colMeans(fit_data$normalized[, top_feat$feature, drop = FALSE], na.rm = TRUE)
top_feat$mean_abund <- mean_abund

# Optionally join extra info if you have a taxonomy table
# top_feat <- left_join(top_feat, top_feat_unique, by = "feature")

# Example: volcano for a subset of features (adjust as needed)
b_rg <- top_feat %>%
  # filter(name.x == "growthregrowth") %>% # if you have a variable for before/regrowth
  filter(!is.na(feature))  # just keep valid features

# Label a few genera if you have that info
# b_rg$label_it <- ifelse(b_rg$genus %in% c("Cupriavidus","Ralstonia","Burkholderia","Pseudomonas","Alcaligenes","Stenotrophomonas"), b_rg$genus, NA)

volc <- ggplot(b_rg, aes(
    x = coef,
    y = -log10(qval_joint),
    color = coef > 0,
    size = mean_abund
)) +
  geom_point() +
  # geom_text(aes(label = label_it), vjust = -1, size = 3, na.rm = TRUE, show.legend = FALSE, fontface = "italic") +
  scale_color_manual(values =c("#6d7ecd", "#c86c69")) +
  scale_size_continuous(range = c(1, 5)) +
  scale_y_continuous(
    breaks = c(0, 2.5, 5, 7.5),
    limits = c(0, 7.5)
  ) +
  labs(
    title = "",
    x = "Effect Size (MaAsLin3 coefficient)",
    y = "-log10(q-value)",
    color = "",
    size = "Mean RPKM"
  ) +
  theme_minimal() +
  coord_cartesian(clip = "off")

volc

#get top classes-------

#label top args in new volcano plot----------
# Get top ARGs (e.g., top 5 by q-value)

top_classes = b_rg %>% 
  #filter(-log10(qval_joint)>8 | (-log10(qval_joint)>5 & coef > 3) | mean_abund > 0.06) %>% 
  dplyr::select(feature) %>% distinct() %>% pull(); top_classes


# Add a column for labels in the volcano data
b_rg$label_it <- ifelse(b_rg$feature %in% top_classes, b_rg$feature, NA)

# Volcano plot with top ARG labels
volc <- ggplot(b_rg, aes(
    x = coef,
    y = -log10(qval_joint),
    color = coef > 0,
    size = mean_abund
)) +
  geom_point() +
    geom_text_repel(
    aes(label = label_it),
    size = 3,                 # smaller text
    max.overlaps = Inf,         # don’t drop labels
    box.padding = 0.3,          # space around text
    point.padding = 0.3,        # space from points
    min.segment.length = 0,     # always draw connecting line
    segment.size = 0.5,
    na.rm = TRUE,
    show.legend = FALSE,
    fontface = "italic"
  ) +
  #geom_text(aes(label = label_it), vjust = -1, size = 3, na.rm = TRUE, show.legend = FALSE, fontface = "italic") +
  scale_color_manual(values = c("#6d7ecd", "#c86c69"),labels = c("Biofilm","Wastewater")) +
  scale_size_continuous(range = c(1, 5)) +
  scale_y_continuous(
    breaks = c(0, 2.5, 5, 7.5),
    limits = c(0, 7.5)
  ) +
  labs(
    title = "",
    x = "Effect Size (MaAsLin3 coefficient)",
    y = "-log10(q-value)",
    color = "",
    size = "Mean RPKM"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey88"),  # <--- major grid color
    panel.grid.minor = element_line(color = "grey88")  # <--- minor grid color
    #panel.border = element_rect(color = "", fill = NA)  # optional border
  ) +
  coord_cartesian(clip = "off")

volc



# Save final heatmap
ggsave(
  plot = volc,
  filename = "Biofilm Project Figures/CARD-MaAsLin3 Volcano Plot by Drug Class.png",
  width = 8,
  height = 4,
  units = "in",
  dpi = 300
)

# Save final heatmap
ggsave(
  plot = volc,
  filename = "Biofilm Project Figures/CARD-MaAsLin3 Volcano Plot by Drug Class.svg",
  width = 8,
  height = 4,
  units = "in",
  dpi = 500
)




```
 
 
 #Time Series
