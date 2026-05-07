#!/bin/bash
############ CARD alignment, Kraken2, MetaPhlAn, Sylph, Centrifuger ############
#SBATCH -c 16
#SBATCH -t 0-06:00          # 6 hours
#SBATCH -p hsph,sapphire,shared
#SBATCH --mem=128G           # Kraken2 PlusPF DB needs ~74 GB
# ---- SLURM log paths (edit if needed; cannot use shell variables) ----
#SBATCH -o /n/home11/avdarling/slurm/%j.fullWorkflow.output
#SBATCH -e /n/home11/avdarling/slurm/%j.fullWorkflow.err

set -euo pipefail

# ---------------- Config (edit these for your environment) ----------------
# Conda
CONDA_SH="/n/home11/avdarling/miniconda3/etc/profile.d/conda.sh"
KRAKEN2_ENV="/n/holylabs/hhealy_lab/Lab/avdarling_conda_envs/kraken2_latest"
DIAMOND_ENV="diamond_env"
CENTRIFUGER_ENV="/n/holylabs/hhealy_lab/Lab/avdarling_conda_envs/centrifuger_env"
SYLPH_ENV="/n/holylabs/hhealy_lab/Lab/avdarling_conda_envs/sylph_env"
METAPHLAN_ENV="metaphlan_env"

# Input
INPUT_DIR="/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/Ginkgo_rpip"

# Reference databases
KRAKEN2_DB="/n/holylabs/hhealy_lab/Lab/databases/kraken2_pluspf"
CARD_DB="/n/home11/avdarling/databases/CARD/card_db"
CENTRIFUGER_DB="/n/holylabs/hhealy_lab/Lab/databases/centrifuger_db/cfr_hpv+gbsarscov2"
SYLPH_DB="/n/holylabs/hhealy_lab/Lab/databases/gtdb-r220-c200-dbv1.syldb"
SYLPH_TAX="GTDB_r220"
METAPHLAN_DB_DIR="/n/home11/avdarling/miniconda3/envs/metaphlan_env/lib/python3.10/site-packages/metaphlan/metaphlan_databases"

# Output directories
KRAKEN_OUT_DIR="/n/netscratch/hhealy_lab/avdarling/kraken_out/kraken_output_ct0_5_min_hit_3"
CARD_OUT_DIR="/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/CARD_output"
CENTRIFUGER_OUT_DIR="/n/netscratch/hhealy_lab/avdarling/kraken_out/centrifuger_output"
SYLPH_OUT_DIR="/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/sylph_output"
SYLPH_TAX_DIR="${SYLPH_OUT_DIR}/sylph_taxonomy_out"
METAPHLAN_OUT_DIR="/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/metaphlan_output_parallel_nonanadama_rerun"

# Tool parameters
THREADS=16
KRAKEN2_CONFIDENCE=0.5
KRAKEN2_MIN_HIT_GROUPS=3
CARD_EVALUE="1e-10"
CARD_ID=80
# --------------------------------------------------------------------------

# Enable conda in this shell
source "$CONDA_SH"

# ---------------- Sample name (replaced per sample by makeScriptsFullWorkflow.sh) ----
SAMPLE="BANANA"
base=$(echo "${SAMPLE}" | sed "s/_R1.fastq.gz//")
echo "Processing sample: ${base}"

cd "$INPUT_DIR" || exit 1

# ---------------- Locate inputs ----------------
input_R1=$(find . -type f -name "${SAMPLE}_R1.fastq.gz" | head -n 1)
input_R2=$(find . -type f -name "${SAMPLE}_R2.fastq.gz" | head -n 1)
if [[ -z "$input_R1" || -z "$input_R2" ]]; then
    echo "ERROR: Could not find FASTQ files for ${SAMPLE}"
    exit 1
fi
echo "Found R1: $input_R1"
echo "Found R2: $input_R2"

# ---------------- Validate gzip ----------------
for f in "$input_R1" "$input_R2"; do
    if gunzip -t "$f" 2>/dev/null; then
        echo "$f is valid GZIP."
    else
        echo "ERROR: $f is not valid GZIP or is corrupted."
        exit 1
    fi
done

# ---------------- Kraken2 ----------------
conda activate "$KRAKEN2_ENV"
mkdir -p "$KRAKEN_OUT_DIR"
kraken_output="${KRAKEN_OUT_DIR}/${base}.kraken2"
kraken_report="${KRAKEN_OUT_DIR}/${base}.k2report"

if [ -s "$kraken_output" ]; then
    echo "Kraken2 output exists. Skipping..."
else
    echo "Running Kraken2 for ${base}..."
    kraken2 \
        --db "$KRAKEN2_DB" \
        --threads "$THREADS" \
        --confidence "$KRAKEN2_CONFIDENCE" \
        --minimum-hit-groups "$KRAKEN2_MIN_HIT_GROUPS" \
        --paired "$input_R1" "$input_R2" \
        --report "$kraken_report" \
        --output "$kraken_output"
fi
conda deactivate

# ---------------- CARD (DIAMOND blastx) ----------------
conda activate "$DIAMOND_ENV"
mkdir -p "$CARD_OUT_DIR"
card_output="${CARD_OUT_DIR}/${base}.card"

if [ -f "$card_output" ]; then
    echo "CARD output exists. Skipping..."
else
    echo "Running CARD for ${base}..."
    temp_fastq="temp_${base}_combined.fastq.gz"
    echo "Combining R1 and R2 into $temp_fastq ..."
    cat "$input_R1" "$input_R2" > "$temp_fastq"

    diamond blastx \
        -d "$CARD_DB" \
        -q "$temp_fastq" \
        -o "$card_output" \
        --max-target-seqs 1 \
        --evalue "$CARD_EVALUE" \
        --id "$CARD_ID"

    echo "Cleaning up temporary file: $temp_fastq"
    rm -f "$temp_fastq"
fi
conda deactivate

# ---------------- Centrifuger ----------------
conda activate "$CENTRIFUGER_ENV"
mkdir -p "$CENTRIFUGER_OUT_DIR"
centrifuger_output="${CENTRIFUGER_OUT_DIR}/${base}.centrifuger"
centrifuger_quant_output="${CENTRIFUGER_OUT_DIR}/${base}.centrifuger_quant.tsv"

if [ -s "$centrifuger_output" ]; then
    echo "Centrifuger output exists. Skipping..."
else
    echo "Running Centrifuger for ${base}..."
    centrifuger \
        -x "$CENTRIFUGER_DB" \
        -1 "$input_R1" -2 "$input_R2" \
        -t "$THREADS" \
        > "$centrifuger_output"
fi

if [ -f "$centrifuger_quant_output" ]; then
    echo "Centrifuger quantification output exists. Skipping..."
else
    echo "Running Centrifuger quantification for ${base}..."
    centrifuger-quant \
        -x "$CENTRIFUGER_DB" \
        -c "$centrifuger_output" \
        --output-format 1 \
        > "$centrifuger_quant_output"
fi
conda deactivate

# ---------------- Sylph ----------------
conda activate "$SYLPH_ENV"
mkdir -p "$SYLPH_OUT_DIR"
sylph_file="${SYLPH_OUT_DIR}/${base}_results.tsv"

if [ -f "$sylph_file" ]; then
    echo "Sylph output exists for ${base}. Skipping..."
else
    echo "Running Sylph for ${base}..."
    sylph profile "$SYLPH_DB" \
        -1 "$input_R1" \
        -2 "$input_R2" \
        -t 10 \
        -o "$sylph_file"
fi

mkdir -p "$SYLPH_TAX_DIR"
sylph_tax_prefix="${SYLPH_TAX_DIR}/${base}"

if ls "${sylph_tax_prefix}"*.sylphmpa 1> /dev/null 2>&1; then
    echo "Sylph tax file exists for ${base}. Skipping..."
else
    echo "Running Sylph taxonomy for ${base}..."
    sylph-tax taxprof "$sylph_file" \
        -t "$SYLPH_TAX" \
        -o "${sylph_tax_prefix}"
fi
conda deactivate

# ---------------- MetaPhlAn ----------------
conda activate "$METAPHLAN_ENV"
mkdir -p "$METAPHLAN_OUT_DIR/sams" "$METAPHLAN_OUT_DIR/bowtie2" "$METAPHLAN_OUT_DIR/profiles"
profile_file="${METAPHLAN_OUT_DIR}/profiles/${base}_profile.txt"
bowtie_file="${METAPHLAN_OUT_DIR}/bowtie2/${base}.bowtie2.bz2"
sam_file="${METAPHLAN_OUT_DIR}/sams/${base}.sam.bz2"

if [ -f "$profile_file" ]; then
    echo "MetaPhlAn output exists for ${base}. Skipping..."
else
    echo "Running MetaPhlAn for ${base}..."
    metaphlan "$input_R1,$input_R2" \
        --input_type fastq \
        --db_dir "$METAPHLAN_DB_DIR" \
        --tax_lev a \
        --nproc "$THREADS" \
        --mapout "$bowtie_file" \
        -s "$sam_file" \
        -o "$profile_file"
    echo "MetaPhlAn finished for ${base}."
fi
conda deactivate
