#!/bin/bash
#SBATCH -c 16
#SBATCH -t 0-06:00          # 6 hours should be plenty for Kraken2 + Bracken
#SBATCH -p hsph,sapphire,shared
#SBATCH --mem=128G           # pluspf database needs ~74GB so keep this
#SBATCH -o /n/home11/avdarling/slurm/%j.kraken_bracken.output
#SBATCH -e /n/home11/avdarling/slurm/%j.kraken_bracken.err

set -euo pipefail

# Enable conda in this shell script
source /n/home11/avdarling/miniconda3/etc/profile.d/conda.sh



# ---------------- Sample name ----------------
SAMPLE="BANANA"  # This will be replaced in the per-sample script
base=$(echo ${SAMPLE} | sed "s/_R1.fastq.gz//")

cd /n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/Ginkgo_rpip || exit 1
echo "Processing sample: ${base}"

# ---------------- Find FASTQ files ----------------
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
source /n/home11/avdarling/miniconda3/etc/profile.d/conda.sh
conda activate /n/holylabs/hhealy_lab/Lab/avdarling_conda_envs/kraken2_latest
kraken_out_dir4="/n/netscratch/hhealy_lab/avdarling/kraken_out/kraken_output_ct0_5_min_hit_3"

mkdir -p "$kraken_out_dir4"
kraken_output4="${kraken_out_dir4}/${base}.kraken2"
kraken_report4="${kraken_out_dir4}/${base}.k2report"

if [ -s "$kraken_output4" ]; then
    echo "Kraken2-recreated output exists. Skipping..."
else
    echo "Running Kraken2-recreated for ${base}..."
    kraken2 \
        --db /n/holylabs/hhealy_lab/Lab/databases/kraken2_pluspf \
        --threads 16 \
        --confidence 0.5 \
        --minimum-hit-groups 3 \
        --paired "$input_R1" "$input_R2" \
        --report "$kraken_report4" \
        --output "$kraken_output4"
fi




# ---------------- CARD ----------------
source /n/home11/avdarling/miniconda3/etc/profile.d/conda.sh​
conda activate diamond_env​

card_out_dir="/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/CARD_output"​
mkdir -p "$card_out_dir"​
card_output="${card_out_dir}/${base}.card"​

if [ -f "$card_output" ]; then​
echo "CARD output exists. Skipping..."​
else​
echo "Running CARD for ${base}..."​
temp_fastq="temp_${base}_combined.fastq.gz"​

# Combine R1 and R2 into one gzipped FASTQ​
echo "Combining R1 and R2 into $temp_fastq ..."​
cat "$input_R1" "$input_R2" > "$temp_fastq"
diamond blastx \​
-d /n/home11/avdarling/databases/CARD/card_db \​
-q "$temp_fastq" \​
-o "$card_output" \​
--max-target-seqs 1 \​
--evalue 1e-10 \​
--id 80​
​
# Clean up temporary combined file​
echo "Cleaning up temporary file: $temp_fastq"​
rm -f "$temp_fastq"​
fi

# ---------------- Centrifuger ----------------
source /n/home11/avdarling/miniconda3/etc/profile.d/conda.sh
conda activate /n/holylabs/hhealy_lab/Lab/avdarling_conda_envs/centrifuger_env

centrifuger_out_dir="/n/netscratch/hhealy_lab/avdarling/kraken_out/centrifuger_output"
mkdir -p "$centrifuger_out_dir"

centrifuger_output="${centrifuger_out_dir}/${base}.centrifuger"

if [ -s "$centrifuger_output" ]; then
    echo "Centrifuger output exists. Skipping..."
else
    echo "Running Centrifuger for ${base}..."

    centrifuger \
        -x /n/holylabs/hhealy_lab/Lab/databases/centrifuger_db/cfr_hpv+gbsarscov2 \
        -1 "$input_R1" -2 "$input_R2" \
        -t 16 \
        > "$centrifuger_output"
fi

centrifuger_quant_output="${centrifuger_out_dir}/${base}.centrifuger_quant.tsv"

if [ -f "$centrifuger_quant_output" ]; then
    echo "Centrifuger quantification output exists. Skipping..."
else
    echo "Running Centrifuger quantification for ${base}..."
    
    centrifuger-quant \
        -x /n/holylabs/hhealy_lab/Lab/databases/centrifuger_db/cfr_hpv+gbsarscov2 \
        -c "$centrifuger_output" \
        --output-format 1 \
        > "$centrifuger_quant_output"
fi

    
#-----------------Sylph-----------------------
source /n/home11/avdarling/miniconda3/etc/profile.d/conda.sh
conda activate /n/holylabs/hhealy_lab/Lab/avdarling_conda_envs/sylph_env

sylph_out_dir="/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/sylph_output"
mkdir -p "$sylph_out_dir"

sylph_file="${sylph_out_dir}/${base}_results.tsv"

if [ -f "$sylph_file" ]; then
    echo "Sylph output exists for ${base}. Skipping..."
else
    echo "Running Sylph for ${base}..."
    sylph profile /n/holylabs/hhealy_lab/Lab/databases/gtdb-r220-c200-dbv1.syldb \
        -1 "$input_R1" \
        -2 "$input_R2" \
        -t 10 \
        -o "$sylph_file"
fi


sylph_taxonomy_dir="/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/sylph_output/sylph_taxonomy_out"
mkdir -p "$sylph_taxonomy_dir"

# Use ONLY the sample name, nothing else
sylph_tax_prefix="${sylph_taxonomy_dir}/${base}"

# Check for ANY file Sylph might create
if ls "${sylph_tax_prefix}"*.sylphmpa 1> /dev/null 2>&1; then
    echo "Sylph tax file exists for ${base}. Skipping..."
else
    echo "Running Sylph taxonomy for ${base}..."

    sylph-tax taxprof "$sylph_file" \
        -t GTDB_r220 \
        -o "${sylph_tax_prefix}"
fi


# ---------------- MetaPhlAn ----------------
source /n/home11/avdarling/miniconda3/etc/profile.d/conda.sh
conda activate metaphlan_env

output_dir="/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/metaphlan_output_parallel_nonanadama_rerun"
mkdir -p "$output_dir"
mkdir -p "$output_dir/sams" "$output_dir/bowtie2" "$output_dir/profiles"

profile_file="${output_dir}/profiles/${base}_profile.txt"
bowtie_file="${output_dir}/bowtie2/${base}.bowtie2.bz2"
sam_file="${output_dir}/sams/${base}.sam.bz2"

# Skip MetaPhlAn if the profile already exists
if [ -f "$profile_file" ]; then
    echo "MetaPhlAn output exists for ${base}. Skipping..."
else
    db_dir="/n/home11/avdarling/miniconda3/envs/metaphlan_env/lib/python3.10/site-packages/metaphlan/metaphlan_databases"

    metaphlan "$input_R1,$input_R2" \
        --input_type fastq \
        --db_dir "$db_dir" \
        --tax_lev a \
        --nproc 16 \
        --mapout "$bowtie_file" \
        -s "$sam_file" \
        -o "$profile_file"

    echo "Full MetaPhlAn run finished for ${base}."
fi
