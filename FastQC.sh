#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-09:00
#SBATCH -p shared,sapphire
#SBATCH --mem=4G
#SBATCH -o /n/home11/avdarling/slurm/%j.FastQC.out.output
#SBATCH -e /n/home11/avdarling/slurm/%j.FastQC.out.err

source /n/home11/avdarling/miniconda3/etc/profile.d/conda.sh
conda activate /n/netscratch/hhealy_lab/avdarling_conda_envs/fastqc_env

cd /n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/Ginkgo_rpip || exit

SAMPLE="BANANA"  # replaced in per-sample script
base="${SAMPLE}"

echo "Base: ${base}"

input_R1=$(find . -type f -name "${SAMPLE}_R1.fastq.gz" | head -n 1)
input_R2=$(find . -type f -name "${SAMPLE}_R2.fastq.gz" | head -n 1)

if [[ -z "$input_R1" || -z "$input_R2" ]]; then
    echo "Could not find input files for ${SAMPLE}"
    exit 1
fi

echo "Found R1: $input_R1"
echo "Found R2: $input_R2"

if gunzip -t "$input_R1" 2>/dev/null; then
    echo "$input_R1 is in correct GZIP format."
else
    echo "$input_R1 is NOT in correct GZIP format or is corrupted."
    exit 1
fi

if gunzip -t "$input_R2" 2>/dev/null; then
    echo "$input_R2 is in correct GZIP format."
else
    echo "$input_R2 is NOT in correct GZIP format or is corrupted."
    exit 1
fi

echo "Running FastQC"

OUTDIR="/n/netscratch/hhealy_lab/avdarling/fastqc"
mkdir -p "$OUTDIR"

if [[ -e "${OUTDIR}/${base}_R1_fastqc.html" && -e "${OUTDIR}/${base}_R1_fastqc.zip" ]]; then
    echo "FastQC already completed for R1"
else
    fastqc "$input_R1" -o "$OUTDIR"
fi

if [[ -e "${OUTDIR}/${base}_R2_fastqc.html" && -e "${OUTDIR}/${base}_R2_fastqc.zip" ]]; then
    echo "FastQC already completed for R2"
else
    fastqc "$input_R2" -o "$OUTDIR"
fi
