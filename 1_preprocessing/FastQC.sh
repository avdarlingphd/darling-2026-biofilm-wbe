#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-09:00
#SBATCH -p shared,sapphire
#SBATCH --mem=4G
# ---- SLURM log paths (edit if needed; cannot use shell variables) ----
#SBATCH -o /n/home11/avdarling/slurm/%j.FastQC.out.output
#SBATCH -e /n/home11/avdarling/slurm/%j.FastQC.out.err

set -euo pipefail

# ---------------- Config (edit these for your environment) ----------------
INPUT_DIR="/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/Ginkgo_rpip"
OUTDIR="/n/netscratch/hhealy_lab/avdarling/fastqc"
CONDA_SH="/n/home11/avdarling/miniconda3/etc/profile.d/conda.sh"
FASTQC_ENV="/n/netscratch/hhealy_lab/avdarling_conda_envs/fastqc_env"
# --------------------------------------------------------------------------

# ---------------- Sample name (replaced per sample by makeScriptsforFastQC.sh) ----
SAMPLE="BANANA"
base="${SAMPLE}"
echo "Base: ${base}"

# ---------------- Activate conda env ----------------
source "$CONDA_SH"
conda activate "$FASTQC_ENV"

# ---------------- Locate inputs ----------------
cd "$INPUT_DIR" || exit 1

input_R1=$(find . -type f -name "${SAMPLE}_R1.fastq.gz" | head -n 1)
input_R2=$(find . -type f -name "${SAMPLE}_R2.fastq.gz" | head -n 1)

if [[ -z "$input_R1" || -z "$input_R2" ]]; then
    echo "ERROR: Could not find input files for ${SAMPLE}"
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

# ---------------- Run FastQC ----------------
mkdir -p "$OUTDIR"
echo "Running FastQC"

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
