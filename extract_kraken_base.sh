#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-04:00
#SBATCH -p hsph,shared,sapphire
#SBATCH --mem=16G
#SBATCH -o /n/home11/avdarling/slurm/%j.extract_BANANA_SAMPLE.out
#SBATCH -e /n/home11/avdarling/slurm/%j.extract_BANANA_SAMPLE.err
set -euo pipefail
source /n/home11/avdarling/miniconda3/etc/profile.d/conda.sh
conda activate /n/netscratch/hhealy_lab/avdarling_conda_envs/seqtk_env
echo "CONDA_PREFIX=${CONDA_PREFIX}"
which seqtk

KRAKENTOOLS="/n/home11/avdarling/KrakenTools/extract_kraken_reads.py"
PYTHON="/n/home11/avdarling/miniconda3/bin/python3"
if [ ! -f "${KRAKENTOOLS}" ]; then
    echo "ERROR: extract_kraken_reads.py not found at ${KRAKENTOOLS}"
    exit 1
fi
if [ ! -f "${PYTHON}" ]; then
    echo "ERROR: base python3 not found at ${PYTHON}"
    exit 1
fi

SAMPLE="BANANA_SAMPLE"
TAXIDS=(BANANA_TAXIDS)   # array of taxids for this sample

KRAKEN_BASE="/n/netscratch/hhealy_lab/avdarling/kraken_out/kraken_output_ct0_5_min_hit_3"
FASTQ_BASE="/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/Ginkgo_rpip_fastqs"
OUTBASE="/n/netscratch/hhealy_lab/avdarling/blast_validation"

KRAKEN_FILE="${KRAKEN_BASE}/${SAMPLE}.kraken"
KREPORT_FILE="${KRAKEN_BASE}/${SAMPLE}.kreport"
FASTQ_R1="${FASTQ_BASE}/${SAMPLE}_R1.fastq.gz"
FASTQ_R2="${FASTQ_BASE}/${SAMPLE}_R2.fastq.gz"

# --- Validate shared inputs ---
for f in "${KRAKEN_FILE}" "${KREPORT_FILE}" "${FASTQ_R1}" "${FASTQ_R2}"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required file not found: $f"
        exit 1
    fi
done

echo "Processing sample: ${SAMPLE}"
echo "Taxids: ${TAXIDS[*]}"

for TAXID in "${TAXIDS[@]}"; do
    echo "========================================"
    echo "  TAXID=${TAXID}"

    OUTDIR="${OUTBASE}/blast_validation_${SAMPLE}_${TAXID}"
    R1_OUT_FASTQ="${OUTDIR}/${SAMPLE}_taxid_${TAXID}_R1.fastq"
    R2_OUT_FASTQ="${OUTDIR}/${SAMPLE}_taxid_${TAXID}_R2.fastq"
    R1_OUT_FASTA="${OUTDIR}/${SAMPLE}_taxid_${TAXID}_R1.fasta"
    R2_OUT_FASTA="${OUTDIR}/${SAMPLE}_taxid_${TAXID}_R2.fasta"
    SUMMARY="${OUTDIR}/${SAMPLE}_taxid_${TAXID}_summary.txt"

    # --- SKIP if output FASTAs already exist and are non-empty ---
    if [ -s "${R1_OUT_FASTA}" ] && [ -s "${R2_OUT_FASTA}" ]; then
        echo "  SKIPPING: FASTAs already exist for taxid ${TAXID}"
        continue
    fi

    mkdir -p "${OUTDIR}"

    # --- Extract reads using KrakenTools (includes child taxon reads) ---
    # extract_kraken_reads.py uses the kreport to walk the taxonomy tree,
    # so reads classified to any strain/subspecies under TAXID are captured.
    echo "  Extracting paired reads for taxid ${TAXID} (including children)"
    "${PYTHON}" "${KRAKENTOOLS}" \
        -k "${KRAKEN_FILE}" \
        -r "${KREPORT_FILE}" \
        -s  "${FASTQ_R1}" \
        -s2 "${FASTQ_R2}" \
        -t  "${TAXID}" \
        --include-children \
        --fastq-output \
        -o  "${R1_OUT_FASTQ}" \
        -o2 "${R2_OUT_FASTQ}"

    # --- Count extracted reads ---
    NREADS_R1=0
    NREADS_R2=0
    if [ -s "${R1_OUT_FASTQ}" ]; then
        NREADS_R1=$(grep -c "^@" "${R1_OUT_FASTQ}" || true)
    fi
    if [ -s "${R2_OUT_FASTQ}" ]; then
        NREADS_R2=$(grep -c "^@" "${R2_OUT_FASTQ}" || true)
    fi
    echo "  R1 extracted reads: ${NREADS_R1}"
    echo "  R2 extracted reads: ${NREADS_R2}"

    {
        echo "sample=${SAMPLE}"
        echo "taxid=${TAXID}"
        echo "nreads_R1=${NREADS_R1}"
        echo "nreads_R2=${NREADS_R2}"
    } > "${SUMMARY}"

    if [ "${NREADS_R1}" -eq 0 ]; then
        echo "  No reads found for taxid ${TAXID}; skipping FASTA conversion."
        continue
    fi

    # --- Convert FASTQ → FASTA (seqtk is fast for this step) ---
    echo "  Converting R1 to FASTA"
    seqtk seq -a "${R1_OUT_FASTQ}" > "${R1_OUT_FASTA}"
    echo "  Converting R2 to FASTA"
    seqtk seq -a "${R2_OUT_FASTQ}" > "${R2_OUT_FASTA}"

    echo "  Done: ${OUTDIR}"
done

echo "========================================"
echo "All taxids processed for ${SAMPLE}."
