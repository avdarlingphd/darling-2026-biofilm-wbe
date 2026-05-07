#!/bin/bash
#SBATCH -c 8
#SBATCH -t 1-00:00
#SBATCH -p hsph,shared,sapphire
#SBATCH --mem=32G
#SBATCH -o /n/home11/avdarling/slurm/%j.blast_BANANA_SAMPLE_BANANA_TAXID.out
#SBATCH -e /n/home11/avdarling/slurm/%j.blast_BANANA_SAMPLE_BANANA_TAXID.err
set -euo pipefail

# Use full paths — avoids conda activation issues in non-interactive SLURM shells
BLASTN="/n/home11/avdarling/conda_envs/blast_env/bin/blastn"
SEQTK="/n/netscratch/hhealy_lab/avdarling_conda_envs/seqtk_env/bin/seqtk"
export BLASTDB="/n/netscratch/informatics/Everyone/external_repos/blast/nt/latest/"

# Sanity check: confirm binaries exist
if [ ! -x "${BLASTN}" ]; then
    echo "ERROR: blastn not found at ${BLASTN}"
    exit 1
fi
if [ ! -x "${SEQTK}" ]; then
    echo "ERROR: seqtk not found at ${SEQTK}"
    exit 1
fi

SAMPLE="BANANA_SAMPLE"
TAXID="BANANA_TAXID"
OUTDIR="/n/netscratch/hhealy_lab/avdarling/blast_validation/blast_validation_${SAMPLE}_${TAXID}"
DB="/n/netscratch/informatics/Everyone/external_repos/blast/nt/latest/nt"
R1_FASTA="${OUTDIR}/${SAMPLE}_taxid_${TAXID}_R1.fasta"
R2_FASTA="${OUTDIR}/${SAMPLE}_taxid_${TAXID}_R2.fasta"
R1_FASTA_SUB="${OUTDIR}/${SAMPLE}_taxid_${TAXID}_R1_sub.fasta"
R2_FASTA_SUB="${OUTDIR}/${SAMPLE}_taxid_${TAXID}_R2_sub.fasta"
R1_OUT="${OUTDIR}/${SAMPLE}_taxid_${TAXID}_R1_blast_results.tsv"
R2_OUT="${OUTDIR}/${SAMPLE}_taxid_${TAXID}_R2_blast_results.tsv"

MAX_READS=500  # cap per read set — plenty for false positive detection

echo "Sample: ${SAMPLE}, TaxID: ${TAXID}"
echo "Using BLAST database: ${DB}"
echo "blastn: ${BLASTN}"

if [ ! -f "${DB}.nal" ]; then
    echo "ERROR: BLAST database not found at ${DB}"
    exit 1
fi

# --- Skip if FASTA files are missing or empty ---
if [ ! -s "${R1_FASTA}" ] || [ ! -s "${R2_FASTA}" ]; then
    echo "SKIPPING: FASTA files missing or empty for ${SAMPLE} taxid ${TAXID}"
    echo "  Expected: ${R1_FASTA}"
    echo "  Expected: ${R2_FASTA}"
    exit 0
fi

# --- Subsample to MAX_READS to keep BLAST runtime manageable ---
N_R1=$(grep -c "^>" "${R1_FASTA}" || true)
echo "R1 total reads: ${N_R1}"

if [ "${N_R1}" -gt "${MAX_READS}" ]; then
    echo "Subsampling R1 to ${MAX_READS} reads..."
    "${SEQTK}" sample -s 42 "${R1_FASTA}" "${MAX_READS}" > "${R1_FASTA_SUB}"
    "${SEQTK}" sample -s 42 "${R2_FASTA}" "${MAX_READS}" > "${R2_FASTA_SUB}"
    R1_QUERY="${R1_FASTA_SUB}"
    R2_QUERY="${R2_FASTA_SUB}"
else
    echo "Read count within limit — using full FASTA"
    R1_QUERY="${R1_FASTA}"
    R2_QUERY="${R2_FASTA}"
fi

# --- BLAST R1 ---
if [ -s "${R1_OUT}" ]; then
    echo "R1 BLAST results already exist. Skipping..."
else
    echo "Running BLAST on R1 ($(grep -c "^>" "${R1_QUERY}") reads)..."
    "${BLASTN}" \
        -query "${R1_QUERY}" \
        -db "${DB}" \
        -task blastn \
        -num_threads 8 \
        -max_target_seqs 20 \
        -outfmt "6 qseqid sscinames pident length evalue stitle" \
        -out "${R1_OUT}"
    echo "R1 BLAST done."
fi

# --- BLAST R2 ---
if [ -s "${R2_OUT}" ]; then
    echo "R2 BLAST results already exist. Skipping..."
else
    echo "Running BLAST on R2 ($(grep -c "^>" "${R2_QUERY}") reads)..."
    "${BLASTN}" \
        -query "${R2_QUERY}" \
        -db "${DB}" \
        -task blastn \
        -num_threads 8 \
        -max_target_seqs 20 \
        -outfmt "6 qseqid sscinames pident length evalue stitle" \
        -out "${R2_OUT}"
    echo "R2 BLAST done."
fi

echo "All done for ${SAMPLE} taxid ${TAXID}"
