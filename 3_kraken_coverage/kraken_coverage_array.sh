#!/bin/bash
#SBATCH -J kraken_cov
#SBATCH -c 8
#SBATCH -t 3-00:00
#SBATCH -p hsph,shared,sapphire
#SBATCH --mem=100G
#SBATCH --array=0-4540%20
#SBATCH -o /n/home11/avdarling/slurm/%A_%a.kraken_coverage.output
#SBATCH -e /n/home11/avdarling/slurm/%A_%a.kraken_coverage.err

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
COMBINATIONS="/n/home11/avdarling/downloadableData/kraken_coverage_combinations_hospital.txt"

KRAKEN_DIR="/n/netscratch/hhealy_lab/avdarling/kraken_out/kraken_output_ct0_5_min_hit_3"
FASTQ_DIR="/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/Ginkgo_rpip_fastqs"
REF_BASE="/n/netscratch/hhealy_lab/avdarling/reference_genomes"
OUT_BASE="/n/netscratch/hhealy_lab/avdarling/kraken_coverage"

PYTHON="/n/home11/avdarling/miniconda3/bin/python"
BWA="/n/holylabs/hhealy_lab/Lab/avdarling_conda_envs/bwa_env/bin/bwa"
SAMTOOLS="/n/holylabs/hhealy_lab/Lab/avdarling_conda_envs/bwa_env/bin/samtools"
MOSDEPTH="/n/holylabs/hhealy_lab/Lab/avdarling_conda_envs/bwa_env/bin/mosdepth"
KRAKENTOOLS="/n/home11/avdarling/KrakenTools/extract_kraken_reads.py"

THREADS=8
# ──────────────────────────────────────────────────────────────────────────────

# Read this task's sample + taxid from the combinations file
COMBO=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$COMBINATIONS")
SAMPLE=$(echo "$COMBO" | awk '{print $1}')
TAXID=$(echo "$COMBO"  | awk '{print $2}')

if [[ -z "$SAMPLE" || -z "$TAXID" ]]; then
    echo "ERROR: empty SAMPLE or TAXID for array index ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

# ── Input paths ────────────────────────────────────────────────────────────────
KRAKEN="${KRAKEN_DIR}/${SAMPLE}.kraken"
REPORT="${KRAKEN_DIR}/${SAMPLE}.kreport"
R1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"
REF="${REF_BASE}/taxid_${TAXID}/reference_genome.fna"

# ── Output paths ───────────────────────────────────────────────────────────────
OUTDIR="${OUT_BASE}/${SAMPLE}_${TAXID}"
mkdir -p "$OUTDIR"

if [[ -f "${OUTDIR}/coverage_metrics.tsv" ]]; then
    echo "Already complete, skipping"
    exit 0
fi

cd "$OUTDIR"

EXTRACT_R1="${OUTDIR}/${SAMPLE}_${TAXID}_R1.fastq"
EXTRACT_R2="${OUTDIR}/${SAMPLE}_${TAXID}_R2.fastq"

echo "=========================================="
echo "Array task : ${SLURM_ARRAY_TASK_ID}"
echo "Sample     : ${SAMPLE}"
echo "TaxID      : ${TAXID}"
echo "=========================================="
date

# ── Sanity-check required files ────────────────────────────────────────────────
for f in "$PYTHON" "$BWA" "$SAMTOOLS" "$MOSDEPTH" "$KRAKEN" "$REPORT" "$R1" "$R2" "$REF" "$KRAKENTOOLS"; do
    if [[ ! -e "$f" ]]; then
        echo "ERROR: required file not found: $f" >&2
        exit 1
    fi
done

# ── Extract Kraken-classified reads ───────────────────────────────────────────
"$PYTHON" "$KRAKENTOOLS" \
    -k "$KRAKEN" \
    -r "$REPORT" \
    -s1 "$R1" \
    -s2 "$R2" \
    -t "$TAXID" \
    --include-children \
    --fastq-output \
    -o "$EXTRACT_R1" \
    -o2 "$EXTRACT_R2"

# ── Handle zero-read case ─────────────────────────────────────────────────────
if [[ ! -s "$EXTRACT_R1" || ! -s "$EXTRACT_R2" ]]; then
    echo "No extracted reads for sample=${SAMPLE}, taxid=${TAXID}" | tee no_reads.txt
    cat > coverage_metrics.tsv <<EOF
sample	taxid	extracted_r1_reads	extracted_r2_reads	mapped_reads	reference_length	breadth_1x	breadth_10x	mean_depth	median_depth
${SAMPLE}	${TAXID}	0	0	0	NA	0	0	0	0
EOF
    exit 0
fi

EXTRACTED_R1_READS=$(awk 'END{print NR/4}' "$EXTRACT_R1")
EXTRACTED_R2_READS=$(awk 'END{print NR/4}' "$EXTRACT_R2")
echo "Extracted R1 reads: $EXTRACTED_R1_READS"
echo "Extracted R2 reads: $EXTRACTED_R2_READS"

# ── Build BWA index if needed ─────────────────────────────────────────────────
if [[ ! -f "${REF}.bwt" ]]; then
    echo "Building BWA index..."
    "$BWA" index "$REF"
fi

# ── Align with BWA ────────────────────────────────────────────────────────────
echo "Aligning extracted reads with BWA..."
"$BWA" mem -t "$THREADS" "$REF" "$EXTRACT_R1" "$EXTRACT_R2" 2> bwa.log \
    | "$SAMTOOLS" view -@ "$THREADS" -b - \
    | "$SAMTOOLS" sort -@ "$THREADS" -o mapped.sorted.bam

"$SAMTOOLS" index mapped.sorted.bam

MAPPED_READS=$("$SAMTOOLS" view -c -F 4 mapped.sorted.bam)

# ── Run mosdepth for coverage metrics ─────────────────────────────────────────
# --quantize 0:1:10: : create coverage bands (0x, 1-9x, 10x+) for breadth
echo "Running mosdepth..."
"$MOSDEPTH" \
    --threads "$THREADS" \
    --quantize 0:1:10: \
    "${OUTDIR}/mosdepth" \
    mapped.sorted.bam

# ── Extract metrics from mosdepth output ──────────────────────────────────────
# Summary file columns: chrom, length, bases, mean, min, max
MEAN_DEPTH=$(awk '$1=="total" {print $4}' "${OUTDIR}/mosdepth.mosdepth.summary.txt")
REFERENCE_LENGTH=$(awk '$1=="total" {print $2}' "${OUTDIR}/mosdepth.mosdepth.summary.txt")

# Breadth: fraction of reference covered at >=1x and >=10x
# Quantized bed columns: chrom, start, end, band (e.g. "1:10" or "10:inf")
BREADTH_1X=$(zcat "${OUTDIR}/mosdepth.quantized.bed.gz" | awk -v reflen="$REFERENCE_LENGTH" '
    $4 != "0:1" { sum += ($3 - $2) }
    END { printf "%.6f", (reflen>0 ? sum/reflen : 0) }')

BREADTH_10X=$(zcat "${OUTDIR}/mosdepth.quantized.bed.gz" | awk -v reflen="$REFERENCE_LENGTH" '
    $4 == "10:inf" { sum += ($3 - $2) }
    END { printf "%.6f", (reflen>0 ? sum/reflen : 0) }')

# Median: computed from per-base depth (efficient for pathogen-sized genomes)
# Note: for sparse metagenomic data, median is often 0 — breadth_1x is more informative
MEDIAN_DEPTH=$(zcat "${OUTDIR}/mosdepth.per-base.bed.gz" \
    | awk '{print $4}' \
    | sort -n \
    | awk 'BEGIN{c=0} {a[c++]=$1} END{
        if (c==0) { print 0 }
        else if (c%2) { print a[int(c/2)] }
        else { print (a[c/2-1]+a[c/2])/2 }
    }')

# ── Write coverage metrics ─────────────────────────────────────────────────────
cat > coverage_metrics.tsv <<EOF
sample	taxid	extracted_r1_reads	extracted_r2_reads	mapped_reads	reference_length	breadth_1x	breadth_10x	mean_depth	median_depth
${SAMPLE}	${TAXID}	${EXTRACTED_R1_READS}	${EXTRACTED_R2_READS}	${MAPPED_READS}	${REFERENCE_LENGTH}	${BREADTH_1X}	${BREADTH_10X}	${MEAN_DEPTH}	${MEDIAN_DEPTH}
EOF

echo "Done."
date
echo "Outputs in: $OUTDIR"
echo "  coverage_metrics.tsv"
echo "  mosdepth.mosdepth.summary.txt"
echo "  mosdepth.per-base.bed.gz"
echo "  mosdepth.quantized.bed.gz"
echo "  bwa.log"
