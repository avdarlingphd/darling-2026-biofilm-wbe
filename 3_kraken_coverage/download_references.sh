#!/bin/bash
#SBATCH -J download_refs
#SBATCH -c 1
#SBATCH -t 4:00:00
#SBATCH -p shared
#SBATCH --mem=4G
#SBATCH -o /n/home11/avdarling/slurm/%j.download_refs.out
#SBATCH -e /n/home11/avdarling/slurm/%j.download_refs.err

# Download NCBI reference genomes for all taxids in taxids.txt
# Usage: sbatch download_references.sh
# Or interactively: bash download_references.sh
#
# Reads taxids from: /n/home11/avdarling/taxids.txt  (one taxid per line)
# Saves each genome to: /n/netscratch/hhealy_lab/avdarling/reference_genomes/taxid_XXXXX/reference_genome.fna

set -euo pipefail

TAXIDS_FILE="/n/home11/avdarling/downloadableData/new_taxids_to_download.txt"
REF_BASE="/n/netscratch/hhealy_lab/avdarling/reference_genomes"

# Activate conda to get the datasets CLI
# (datasets is a compiled binary — works even though Python in this env is broken)
source /n/home11/avdarling/miniconda3/etc/profile.d/conda.sh
conda activate /n/netscratch/hhealy_lab/avdarling_conda_envs/kraken2_2.17

if [[ ! -f "$TAXIDS_FILE" ]]; then
    echo "ERROR: taxids file not found: $TAXIDS_FILE" >&2
    echo "Create a file with one taxid per line, e.g.:" >&2
    echo "  13373" >&2
    echo "  2697049" >&2
    exit 1
fi

echo "Downloading reference genomes for taxids in: $TAXIDS_FILE"
echo ""

while IFS= read -r taxid || [[ -n "$taxid" ]]; do
    # Skip blank lines and comments
    [[ -z "$taxid" || "$taxid" == \#* ]] && continue

    REFDIR="${REF_BASE}/taxid_${taxid}"
    REF="${REFDIR}/reference_genome.fna"

    if [[ -f "$REF" ]]; then
        echo "[SKIP] taxid ${taxid} — reference already exists"
        continue
    fi

    echo "[DOWNLOAD] taxid ${taxid}..."
    mkdir -p "$REFDIR"
    cd "$REFDIR"

    datasets download genome taxon "$taxid" \
        --reference \
        --include genome \
        --filename genome.zip || {
            echo "  WARNING: download failed for taxid ${taxid}, skipping" >&2
            cd - > /dev/null
            continue
        }

    unzip genome.zip -d . > /dev/null

    # Use cat instead of mv — handles genomes with multiple .fna files (e.g. multiple chromosomes)
    FNA_FILES=(ncbi_dataset/data/*/*_genomic.fna)
    if [[ ${#FNA_FILES[@]} -eq 0 ]]; then
        echo "  WARNING: no *_genomic.fna found for taxid ${taxid}, skipping" >&2
        cd - > /dev/null
        continue
    fi
    cat "${FNA_FILES[@]}" > reference_genome.fna

    # Clean up zip and temp files
    rm -rf genome.zip ncbi_dataset README.md

    echo "  -> Saved: $REF"
    cd - > /dev/null

done < "$TAXIDS_FILE"

echo ""
echo "All done. Reference genomes in: $REF_BASE"
