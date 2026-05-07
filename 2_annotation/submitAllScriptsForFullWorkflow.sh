#!/bin/bash
set -euo pipefail

# ---------------- Config (edit for your environment) ----------------
SCRIPT_DIR="/n/home11/avdarling/scripts/workflow_scripts_per_sample"
# --------------------------------------------------------------------

cd "$SCRIPT_DIR" || exit 1

# Loop through all script files in the directory
for script_file in *.sh; do
    if [ -f "$script_file" ]; then
        echo "Submitting $script_file to Slurm..."
        sbatch "$script_file"
    fi
done
