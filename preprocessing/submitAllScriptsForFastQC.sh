#!/bin/bash

cd /n/home11/avdarling/scripts/FastQCscripts

# Loop through all script files in the directory
for script_file in *.sh; do
    if [ -f "$script_file" ]; then
        echo "Submitting $script_file to Slurm..."
        sbatch "$script_file"
    fi
done
