#!/bin/bash

### Using the base script "metaphlan.sh", write a script for each sample
#substitute the word banana


# Set the base directory where all subfolders are located
cd /n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/Ginkgo_rpip || exit

# Find all R1 files recursively that match the pattern
find . -type f -name "*_R1.fastq.gz" | while read -r filepath; do
    # Extract the filename from the path
    filename=$(basename "$filepath" _R1.fastq.gz)
    sample="$filename"

    echo "Sample: $sample"

    # Create the merge script for the sample by replacing BANANA in the template
    sed 's@BANANA@'"${sample}"'@' /n/home11/avdarling/scripts/fullWorkflow.sh > "/n/home11/avdarling/scripts/workflow_scripts_per_sample/${sample}.full_workflow.sh"
done
