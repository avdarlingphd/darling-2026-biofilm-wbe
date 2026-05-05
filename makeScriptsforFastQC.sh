#!/bin/bash

cd /n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024/Ginkgo_rpip || exit

mkdir -p /n/home11/avdarling/scripts/FastQCscripts

find . -type f -name "*_R1.fastq.gz" | while read -r filepath; do
    sample=$(basename "$filepath" _R1.fastq.gz)

    echo "Sample: $sample"

    sed "s@BANANA@${sample}@g" \
        /n/home11/avdarling/scripts/FastQC.sh \
        > "/n/home11/avdarling/scripts/FastQCscripts/${sample}.fastqc.sh"
done

