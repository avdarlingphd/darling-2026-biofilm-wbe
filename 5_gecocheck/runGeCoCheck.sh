#!/bin/bash
#SBATCH -c 8
#SBATCH -t 7-00:00
#SBATCH -p intermediate
#SBATCH --mem=100G
#SBATCH -o /n/home11/avdarling/slurm/%j.gecocheck.output
#SBATCH -e /n/home11/avdarling/slurm/%j.gecocheck.err

set -euo pipefail
source /n/home11/avdarling/miniconda3/etc/profile.d/conda.sh
conda activate /n/holylabs/hhealy_lab/Lab/avdarling_conda_envs/coverage_checker
cd /n/netscratch/hhealy_lab/avdarling/gcc_run

DEST=/n/holylabs/hhealy_lab/Lab/ynhh_ww_rpip_2024

mkdir -p /n/netscratch/hhealy_lab/avdarling/gcc_run/Ginkgo_rpip_coverage_results_v10

coverage_pipeline.py \
    --processors 8 \
    --sample_metadata /n/home11/avdarling/sample_metadata_nonmunicipal.csv \
    --project_name Ginkgo_rpip \
    --fastq_dir /n/netscratch/hhealy_lab/avdarling/gcc_run/gcc_merged_fastq \
    --kraken_kreport_dir /n/netscratch/hhealy_lab/avdarling/kraken_out/kraken_output_ct0_5_min_hit_3 \
    --kraken_outraw_dir /n/netscratch/hhealy_lab/avdarling/kraken_out/kraken_output_ct0_5_min_hit_3 \
    --output_dir /n/netscratch/hhealy_lab/avdarling/gcc_run/Ginkgo_rpip_coverage_results_v10 \
    --genome_dir ${DEST}/gcc_genome_cache \
    --bowtie2_db_dir ${DEST}/gcc_bowtie2_cache \
    --coverage_program Bowtie2 \
    --no_grouped_samples \
    --read_lim 50000 \
    --skip_cleanup
