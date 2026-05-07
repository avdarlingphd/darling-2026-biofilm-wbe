#!/usr/bin/env python3
"""
Generate per-sample/taxid SLURM scripts for BLASTN validation.

One job per sample/taxid pair. Jobs gracefully skip if FASTA files
are missing (no reads extracted) or BLAST results already exist.

Usage:
    python makeScripts_blast.py \
        --taxid_list sample_taxid_list_all_pathogens.tsv \
        --base_script blast_validation_base.sh \
        --outdir /path/to/blast_scripts \
        --slurm_logdir /path/to/slurm/logs

Then submit with:
    bash /path/to/blast_scripts/submitAll_blast.sh
"""

import os
import argparse


def main():
    parser = argparse.ArgumentParser(
        description="Generate per-sample/taxid SLURM BLAST scripts"
    )
    parser.add_argument(
        "--taxid_list", required=True,
        help="TSV with columns: sample, taxid (no header)"
    )
    parser.add_argument(
        "--base_script", required=True,
        help="Path to blast_validation_base.sh template"
    )
    parser.add_argument(
        "--outdir", required=True,
        help="Directory to write generated scripts"
    )
    parser.add_argument(
        "--slurm_logdir", required=True,
        help="Directory for SLURM .out/.err log files"
    )
    args = parser.parse_args()

    # --- Read base template ---
    with open(args.base_script) as f:
        template = f.read()

    # --- Read sample/taxid pairs ---
    pairs = []
    with open(args.taxid_list) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            pairs.append((parts[0], parts[1]))

    os.makedirs(args.outdir, exist_ok=True)

    script_paths = []
    for sample, taxid in pairs:
        script = template.replace("BANANA_SAMPLE", sample)
        script = script.replace("BANANA_TAXID", taxid)
        script = script.replace("BANANA_SLURM_LOGDIR", args.slurm_logdir)

        script_path = os.path.join(args.outdir, f"blast_{sample}_{taxid}.sh")
        with open(script_path, "w") as f:
            f.write(script)
        script_paths.append(script_path)

    # --- Write submitAll script ---
    submit_path = os.path.join(args.outdir, "submitAll_blast.sh")
    with open(submit_path, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("# Submit all BLAST validation jobs\n")
        f.write(f"# Generated for {len(script_paths)} sample/taxid pairs\n\n")
        for path in script_paths:
            f.write(f"sbatch {path}\n")
        f.write('\necho "Submitted all BLAST jobs."\n')

    print(f"Generated {len(script_paths)} scripts in: {args.outdir}")
    print(f"Submit script: {submit_path}")
    print(f"\nTo submit all jobs:")
    print(f"  bash {submit_path}")


if __name__ == "__main__":
    main()
