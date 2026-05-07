#!/usr/bin/env python3
"""
Generate per-sample SLURM scripts for Kraken2 read extraction.

Each job processes ONE sample and ALL its taxids in a single pass over the
.kraken file — much faster than one job per sample/taxid pair.

Usage:
    python makeScripts_extract.py \
        --taxid_list sample_taxid_list_all_pathogens.tsv \
        --base_script extract_kraken_base.sh \
        --outdir /path/to/extract_scripts \
        --slurm_logdir /path/to/slurm/logs

Then submit with:
    bash /path/to/extract_scripts/submitAll_extract.sh
"""

import os
import argparse
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(
        description="Generate per-sample SLURM extraction scripts"
    )
    parser.add_argument(
        "--taxid_list", required=True,
        help="TSV with columns: sample, taxid (no header)"
    )
    parser.add_argument(
        "--base_script", required=True,
        help="Path to extract_kraken_base.sh template"
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

    # --- Group taxids by sample ---
    sample_taxids = defaultdict(list)
    with open(args.taxid_list) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            sample, taxid = parts[0], parts[1]
            sample_taxids[sample].append(taxid)

    os.makedirs(args.outdir, exist_ok=True)

    script_paths = []
    for sample, taxids in sorted(sample_taxids.items()):
        taxid_array = " ".join(taxids)  # bash array elements: "570 573 562 ..."
        script = template.replace("BANANA_SAMPLE", sample)
        script = script.replace("BANANA_TAXIDS", taxid_array)
        script = script.replace("BANANA_SLURM_LOGDIR", args.slurm_logdir)

        script_path = os.path.join(args.outdir, f"extract_{sample}.sh")
        with open(script_path, "w") as f:
            f.write(script)
        script_paths.append(script_path)

    # --- Write submitAll script ---
    submit_path = os.path.join(args.outdir, "submitAll_extract.sh")
    with open(submit_path, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("# Submit all per-sample extraction jobs\n")
        f.write(f"# Generated for {len(script_paths)} samples\n\n")
        for path in script_paths:
            f.write(f"sbatch {path}\n")
        f.write('\necho "Submitted all extraction jobs."\n')

    print(f"Generated {len(script_paths)} scripts in: {args.outdir}")
    print(f"Submit script: {submit_path}")
    print(f"\nTo submit all jobs:")
    print(f"  bash {submit_path}")
    print(f"\nSample breakdown:")
    counts = sorted([(s, len(t)) for s, t in sample_taxids.items()],
                    key=lambda x: -x[1])
    print(f"  Max taxids per sample:  {counts[0][1]}  ({counts[0][0]})")
    print(f"  Min taxids per sample:  {counts[-1][1]}  ({counts[-1][0]})")
    total = sum(len(t) for t in sample_taxids.values())
    print(f"  Total sample/taxid pairs: {total}")


if __name__ == "__main__":
    main()
