#!/bin/bash
# Aggregate all per-sample/taxid coverage_metrics.tsv files into one table
# Usage: bash aggregate_results.sh
#
# Output: all_coverage_metrics.tsv in the current directory

set -euo pipefail

RESULTS_BASE="/n/netscratch/hhealy_lab/avdarling/kraken_coverage"
OUTPUT="all_coverage_metrics.tsv"

# Write header from the first file found
HEADER_WRITTEN=0
> "$OUTPUT"

for tsv in "${RESULTS_BASE}"/*/*/coverage_metrics.tsv "${RESULTS_BASE}"/*/coverage_metrics.tsv; do
    [[ -f "$tsv" ]] || continue

    if [[ $HEADER_WRITTEN -eq 0 ]]; then
        head -1 "$tsv" > "$OUTPUT"
        HEADER_WRITTEN=1
    fi

    # Append data rows only (skip header line)
    tail -n +2 "$tsv" >> "$OUTPUT"
done

NROWS=$(( $(wc -l < "$OUTPUT") - 1 ))
echo "Aggregated $NROWS rows into: $OUTPUT"
