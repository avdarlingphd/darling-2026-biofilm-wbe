#!/usr/bin/env python3
"""
BLAST False Positive Detection Script
======================================
Analyzes BLASTN results to classify Kraken2 hits as TRUE_POSITIVE,
FALSE_POSITIVE, UNCERTAIN, UNCULTURED_DOMINANT, or BLAST_NOT_RUN
for each sample/taxid pair.

Logic:
  - For each read, consider ALL 20 BLAST hits (not just the top hit).
    Count how many reads had any hit to each organism — this gives a
    per-organism read count across all hits returned.
  - pct_match = % of reads that had ANY hit matching the expected organism
    across all 20 hits.
  - Uncultured check uses the best (top) hit per read.
  - Genus-level fallback matching only used when expected organism is
    genus-only (one word). Species-level expected organisms require a
    species-level hit (prevents "Naegleria sp." matching "Naegleria fowleri").

TRUE_POSITIVE requires ALL THREE:
  1. avg pct_match >= threshold (default 80%) across R1 + R2
  2. Expected organism is ranked #1 or #2 by read count in BOTH R1 and R2
  3. n_match (reads matching expected organism) > 2 in BOTH R1 and R2

Other classifications:
  - UNCULTURED_DOMINANT: >50% of best hits are uncultured/environmental
  - FALSE_POSITIVE: avg pct_match < 10%
  - UNCERTAIN: everything else

Usage:
    python false_positive_detection.py \
        --blast_dir /n/netscratch/hhealy_lab/avdarling/blast_validation \
        --taxid_list /n/home11/avdarling/downloadableData/sample_taxid_list_all_pathogens.tsv \
        --taxid_names /n/home11/avdarling/downloadableData/taxid_to_name_all_pathogens.csv \
        --output /n/home11/avdarling/downloadableData/blast_false_positive_report.tsv

Requirements:
    pip install pandas
"""

import os
import re
import argparse
import pandas as pd
from collections import Counter


# ---------------------------------------------------------------------------
# Taxid-to-name lookup (fallback if no CSV provided)
# ---------------------------------------------------------------------------

BUILTIN_TAXID_NAMES = {
    "1764": "Mycobacterium avium",
    "1773": "Mycobacterium tuberculosis",
    "1314": "Streptococcus pyogenes",
    "666":  "Vibrio cholerae",
    "5763": "Acanthamoeba",
    "13373": "Burkholderia cepacia complex",
}


# ---------------------------------------------------------------------------
# Organism name extraction from stitle
# ---------------------------------------------------------------------------

_STOP_PATTERN = re.compile(
    r"\b(strain|subsp\.|subsp\b|var\.|pv\.|chromosome|plasmid|DNA|RNA|"
    r"scaffold|contig|complete|partial|whole|genome|sequence|assembly|"
    r"node|clone|isolate|sp\.\s|cf\.\s)\b",
    re.IGNORECASE,
)


def extract_organism_from_stitle(stitle: str) -> str:
    """
    Extract the organism name from a BLAST stitle string.

    Example inputs → outputs:
      "Mycobacterium avium subsp. hominissuis JP-H-1 DNA, complete genome"
        → "Mycobacterium avium"
      "Streptococcus pyogenes strain SF370 chromosome, complete genome"
        → "Streptococcus pyogenes"
      "Burkholderia cepacia strain BC7 whole genome shotgun"
        → "Burkholderia cepacia"
    """
    if not stitle or str(stitle).strip() in ("", "N/A", "nan"):
        return "Unknown"

    stitle = str(stitle).strip()

    match = _STOP_PATTERN.search(stitle)
    if match:
        organism = stitle[: match.start()].strip().rstrip(",;")
    else:
        organism = stitle

    words = organism.split()
    return " ".join(words[:3]) if len(words) > 3 else organism


def organism_matches(hit_organism: str, expected_organism: str) -> bool:
    """
    Return True if hit_organism is consistent with expected_organism.

    Matching strategy:
      1. Direct substring match (case-insensitive)
      2. Genus-level fallback ONLY when expected organism is a single word
         (genus only). If expected has a species name, genus match alone
         is insufficient — prevents e.g. "Naegleria sp." matching
         "Naegleria fowleri".
    """
    if not hit_organism or hit_organism.lower() in ("unknown", "n/a", ""):
        return False

    hit = hit_organism.lower()
    exp = expected_organism.lower()

    if exp in hit or hit in exp:
        return True

    # Genus-level fallback only when expected is genus-only (one word)
    exp_words = exp.split()
    if len(exp_words) == 1:
        hit_genus = hit.split()[0] if hit.split() else ""
        return bool(hit_genus and hit_genus == exp_words[0])

    return False


def is_uncultured(organism: str) -> bool:
    """Return True if the hit organism is an uncultured/environmental sequence."""
    if not organism:
        return False
    lower = organism.lower()
    return any(kw in lower for kw in (
        "uncultured", "unclassified", "environmental sample",
        "metagenome", "synthetic construct"
    ))


# ---------------------------------------------------------------------------
# Per-file analysis
# ---------------------------------------------------------------------------

COLS = ["qseqid", "sscinames", "pident", "length", "evalue", "stitle"]


def analyze_tsv(filepath: str, expected_organism: str) -> dict:
    """
    Analyze a single BLAST result TSV. Considers ALL hits per read (up to 20),
    counting how many reads had any hit to each organism.

    pct_match = % of reads with ANY hit matching the expected organism.
    expected_rank = rank of expected organism by read count (1 = most reads).
    Uncultured check uses only the best (first) hit per read.
    """
    if not filepath or not os.path.exists(filepath):
        return {"status": "missing"}
    if os.path.getsize(filepath) == 0:
        return {"status": "empty"}

    try:
        df = pd.read_csv(filepath, sep="\t", header=None, names=COLS,
                         dtype=str, on_bad_lines="skip")
    except Exception as e:
        return {"status": f"read_error: {e}"}

    if df.empty:
        return {"status": "empty"}

    df["pident"] = pd.to_numeric(df["pident"], errors="coerce")
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")

    # Assign organism name to every hit row (all 20 per read)
    df["hit_organism"] = df.apply(
        lambda r: r["sscinames"]
        if pd.notna(r["sscinames"]) and str(r["sscinames"]).strip() not in ("N/A", "")
        else extract_organism_from_stitle(r["stitle"]),
        axis=1,
    )

    # Total unique reads blasted
    n_reads = df["qseqid"].nunique()

    # --- All-hits organism counts ---
    # For each read, collect the set of all organisms it hit across all 20 results.
    # Count how many reads had any hit to each organism (a read counted once per org).
    read_to_orgs = df.groupby("qseqid")["hit_organism"].apply(set)

    org_read_counts = Counter()
    for orgs in read_to_orgs:
        for org in orgs:
            org_read_counts[org] += 1

    org_counts_sorted = org_read_counts.most_common()  # list of (org, count)

    # pct_match: reads with any hit matching expected organism (across all 20 hits)
    n_match = sum(
        1 for orgs in read_to_orgs
        if any(organism_matches(o, expected_organism) for o in orgs)
    )
    pct_match = round(n_match / n_reads * 100, 1) if n_reads > 0 else 0.0

    # Rank of expected organism by read count (1 = highest)
    expected_rank = None
    for rank, (org, _) in enumerate(org_counts_sorted, 1):
        if organism_matches(org, expected_organism):
            expected_rank = rank
            break

    # --- Uncultured check (best hit per read only) ---
    top_per_read = df.groupby("qseqid", sort=False).first().reset_index()
    n_uncultured = int(top_per_read["hit_organism"].apply(is_uncultured).sum())
    pct_uncultured = round(n_uncultured / n_reads * 100, 1) if n_reads > 0 else 0.0

    # Mean % identity across all hits
    mean_pident = round(df["pident"].mean(), 2) if not df["pident"].isna().all() else None

    # Top 5 organisms by read count (all-hits)
    top_orgs = "; ".join(f"{org}({cnt})" for org, cnt in org_counts_sorted[:5])

    return {
        "status": "complete",
        "n_reads": n_reads,
        "n_match": n_match,
        "pct_match": pct_match,
        "pct_uncultured": pct_uncultured,
        "expected_rank": expected_rank,
        "mean_pident": mean_pident,
        "top_organisms": top_orgs,
    }


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def classify(r1: dict, r2: dict, threshold: float = 80.0) -> tuple:
    """
    Classify as TRUE_POSITIVE, FALSE_POSITIVE, UNCERTAIN,
    UNCULTURED_DOMINANT, or NO_DATA.

    TRUE_POSITIVE requires ALL THREE conditions in both R1 and R2:
      1. avg pct_match >= threshold (default 80%)
      2. expected organism is ranked #1 or #2 by read count
      3. n_match > 2

    UNCULTURED_DOMINANT: >50% of best hits per read are uncultured/environmental
    — takes priority over other classifications.
    """
    valid = [
        s for s in (r1, r2)
        if s and s.get("status") == "complete" and s.get("n_reads", 0) > 0
    ]
    if not valid:
        return "NO_DATA", None, None, None

    avg_pct = round(sum(s["pct_match"] for s in valid) / len(valid), 1)
    avg_pct_uncultured = round(
        sum(s.get("pct_uncultured", 0) for s in valid) / len(valid), 1
    )

    # Collect per-read-set ranks and n_match for criteria checks
    ranks = [s.get("expected_rank") for s in valid]
    n_matches = [s.get("n_match", 0) for s in valid]

    # Uncultured check takes priority
    if avg_pct_uncultured > 50.0:
        return "UNCULTURED_DOMINANT", avg_pct, avg_pct_uncultured, ranks

    if avg_pct >= threshold:
        # Criterion 2: expected organism must be rank 1 or 2 in every read set
        rank_ok = all(r is not None and r <= 2 for r in ranks)
        # Criterion 3: n_match > 2 in every read set
        reads_ok = all(n > 2 for n in n_matches)

        if rank_ok and reads_ok:
            label = "TRUE_POSITIVE"
        else:
            label = "UNCERTAIN"
    elif avg_pct < 10.0:
        label = "FALSE_POSITIVE"
    else:
        label = "UNCERTAIN"

    return label, avg_pct, avg_pct_uncultured, ranks


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Classify BLAST results as true/false positives per sample/taxid"
    )
    parser.add_argument("--blast_dir", required=True,
                        help="Root directory containing blast_validation_SAMPLE_TAXID/ folders")
    parser.add_argument("--taxid_list", required=True,
                        help="TSV with columns: sample, taxid (no header)")
    parser.add_argument("--taxid_names", default=None,
                        help="CSV with columns: taxid, organism_name. "
                             "If omitted, built-in lookup is used.")
    parser.add_argument("--output", required=True,
                        help="Output TSV report path")
    parser.add_argument("--threshold", type=float, default=80.0,
                        help="Min %% of reads matching expected organism "
                             "to call TRUE_POSITIVE (default: 80)")
    args = parser.parse_args()

    # --- Load taxid-to-name mapping ---
    if args.taxid_names and os.path.exists(args.taxid_names):
        tn = pd.read_csv(args.taxid_names)
        taxid_to_name = dict(zip(tn["taxid"].astype(str), tn["organism_name"]))
        print(f"Loaded {len(taxid_to_name)} taxid names from {args.taxid_names}")
    else:
        taxid_to_name = BUILTIN_TAXID_NAMES
        print(f"No taxid_names CSV provided — using built-in lookup ({len(taxid_to_name)} entries)")

    # --- Load sample/taxid list ---
    samples_df = pd.read_csv(args.taxid_list, sep="\t", header=None,
                              names=["sample", "taxid"])
    samples_df["taxid"] = samples_df["taxid"].astype(str)
    total = len(samples_df)
    print(f"Processing {total} sample/taxid combinations...\n")

    rows = []
    for i, (_, row) in enumerate(samples_df.iterrows(), 1):
        sample = row["sample"]
        taxid = row["taxid"]
        expected = taxid_to_name.get(taxid, f"taxid_{taxid}")

        outdir = os.path.join(args.blast_dir, f"blast_validation_{sample}_{taxid}")
        r1_path = os.path.join(outdir, f"{sample}_taxid_{taxid}_R1_blast_results.tsv")
        r2_path = os.path.join(outdir, f"{sample}_taxid_{taxid}_R2_blast_results.tsv")

        r1 = analyze_tsv(r1_path, expected)
        r2 = analyze_tsv(r2_path, expected)

        label, avg_pct, avg_pct_uncultured, ranks = classify(r1, r2, args.threshold)

        # Override if BLAST was never run (no directory at all)
        if not os.path.isdir(outdir):
            label = "BLAST_NOT_RUN"

        r1_rank = ranks[0] if ranks and len(ranks) > 0 else None
        r2_rank = ranks[1] if ranks and len(ranks) > 1 else None

        if i % 500 == 0 or i == total:
            print(f"  {i}/{total} done...")

        rows.append({
            "sample":              sample,
            "taxid":               taxid,
            "expected_organism":   expected,
            "classification":      label,
            "avg_pct_match":       avg_pct,
            "avg_pct_uncultured":  avg_pct_uncultured,
            "R1_n_reads":          r1.get("n_reads"),
            "R1_n_match":          r1.get("n_match"),
            "R1_pct_match":        r1.get("pct_match"),
            "R1_pct_uncultured":   r1.get("pct_uncultured"),
            "R1_expected_rank":    r1.get("expected_rank"),
            "R1_mean_pident":      r1.get("mean_pident"),
            "R1_top_organisms":    r1.get("top_organisms"),
            "R2_n_reads":          r2.get("n_reads"),
            "R2_n_match":          r2.get("n_match"),
            "R2_pct_match":        r2.get("pct_match"),
            "R2_pct_uncultured":   r2.get("pct_uncultured"),
            "R2_expected_rank":    r2.get("expected_rank"),
            "R2_mean_pident":      r2.get("mean_pident"),
            "R2_top_organisms":    r2.get("top_organisms"),
        })

    out_df = pd.DataFrame(rows)

    # --- Summary ---
    print("\n" + "=" * 50)
    print("CLASSIFICATION SUMMARY")
    print("=" * 50)
    counts = out_df["classification"].value_counts()
    for label, count in counts.items():
        print(f"  {label:<25} {count:>5}  ({count/total*100:.1f}%)")
    print(f"  {'TOTAL':<25} {total:>5}")

    fp = out_df[out_df["classification"] == "FALSE_POSITIVE"]
    if not fp.empty:
        print(f"\nFALSE POSITIVES ({len(fp)}):")
        for _, r in fp.iterrows():
            print(f"  {r['sample']}  taxid {r['taxid']}  "
                  f"({r['expected_organism']})  "
                  f"avg_match={r['avg_pct_match']}%")
            if r["R1_top_organisms"]:
                print(f"    R1 top hits: {r['R1_top_organisms']}")
            if r["R2_top_organisms"]:
                print(f"    R2 top hits: {r['R2_top_organisms']}")

    unc_dom = out_df[out_df["classification"] == "UNCULTURED_DOMINANT"]
    if not unc_dom.empty:
        print(f"\nUNCULTURED_DOMINANT ({len(unc_dom)}) — top hits are uncultured/environmental:")
        for _, r in unc_dom.iterrows():
            print(f"  {r['sample']}  taxid {r['taxid']}  "
                  f"({r['expected_organism']})  "
                  f"avg_match={r['avg_pct_match']}%  "
                  f"avg_uncultured={r['avg_pct_uncultured']}%")

    uncertain = out_df[out_df["classification"] == "UNCERTAIN"]
    if not uncertain.empty:
        print(f"\nUNCERTAIN ({len(uncertain)}) — review manually:")
        for _, r in uncertain.iterrows():
            print(f"  {r['sample']}  taxid {r['taxid']}  "
                  f"({r['expected_organism']})  "
                  f"avg_match={r['avg_pct_match']}%  "
                  f"R1_rank={r['R1_expected_rank']}  R2_rank={r['R2_expected_rank']}")

    out_df.to_csv(args.output, sep="\t", index=False)
    print(f"\nFull report written to: {args.output}")


if __name__ == "__main__":
    main()
