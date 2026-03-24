#!/usr/bin/env python3
"""
Merge SWI, struct-evo, and scfvtools scores BEFORE running AF2.
No AF2 JSON fields included.

Outputs: merged_preAF2_scores.csv
"""

import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--swi", required=True, help="SWI CSV file")
    parser.add_argument("--evo", required=True, help="struct-evo CSV file")
    parser.add_argument("--scfv", required=True, help="scfvtools CSV file")
    parser.add_argument("--out", required=True, help="Merged output CSV")
    parser.add_argument(
        "--sort_by",
        default="scfvtools_blosum_diff_score",
        help="Column to sort by (default: scfvtools_blosum_diff_score)"
    )
    args = parser.parse_args()

    # Load
    swi = pd.read_csv(args.swi)          # SWI should have Accession
    evo = pd.read_csv(args.evo)          # struct-evo should have seqid
    scfv = pd.read_csv(args.scfv)        # scfvtools_score + Accession

    # Merge SWI + struct-evo by Accession / seqid
    merged = swi.merge(evo, left_on="Accession", right_on="seqid", how="inner")


    # Merge in scfvtools scores (hard + soft, if present)
    wanted = [
        "Accession",
        "scfvtools_score",
        "scfvtools_blosum_score",
        "scfvtools_blosum_diff_score",
    ]
    have = [c for c in wanted if c in scfv.columns]

    if "Accession" in have and "scfvtools_score" in have:
        merged = merged.merge(
            scfv[have],
            on="Accession",
            how="left",
        )
    else:
        print(f"[WARNING] scfvtools file missing expected columns. Have: {list(scfv.columns)}")

    # Explicit sort
    if args.sort_by in merged.columns:
        merged = merged.sort_values(
            by=args.sort_by,
            ascending=False,
            na_position="last"
        )
    else:
        print(f"[WARNING] sort column not found: {args.sort_by}")

    merged.to_csv(args.out, index=False)
    print(f"[merge_preAF2] Wrote merged scores → {args.out}")

if __name__ == "__main__":
    main()
