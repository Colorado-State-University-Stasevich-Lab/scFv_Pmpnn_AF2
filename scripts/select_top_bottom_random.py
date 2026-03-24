#!/usr/bin/env python3
"""
Select top-N, bottom-N, and random-Nr sequences based on a chosen score column.

Outputs:
    selected_top.csv
    selected_bottom.csv
    selected_random.csv
    selected_for_af2.txt   (list of Accessions)
"""

import pandas as pd
import argparse
import numpy as np

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="merged_preAF2_scores.csv")
    p.add_argument("--sortcol", required=True, help="Column to sort by (default=scfvtools_score)")
    p.add_argument("--Nt", type=int, required=True)
    p.add_argument("--Nb", type=int, required=True)
    p.add_argument("--Nr", type=int, required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    df = pd.read_csv(args.input)

    if args.sortcol not in df.columns:
        raise ValueError(f"sortcol '{args.sortcol}' not found in columns: {df.columns}")

    # Ensure no NA in sort column (rank misses last)
    df = df.copy()
    df[args.sortcol] = pd.to_numeric(df[args.sortcol], errors="coerce")

    # Sort
    df_sorted = df.sort_values(by=args.sortcol, ascending=False)  # Higher = better

    # Top Nt
    top = df_sorted.head(args.Nt)
    top.to_csv(f"{args.outdir}/selected_top.csv", index=False)

    # Bottom Nb
    bottom = df_sorted.tail(args.Nb)
    bottom.to_csv(f"{args.outdir}/selected_bottom.csv", index=False)

    # Random Nr
    random = df.sample(n=args.Nr, replace=False, random_state=37)
    random.to_csv(f"{args.outdir}/selected_random.csv", index=False)

    # Combined list for AF2
    combined = pd.concat([top, bottom, random], ignore_index=True)
    combined["Accession"].drop_duplicates().to_csv(
        f"{args.outdir}/selected_for_af2.txt",
        index=False,
        header=False
    )

    print("[select_top_bottom_random] Selection complete.")
    print(f"Top → {args.outdir}/selected_top.csv")
    print(f"Bottom → {args.outdir}/selected_bottom.csv")
    print(f"Random → {args.outdir}/selected_random.csv")
    print(f"Selected list → {args.outdir}/selected_for_af2.txt")

if __name__ == "__main__":
    main()
