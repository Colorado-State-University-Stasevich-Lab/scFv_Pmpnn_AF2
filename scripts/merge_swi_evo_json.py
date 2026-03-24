import os
import json
import pandas as pd
import argparse
import numpy as np

def load_swi_file(swi_file):
    """Loads the SWI CSV file and removes the Sequence column."""
    swi_df = pd.read_csv(swi_file)
    if "Sequence" in swi_df.columns:
        swi_df = swi_df.drop(columns=["Sequence"])  # Remove the sequence column
    return swi_df

def load_evo_file(evo_file):
    """Loads the evolutionary scores CSV file."""
    return pd.read_csv(evo_file)

def _mean_safe(x):
    try:
        arr = np.array(x, dtype=float)
        if arr.size == 0:
            return None
        return float(np.nanmean(arr))
    except Exception:
        return None

from pathlib import Path

def extract_json_data(json_dir, seqids):
    """
    Extract mean pLDDT, mean pAE, pTM, and ipTM from rank_001 JSON files in json_dir.
    Uses EXACT seqid matching by parsing filenames as:
        <seqid>_scores_rank_001_*.json
    This avoids collisions like design_42 vs design_427.
    """
    json_dir = Path(json_dir)

    # Build an exact seqid -> json_path map
    json_map = {}
    duplicates = []

    for p in json_dir.glob("*_scores_rank_001_*.json"):
        fname = p.name
        seqid_from_file = fname.split("_scores_rank_")[0]

        if seqid_from_file in json_map:
            duplicates.append(seqid_from_file)
            # Keep the first one deterministically (or you could overwrite / warn)
            continue

        json_map[seqid_from_file] = p

    if duplicates:
        # Not fatal, but good to know
        print(f"[WARN] Multiple rank_001 JSONs found for {len(set(duplicates))} seqids (kept first).")

    json_data = []
    for seqid in seqids:
        p = json_map.get(str(seqid))

        if p is None:
            json_data.append({"seqid": seqid, "pLDDT": None, "pAE_mean": None, "pTM": None, "ipTM": None})
            continue

        plddt_mean = pae_mean = ptm = iptm = None

        try:
            with open(p, "r") as f:
                data = json.load(f)

            # pLDDT: list of per-residue scores
            plddt_vals = data.get("plddt", [])
            plddt_mean = _mean_safe(plddt_vals)

            # pAE: often a 2D list (matrix). We report mean of all entries.
            pae_vals = data.get("pae", [])
            if isinstance(pae_vals, list) and len(pae_vals) > 0 and isinstance(pae_vals[0], list):
                pae_flat = [v for row in pae_vals for v in row]
            else:
                pae_flat = pae_vals
            pae_mean = _mean_safe(pae_flat)

            # pTM / ipTM
            ptm = data.get("ptm", None)
            iptm = data.get("iptm", None)
            ptm = float(ptm) if ptm is not None else None
            iptm = float(iptm) if iptm is not None else None

        except Exception as e:
            print(f"Error processing {p}: {e}")

        json_data.append({
            "seqid": seqid,
            "pLDDT": plddt_mean,
            "pAE_mean": pae_mean,
            "pTM": ptm,
            "ipTM": iptm
        })

    return pd.DataFrame(json_data)

def _load_accession_set(path):
    """
    Load a CSV and return a set of accessions.
    Accepts either an 'Accession' column or uses the first column.
    """
    if path is None:
        return set()
    df = pd.read_csv(path)
    if df.empty:
        return set()
    col = "Accession" if "Accession" in df.columns else df.columns[0]
    return set(df[col].astype(str))


def main():
    parser = argparse.ArgumentParser(
        description="Merge SWI, evolutionary scores, and AF2/ColabFold JSON metrics (pLDDT, pAE mean, pTM, ipTM)."
    )
    parser.add_argument("--swi_file", type=str, required=True, help="Path to the SWI CSV file.")
    parser.add_argument("--evo_file", type=str, required=True, help="Path to the evolutionary scores CSV file.")
    parser.add_argument("--json_dir", type=str, required=True, help="Directory containing LocalColabFold JSON files.")
    parser.add_argument("--outpath", type=str, required=True, help="Path for the merged output CSV file.")
    parser.add_argument("--scfvtools_file", type=str, required=False,
                    help="Optional: CSV file containing scfvtools per-sequence scores.")
    parser.add_argument("--top", type=str, default=None,
                    help="Optional: CSV listing selected top designs (e.g., selected_top.csv).")
    parser.add_argument("--bottom", type=str, default=None,
                    help="Optional: CSV listing selected bottom designs (e.g., selected_bottom.csv).")
    parser.add_argument("--random", type=str, default=None,
                    help="Optional: CSV listing selected random designs (e.g., selected_random.csv).")

    args = parser.parse_args()

    # Load inputs
    swi_df = load_swi_file(args.swi_file)
    evo_df = load_evo_file(args.evo_file)

    # Merge on Accession (SWI) and seqid (evo)
    merged_df = swi_df.merge(evo_df, left_on="Accession", right_on="seqid", how="inner")

    # ---------------------------------------------------------------
    # Optional: merge SCFVTOOLS scores
    # ---------------------------------------------------------------
    if args.scfvtools_file:
        try:
            scfv_df = pd.read_csv(args.scfvtools_file)

            # rename to match merge key
            if "name" in scfv_df.columns:
                scfv_df = scfv_df.rename(columns={"name": "Accession"})

            # rename Score → scfvtools_score
            if "Score" in scfv_df.columns:
                scfv_df = scfv_df.rename(columns={"Score": "scfvtools_score"})

            merged_df = merged_df.merge(
                scfv_df[["Accession", "scfvtools_score", "scfvtools_blosum_score",
    "scfvtools_blosum_diff_score"]],
                on="Accession",
                how="left"
            )

            print("Merged SCFVTOOLS scores.")

        except Exception as e:
            print(f"Warning: Failed to merge scfvtools scores: {e}")



    # Pull JSON metrics
    json_df = extract_json_data(args.json_dir, merged_df["seqid"])

    # Final merge
    final_df = merged_df.merge(json_df, on="seqid", how="left")

    # ---------------------------------------------------------------
    # Optional: add selection labels (TOP/BOTTOM/RANDOM)
    # ---------------------------------------------------------------
    final_df["selection"] = ""

    top_set = _load_accession_set(args.top)
    bottom_set = _load_accession_set(args.bottom)
    random_set = _load_accession_set(args.random)

    # Precedence: TOP > BOTTOM > RANDOM
    if top_set:
        final_df.loc[final_df["Accession"].isin(top_set), "selection"] = "TOP"
    if bottom_set:
        final_df.loc[final_df["Accession"].isin(bottom_set), "selection"] = "BOTTOM"

    # Only fill RANDOM where still blank
    if random_set:
        mask_blank = final_df["selection"].eq("")
        final_df.loc[mask_blank & final_df["Accession"].isin(random_set), "selection"] = "RANDOM"


    # Sort by scfvtools_score
    if "scfvtools_blosum_diff_score" in final_df.columns:
        final_df = final_df.sort_values("scfvtools_blosum_diff_score", ascending=False)
    

    # Save
    final_df.to_csv(args.outpath, index=False)
    print(f"Merged results saved to {args.outpath}")

if __name__ == "__main__":
    main()
