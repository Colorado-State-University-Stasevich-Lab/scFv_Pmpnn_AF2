import pandas as pd
import argparse

def reorder_scores(input_file, output_file):
    """Reorders all_scores.csv based on pLDDT, log_likelihood, and Prob. of Solubility."""
    # Load the CSV
    df = pd.read_csv(input_file)

    # Convert log_likelihood to numerical (if stored as string)
    df["log_likelihood"] = pd.to_numeric(df["log_likelihood"], errors="coerce")
    df["pLDDT"] = pd.to_numeric(df["pLDDT"], errors="coerce")
    df["Prob. of Solubility"] = pd.to_numeric(df["Prob. of Solubility"], errors="coerce")

    # Sort by pLDDT (highest first), log_likelihood (highest first), and Prob. of Solubility (highest first)
    df = df.sort_values(by=["pLDDT", "log_likelihood", "Prob. of Solubility"], ascending=[False, False, False])

    # Save the reordered file
    df.to_csv(output_file, index=False)
    print(f"Reordered file saved as {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reorder all_scores.csv by pLDDT, log_likelihood, and Prob. of Solubility.")
    parser.add_argument("--input", type=str, required=True, help="Path to input CSV file (all_scores.csv).")
    parser.add_argument("--output", type=str, required=True, help="Path to output reordered CSV file.")
    
    args = parser.parse_args()
    reorder_scores(args.input, args.output)
