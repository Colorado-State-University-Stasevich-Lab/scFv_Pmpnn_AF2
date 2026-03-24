#!/usr/bin/env python3
import sys
from pathlib import Path

if len(sys.argv) != 4:
    print("Usage: make_combined_multimer_fasta.py <seqs_dir> <epitope_seq> <output_fasta>")
    sys.exit(1)

seqs_dir = Path(sys.argv[1])
epitope_raw = sys.argv[2]
outfile = Path(sys.argv[3])

# Normalize epitope string
epitope = epitope_raw.strip()

if epitope.lower() == "none" or epitope == "":
    epitope = None
else:
    # allow comma-separated input as well
    epitope = epitope.replace(" ", "")
    epitope = epitope.replace(",", ":")

# Find the only .fa file in seqs_dir
fas = list(seqs_dir.glob("*.fa"))
if len(fas) == 0:
    raise FileNotFoundError(f"No FASTA files found in {seqs_dir}")
if len(fas) > 1:
    raise RuntimeError(f"Expected exactly one FASTA in {seqs_dir}, found {len(fas)}")

fa = fas[0]
base = fa.stem  # basename for labeling

def format_seq(s):
    if epitope is None:
        return s
    else:
        return f"{s}:{epitope}"

with open(fa) as f, open(outfile, "w") as out:
    lines = [l.strip() for l in f.readlines() if l.strip()]

    i = 0
    record_index = 0

    while i < len(lines):
        header = lines[i]
        seq = lines[i+1] if i+1 < len(lines) else ""
        i += 2

        if not header.startswith(">"):
            continue

        if record_index == 0:
            new_header = f">{base}_WT"
        else:
            new_header = f">{base}_design_{record_index}"

        out.write(new_header + "\n")
        out.write(format_seq(seq) + "\n")

        record_index += 1

print(f"Wrote combined FASTA to {outfile}")
