#!/usr/bin/env python3
import re
import argparse
from textwrap import wrap

# --- Martin CDR ranges (including insertions) -------------------------------
# Heavy: H1 = 26–32 ; H2 = 52–56 ; H3 = 95–102
# Light: L1 = 24–34 ; L2 = 50–56 ; L3 = 89–97
CDR_RANGES_MARTIN = {
    "H": [("H1", 26, 32), ("H2", 52, 56), ("H3", 95, 102)],
    "L": [("L1", 24, 34), ("L2", 50, 56), ("L3", 89, 97)],
}

ROW_WIDTH = 80  # characters per printed row


def parse_anarci_txt(path):
    """Parse ANARCI text output (like you pasted) into per-chain lists of residues."""
    chains = {"H": [], "L": []}
    seq_idx = {"H": 0, "L": 0}

    line_re = re.compile(r'^\s*([HL])\s+(\d+)(?:\s+([A-Z]))?\s+([A-Z])\s*$')

    with open(path) as f:
        for line in f:
            m = line_re.match(line)
            if not m:
                continue
            chain, num, ins, aa = m.groups()
            seq_idx[chain] += 1
            num_str = f"{num}{ins or ''}"
            chains[chain].append((seq_idx[chain], num_str, aa))
    return chains


def in_cdr(chain, num_str):
    """Return True if a Martin-numbered position is in a CDR for that chain."""
    base_num = int(re.match(r'(\d+)', num_str).group(1))
    for _, start, end in CDR_RANGES_MARTIN[chain]:
        if start <= base_num <= end:
            return True
    return False


def build_rows(chain_residues, chain_letter):
    """Build formatted rows for sequential index, Martin numbering, and CDR mask."""
    def fmt(x):
        return f"{x:>4}"

    row1 = "".join(fmt(i) for (i, num, aa) in chain_residues)
    row2 = "".join(fmt(num) for (i, num, aa) in chain_residues)
    row3 = "".join(fmt("X" if in_cdr(chain_letter, num) else "") for (i, num, aa) in chain_residues)
    rowAA = "".join(fmt(aa) for (i, num, aa) in chain_residues)

    return row1, row2, row3, rowAA


def print_blocked(title, rows, output_file=None):
    header = f"\n=== {title} ===\n"
    print(header)
    if output_file:
        output_file.write(header)

    for r in rows:
        for chunk in wrap(r, ROW_WIDTH):
            print(chunk)
            if output_file:
                output_file.write(chunk + "\n")
        print()
        if output_file:
            output_file.write("\n")


def main():
    parser = argparse.ArgumentParser(description="Show seq index / Martin numbering / CDR mask from ANARCI text output.")
    parser.add_argument("anarci_txt", help="Path to ANARCI text output (like the one you pasted).")
    parser.add_argument("-o", "--out", default="anarci_cdrs.txt", help="Path to save formatted output.")
    parser.add_argument("--show-aa", action="store_true", help="Also include amino-acid row.")
    args = parser.parse_args()

    chains = parse_anarci_txt(args.anarci_txt)

    with open(args.out, "w") as fout:
        for chain_letter in ("H", "L"):
            if chains[chain_letter]:
                row1, row2, row3, rowAA = build_rows(chains[chain_letter], chain_letter)
                rows = [f"SeqIdx:{row1}", f"Martin:{row2}", f"CDRs  :{row3}"]
                if args.show_aa:
                    rows.append(f"AA    :{rowAA}")
                print_blocked(f"Chain {chain_letter}", rows, fout)
            else:
                msg = f"\n=== Chain {chain_letter} ===\n(no residues parsed)\n"
                print(msg)
                fout.write(msg)

    print(f"\nOutput saved to {args.out}")


if __name__ == "__main__":
    main()