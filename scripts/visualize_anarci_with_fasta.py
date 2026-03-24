#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
visualize_anarci_with_fasta.py

CLI tool to print an aligned, monospaced visualization of antibody numbering schemes
(Kabat/Chothia/Martin/IMGT) overlaid on the ORIGINAL full sequence (FASTA).
- Uses the FASTA as the backbone so linkers/unnumbered residues are shown.
- Accepts a single ANARCI text output that may contain multiple "# Scheme = ..." blocks
  (concatenate multiple ANARCI runs with `cat` to combine schemes).
- Prints labeled rows and supports coordinated wrapping.

Example:
  # Make a multi-scheme ANARCI file
  ANARCI -i my.fasta --scheme kabat  -o kabat.txt
  ANARCI -i my.fasta --scheme martin -o martin.txt
  cat kabat.txt martin.txt > my_multi_schemes.txt

  # Visualize (wrap at 50 residues per block, save to file too)
  python visualize_anarci_with_fasta.py --anarci my_multi_schemes.txt --fasta my.fasta \
      --wrap 50 --outfile out.txt
"""

import re
import argparse
from typing import Dict, List, Tuple, Optional
from contextlib import contextmanager

# ---- CDR boundaries (numeric part; insertions like 52A count as 52) ----
CDR_RANGES = {
    "kabat":   {"H": [(31, 35), (50, 65), (95, 102)], "L": [(24, 34), (50, 56), (89, 97)]},
    "chothia": {"H": [(26, 32), (52, 56), (95, 102)], "L": [(24, 34), (50, 56), (89, 97)]},
    "martin":  {"H": [(26, 32), (52, 56), (95, 102)], "L": [(24, 34), (50, 56), (89, 97)]},
#wrong i think    "imgt":    {"H": [(27, 38), (56, 65), (105, 117)], "L": [(27, 38), (56, 65), (105, 117)]},
}

SCHEME_PRINT_ORDER = ["kabat", "chothia", "martin", "imgt"]
SCHEME_MASK_LETTER = {"kabat": "K", "chothia": "C", "martin": "M", "imgt": "I"}


def _fmtw(x, w=7):
    return f"{str(x):>{w}}"


def _in_cdr(scheme: str, chain: str, num_str: str) -> bool:
    scheme = scheme.lower()
    if (scheme not in CDR_RANGES) or (chain not in CDR_RANGES[scheme]) or (not num_str):
        return False
    m = re.match(r"(\d+)", str(num_str))
    if not m:
        return False
    base_num = int(m.group(1))
    for start, end in CDR_RANGES[scheme][chain]:
        if start <= base_num <= end:
            return True
    return False


def _read_fasta_single(path: str) -> str:
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq).replace(" ", "").replace("\t", "")


def visualize_anarci_with_fasta(
    anarci_txt_path: str,
    fasta_path: str,
    width: int = 7,             # per-residue column width (monospaced)
    row_wrap: Optional[int] = None,  # residues per block (coordinated wrapping)
    outfile: Optional[str] = None,
    label_width: int = 16,      # left label column width
) -> str:
    """
    Use the ORIGINAL full sequence (FASTA) as the backbone (positions 1..N).
    Overlay ANARCI numbering for any present schemes (Kabat/Chothia/Martin/IMGT).
    Intervening/linker residues remain with blank numbering/masks.

    Prints labeled rows:
      1) "Residue #" (sequential index)
      2..(1+S) "<Scheme> #" per present scheme, in SCHEME_PRINT_ORDER
      next) "Residue" (AA sequence from FASTA)
      last S) "<Scheme> CDR?" per scheme (K/C/M/I letters mark CDRs)
    """
    # --- Parse ANARCI text into data[scheme][chain] = list of (num_str, aa) ---
    data: Dict[str, Dict[str, List[Tuple[str, str]]]] = {}
    scheme = None
    # Lines look like: "H 52      N" or "H 52    A L"
    line_re = re.compile(r'^\s*([A-Za-z])\s+(\d+)(?:\s+([A-Z]))?\s+([A-Z])\s*$')
    scheme_re = re.compile(r"Scheme\s*=\s*([A-Za-z]+)", re.IGNORECASE)

    with open(anarci_txt_path, "r") as f:
        for raw in f:
            line = raw.rstrip("\n")
            msch = scheme_re.search(line)
            if msch:
                scheme = msch.group(1).strip().lower()
                data.setdefault(scheme, {"H": [], "L": []})  # non-IG rows ignored here
                continue
            mres = line_re.match(line)
            if mres and scheme is not None:
                chain_label, num, ins, aa = mres.groups()
                ch = chain_label.upper()
                if ch not in ("H", "L"):
                    # Ignore non-IG chains in ANARCI output. Linkers come from FASTA.
                    continue
                num_str = f"{num}{ins or ''}"
                data[scheme][ch].append((num_str, aa))

    present_schemes = [s for s in SCHEME_PRINT_ORDER if s in data and (data[s]["H"] or data[s]["L"])]

    # --- Read full original sequence (backbone) ---
    full_seq = _read_fasta_single(fasta_path)
    N = len(full_seq)
    if N == 0:
        out = "Original FASTA sequence is empty.\n"
        print(out)
        if outfile:
            with open(outfile, "w") as h:
                h.write(out)
        return out

    # --- Master H/L sequence from first scheme that has them (order definition) ---
    def first_chain_seq(ch: str) -> str:
        for s in present_schemes:
            if data[s][ch]:
                return "".join(aa for (_, aa) in data[s][ch])
        return ""

    H_seq, L_seq = first_chain_seq("H"), first_chain_seq("L")

    # --- Place H and L inside full sequence (VH…VL or VL…VH) ---
    def _find_nonoverlap(a: str, b: str, text: str):
        if not a or not b:
            return None
        pA = text.find(a)
        if pA >= 0:
            pB = text.find(b, pA + len(a))
            if pB >= 0:
                return (pA, pB)
        return None

    posH, posL = None, None
    if H_seq and L_seq:
        order1 = _find_nonoverlap(H_seq, L_seq, full_seq)  # VH…VL
        order2 = _find_nonoverlap(L_seq, H_seq, full_seq)  # VL…VH
        if order1:
            posH, posL = order1
        elif order2:
            posL, posH = order2
        else:
            posH = full_seq.find(H_seq) if H_seq else -1
            posL = full_seq.find(L_seq) if L_seq else -1
    elif H_seq:
        posH = full_seq.find(H_seq)
    elif L_seq:
        posL = full_seq.find(L_seq)

    # Build base rows
    idx_row = "".join(_fmtw(i + 1, width) for i in range(N))
    seq_row = "".join(_fmtw(a, width) for a in full_seq)

    numbering_rows: Dict[str, List[str]] = {s: [""] * N for s in present_schemes}
    mask_rows: Dict[str, List[str]]      = {s: [" "] * N for s in present_schemes}

    def overlay_chain(chain: str, start_pos: int):
        if start_pos is None or start_pos < 0:
            return
        master_len = len(H_seq if chain == "H" else L_seq)
        if master_len == 0:
            return
        for s in present_schemes:
            nums = [num for (num, _) in data[s][chain]]
            aas  = [aa  for (num, _) in data[s][chain]]
            Ls = min(master_len, len(nums), len(aas))
            for i in range(Ls):
                pos = start_pos + i
                if pos >= N:
                    break
                numbering_rows[s][pos] = nums[i]
                if _in_cdr(s, chain, nums[i]):
                    mask_rows[s][pos] = SCHEME_MASK_LETTER[s]

    if H_seq and posH is not None and posH >= 0:
        overlay_chain("H", posH)
    if L_seq and posL is not None and posL >= 0:
        overlay_chain("L", posL)

    numbering_rows_str = {s: "".join(_fmtw(x, width) for x in numbering_rows[s]) for s in present_schemes}
    mask_rows_str      = {s: "".join(_fmtw(x, width) for x in mask_rows[s])      for s in present_schemes}

    # ---- Build labeled rows (labels appear at the start of every wrapped block) ----
    rows_labels = []
    rows_data   = []

    rows_labels.append("Residue #")
    rows_data.append(idx_row)

    for s in present_schemes:
        rows_labels.append(f"{s.title()} #")
        rows_data.append(numbering_rows_str[s])

    rows_labels.append("Residue")
    rows_data.append(seq_row)

    for s in present_schemes:
        rows_labels.append(f"{s.title()} CDR?")
        rows_data.append(mask_rows_str[s])

    # Output
    @contextmanager
    def maybe_open(path):
        if path:
            with open(path, "w") as h:
                yield h
        else:
            yield None

    with maybe_open(outfile) as fout:
        def _outln(text=""):
            print(text)
            if fout:
                fout.write(text + "\n")

        header = f"=== Original sequence backbone (schemes: {', '.join(present_schemes)}) ==="
        _outln(header)

        # Coordinated wrapping across all rows (two blank lines between blocks)
        if row_wrap and row_wrap > 0:
            for start in range(0, N, row_wrap):
                end = min(start + row_wrap, N)
                a, b = start * width, end * width
                for label, row in zip(rows_labels, rows_data):
                    _outln(f"{label:<{label_width}}{row[a:b]}")
                _outln("")
                _outln("")
        else:
            for label, row in zip(rows_labels, rows_data):
                _outln(f"{label:<{label_width}}{row}")

    return "OK"


def main():
    ap = argparse.ArgumentParser(
        description="Visualize ANARCI numbering (Kabat/Chothia/Martin/IMGT) over a full FASTA sequence."
    )
    ap.add_argument("--anarci", required=True, help="ANARCI text output (may contain multiple '# Scheme =' blocks)")
    ap.add_argument("--fasta",  required=True, help="Original full sequence (single-entry FASTA)")
    ap.add_argument("--wrap",   type=int, default=0, help="Residues per wrapped block (0 = no wrap)")
    ap.add_argument("--width",  type=int, default=7, help="Per-residue column width (monospaced)")
    ap.add_argument("--label-width", type=int, default=16, help="Left label column width")
    ap.add_argument("--outfile", default=None, help="Optional path to also write the output text")
    args = ap.parse_args()

    visualize_anarci_with_fasta(
        anarci_txt_path=args.anarci,
        fasta_path=args.fasta,
        width=args.width,
        row_wrap=(args.wrap if args.wrap and args.wrap > 0 else None),
        outfile=args.outfile,
        label_width=args.label_width,
    )


if __name__ == "__main__":
    main()
