#!/usr/bin/env python3
"""
extract_chain_seq.py

Extracts the 1-letter amino acid sequence for a specified chain
from a PDB or mmCIF file.

Usage:
    python extract_chain_seq.py <input_file> <chain_id>

Outputs:
    The amino acid sequence on stdout (no header).
"""

import sys
import gemmi
from Bio.PDB import PDBParser
# --- Robust import for three_to_one across Biopython versions ---
try:
    from Bio.PDB.Polypeptide import three_to_one
except ImportError:
    try:
        from Bio.SeqUtils import seq1 as three_to_one
    except ImportError:
        from Bio.PDB.Polypeptide import protein_letters_3to1

        def three_to_one(resname):
            resname = (resname or "").strip().upper()
            if resname == "MSE":
                resname = "MET"
            return protein_letters_3to1.get(resname.title(), "X")

if len(sys.argv) != 3:
    print("Usage: python extract_chain_seq.py <input_file> <chain_id>")
    sys.exit(1)

path = sys.argv[1]
chain_id = sys.argv[2]
ext = path.lower().split(".")[-1]


def seq_from_pdb(pdb_path, chain):
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("x", pdb_path)
    except Exception as e:
        sys.stderr.write(f"Error reading PDB: {e}\n")
        return ""

    for model in structure:
        for c in model:
            if c.id == chain:
                seq = ""
                for res in c:
                    try:
                        seq += three_to_one(res.get_resname())
                    except KeyError:
                        pass
                return seq
    return ""


def seq_from_cif(cif_path, chain):
    try:
        doc = gemmi.cif.read_file(cif_path)
        structure = gemmi.make_structure_from_block(doc.sole_block())
    except Exception as e:
        sys.stderr.write(f"Error reading CIF: {e}\n")
        return ""

    for model in structure:
        for c in model:
            if c.name == chain:
                seq = ""
                for res in c:
                    if res.is_amino_acid():
                        seq += gemmi.residue_name_3to1.get(res.name, "")
                return seq
    return ""


# Select correct parser
if ext == "pdb":
    seq = seq_from_pdb(path, chain_id)
elif ext in ("cif", "mmcif", "mcif"):
    seq = seq_from_cif(path, chain_id)
else:
    sys.stderr.write(f"Unsupported file extension: {ext}\n")
    seq = ""

# Output sequence
print(seq)
