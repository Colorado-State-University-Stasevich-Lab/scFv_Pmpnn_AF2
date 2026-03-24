#!/usr/bin/env python3
"""
make_scfv_vernier.py  —  simplified “linear map” version

Key idea:
- When reading a PDB/CIF, build a linear FASTA string *and* a parallel list of Residue objects.
- That gives a perfect 1:1 map: fasta_idx i  <->  pdb_reslist[i].
- Later, map VH/VL back to PDB by exact substring (optionally after slicing by the linker).

Features kept:
- Run ANARCI, extract first VH & VL (with numbering)
- Colorize CDRs in terminal
- Optional scFv build: MA + VH + linker + VL + RA
- Fixed positions (CDRs, LINKER, FRAMEWORK)
- Vernier shell from structure using Bio.PDB NeighborSearch
- TSV dumps: pdb_linear_map.tsv, VH_mapping.tsv, VL_mapping.tsv (optional)

Examples:
  python make_scfv_vernier.py input.pdb   --scheme martin --combine --vernier-cutoff 3.0 --dump-mapping
  python make_scfv_vernier.py input.fasta --scheme kabat  --combine
"""

import argparse, subprocess, tempfile, os, sys, itertools
from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa
try:
    from Bio.PDB.PPBuilder import PPBuilder      # some installs
except ImportError:
    from Bio.PDB.Polypeptide import PPBuilder    # others (incl. yours)


# ---- ANSI colors for terminal preview ----
C_GREEN  = "\033[92m"   # L CDRs
C_MAGENTA= "\033[95m"   # H CDRs
C_GRAY   = "\033[90m"   # Linker
C_BLUE   = "\033[94m"   # Vernier
C_RESET  = "\033[0m"


# ---- Scheme → canonical (inclusive) CDR base-number ranges ----
CDR_MAP = {
    "kabat":   {"H1": (30,35), "H2": (50,65), "H3": (95,102), "L1": (24,34), "L2": (50,56), "L3": (89,97)},
    "chothia": {"H1": (26,32), "H2": (52,56), "H3": (95,102), "L1": (26,32), "L2": (50,52), "L3": (91,96)},
    "martin":  {"H1": (30,35), "H2": (47,58), "H3": (95,101), "L1": (30,35), "L2": (46,55), "L3": (89,97)},
    "imgt":    {"H1": (27,38), "H2": (56,65), "H3": (105,117), "L1": (27,38), "L2": (56,65), "L3": (105,117)},
    "aho":     {"H1": (25,40), "H2": (58,77), "H3": (109,137), "L1": (25,40), "L2": (58,77), "L3": (109,137)},
}

# ---- 3→1 map (includes MSE→M normalization) ----
_ONE = {
    "ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G","HIS":"H",
    "ILE":"I","LYS":"K","LEU":"L","MET":"M","ASN":"N","PRO":"P","GLN":"Q",
    "ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y","MSE":"M"
}
def res3_to_1(resname: str) -> str:
    return _ONE.get((resname or "").strip().upper(), "?")

# -------------------------
# Utilities
# -------------------------
def _collapse_ranges(idxs):
    if not idxs:
        return []
    idxs = sorted(set(idxs))
    out, start = [], idxs[0]
    prev = start
    for x in idxs[1:]:
        if x == prev + 1:
            prev = x
        else:
            out.append((start, prev))
            start = prev = x
    out.append((start, prev))
    return out

def _format_ranges(ranges):
    parts = []
    for a, b in ranges:
        parts.append(str(a) if a == b else f"{a}-{b}")
    return ", ".join(parts)

def colorize(seq, nums, scheme, chain_letter):
    cdrs = CDR_MAP.get(scheme.lower(), {})
    colored = []
    for i, aa in enumerate(seq):
        num = nums[i] if i < len(nums) else None
        color = ""
        if num is not None:
            for loop, (a, b) in cdrs.items():
                if loop.startswith(chain_letter) and a <= num <= b:
                    color = C_MAGENTA if loop[0] == "H" else C_GREEN
                    break
        colored.append(f"{color}{aa}{C_RESET}")
    return "".join(colored)

# -------------------------
# Structure → linear FASTA + residue array
# -------------------------
def load_structure(struct_path):
    ext = os.path.splitext(struct_path)[1].lower()
    parser = PDBParser(QUIET=True) if ext == ".pdb" else MMCIFParser(QUIET=True)
    return parser.get_structure("struct", struct_path)

def extract_linear_from_structure(structure, chains_filter=None):
    """
    Walk selected chains (or all), collect only protein residues (res.id[0] == ' ' and is_aa),
    return a single concatenated one-letter sequence and a parallel list of Residue objects.
    """
    reslist = []
    ppb = PPBuilder()
    for model in structure:
        for chain in model:
            if chains_filter and chain.id not in chains_filter:
                continue
            # Use PPBuilder to respect peptide order and skip breaks
            peptides = ppb.build_peptides(chain)
            for pp in peptides:
                for res in pp:
                    if res.id[0] == " " and is_aa(res.get_resname(), standard=True):
                        reslist.append(res)
    seq = "".join(res3_to_1(r.get_resname()) for r in reslist)
    return seq, reslist

def write_linear_map_tsv(fname, pdb_seq, reslist):
    with open(fname, "w") as f:
        f.write("fasta_idx\tchain\tresseq\ticode\tresname\tAA\n")
        for i, r in enumerate(reslist):
            chain = r.get_parent().id
            _, resseq, icode = r.id
            f.write(f"{i}\t{chain}\t{resseq}\t{icode or ''}\t{r.get_resname().strip()}\t{pdb_seq[i]}\n")

# -------------------------
# ANARCI
# -------------------------
def run_anarci_on_seq(seq, scheme):
    tmp_fa = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta").name
    with open(tmp_fa, "w") as f:
        f.write(">linear\n")
        # wrap to avoid insanely long lines in some tools (ANARCI is fine either way)
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")
    tmp_out = tempfile.NamedTemporaryFile(delete=False, suffix=".txt").name
    cmd = ["ANARCI", "-i", tmp_fa, "--scheme", scheme, "--outfile", tmp_out]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    with open(tmp_out) as f:
        text = f.read()
    os.remove(tmp_fa)
    os.remove(tmp_out)
    return text

def run_anarci_on_fasta(path, scheme):
    tmp_out = tempfile.NamedTemporaryFile(delete=False, suffix=".txt").name
    cmd = ["ANARCI", "-i", path, "--scheme", scheme, "--outfile", tmp_out]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    with open(tmp_out) as f:
        text = f.read()
    os.remove(tmp_out)
    return text

def extract_domains(text):
    """
    Extract FIRST VH and FIRST VL(VK) domains from ANARCI output.
    Returns: (vh_seq, vh_nums), (vl_seq, vl_nums)
    nums[i] holds the base number for seq[i] (95A -> 95).
    """
    vh_seq, vh_nums = "", []
    vl_seq, vl_nums = "", []
    have_vh = have_vl = False
    cur_type = None
    buf_seq, buf_nums = [], []

    def flush():
        nonlocal have_vh, have_vl, vh_seq, vh_nums, vl_seq, vl_nums, cur_type, buf_seq, buf_nums
        if not buf_seq or cur_type is None:
            cur_type = None; buf_seq = []; buf_nums = []
            return
        if cur_type == 'H' and not have_vh:
            vh_seq, vh_nums = "".join(buf_seq), buf_nums[:]; have_vh = True
        elif cur_type in ('L','K') and not have_vl:
            vl_seq, vl_nums = "".join(buf_seq), buf_nums[:]; have_vl = True
        cur_type = None; buf_seq = []; buf_nums = []

    for raw in text.splitlines():
        line = raw.strip()
        if line == "//":
            flush(); continue
        if not line or line[0] not in "HLK":
            continue
        toks = line.split()
        if len(toks) < 3:
            continue
        c = toks[0]
        try:
            base = int(toks[1])
        except ValueError:
            continue
        aa = toks[2] if len(toks) == 3 else toks[3]
        if aa == "-":
            continue
        c = 'H' if c == 'H' else 'L'
        if cur_type is None:
            cur_type = c
        elif cur_type != c:
            flush(); cur_type = c
        buf_seq.append(aa)
        buf_nums.append(base)
    flush()
    return (vh_seq, vh_nums), (vl_seq, vl_nums)

# -------------------------
# Exact substring mapping (no aligners)
# -------------------------
def slice_by_linker(host_seq, host_res, linker):
    p = host_seq.find(linker)
    if p == -1:
        return (host_seq, host_res, host_seq, host_res, None, None)
    L = len(linker)
    return (host_seq[:p], host_res[:p], host_seq[p+L:], host_res[p+L:], p, L)

def map_domain_by_substring(dom_seq, host_seq, host_res, label=""):
    if not dom_seq:
        return {}
    start = host_seq.find(dom_seq)
    if start == -1:
        raise ValueError(f"[{label}] domain sequence not found in host sequence.")
    # (Optionally check uniqueness)
    # if host_seq.find(dom_seq, start+1) != -1:
    #     raise ValueError(f"[{label}] domain sequence occurs multiple times; refine search.")
    return {j: host_res[start + j] for j in range(len(dom_seq))}

def write_anarci_mapping_tsv(fname, label, dom_seq, mapping):
    with open(fname, "w") as f:
        f.write("label\tanarci_idx\tAA\tchain\tresseq\ticode\tresname\n")
        for i in sorted(mapping):
            r = mapping[i]
            chain = r.get_parent().id
            _, resseq, icode = r.id
            f.write(f"{label}\t{i}\t{dom_seq[i]}\t{chain}\t{resseq}\t{icode or ''}\t{r.get_resname().strip()}\n")

# -------------------------
# Vernier (NeighborSearch on structure)
# -------------------------
def compute_vernier(structure, cdr_residues, chains_filter, cutoff):
    if not cdr_residues:
        return set()
    atoms = []
    allowed = set()
    for model in structure:
        for chain in model:
            if chains_filter and chain.id not in chains_filter:
                continue
            for res in chain:
                if res.id[0] == " " and is_aa(res.get_resname(), standard=True):
                    allowed.add(res)
                    atoms.extend(res.get_atoms())
    ns = NeighborSearch(list(atoms))
    cdr_set = set(cdr_residues)
    out = set()
    for res in cdr_residues:
        for atom in res.get_atoms():
            for a in ns.search(atom.coord, cutoff + 1e-6):
                rr = a.get_parent()
                if rr in allowed and rr not in cdr_set:
                    out.add(rr)
    return out

# -------------------------
# Reporting
# -------------------------
import os, re, html
from typing import Optional, List, Set, Dict, Any

# -------------------------------------------------------------
# Numbered block generator (safe for HTML tokens)
# -------------------------------------------------------------
def make_numbered_block(
    seq_html: Optional[str],
    start: int = 1,
    wrap: int = 20,     # default wrap size
    colwidth: int = 5,
) -> str:
    if not seq_html:
        return ""

    # Remove tags to determine number of residues
    plain = re.sub(r"<[^>]+>", "", seq_html)

    # Tokenize: match exactly one residue: either a <span> tag or a character
    tokens = []
    html_rest = seq_html
    for _ in plain:
        m = re.match(r'(<span[^>]*>[^<]</span>|.)', html_rest)
        if not m:
            break
        tok = m.group(0)
        tokens.append(tok)
        html_rest = html_rest[len(tok):]

    lines = []
    n = len(tokens)

    for block_start in range(0, n, wrap):
        block_end = min(block_start + wrap, n)
        toks = tokens[block_start:block_end]

        number_line = "".join(f"{i:>{colwidth}}" for i in range(start+block_start, start+block_end))
        aa_line = "".join(" "*(colwidth-1) + tok for tok in toks)

        lines.append(number_line + "\n" + aa_line + "\n")

    return "<pre class='mono'>" + "".join(lines).rstrip() + "</pre>"

# -------------------------------------------------------------
# FULL SUMMARY GENERATION 
# -------------------------------------------------------------
def write_full_html_summary(
    *,
    base: str,
    scheme: str,
    vh: Optional[str],
    vh_nums: Optional[List[int]],
    vl: Optional[str],
    vl_nums: Optional[List[int]],
    linker: str,
    vernier_positions: Optional[Set[int]],
    fixed: Optional[Dict[str, Any]],
    input_path: str,
    save_anarci: bool,
    pdb_seq: Optional[str] = None  # <-- ADD: single combined PDB sequence OR None
) -> str:

    vernier_positions = vernier_positions or set()

    # Build complete scFv
    MA, RA = "MA", "RA"
    scfv_seq = MA + (vh or "") + (linker or "") + (vl or "") + RA

    MA_LEN, VH_LEN, LK_LEN, VL_LEN, RA_LEN = (
        len(MA), len(vh or ""), len(linker or ""), len(vl or ""), len(RA)
    )

    off_vh = MA_LEN
    off_lk = off_vh + VH_LEN
    off_vl = off_lk + LK_LEN

    total_len = len(scfv_seq)

    # CDR definitions
    cdrs = CDR_MAP.get(scheme.lower(), {})

    # ---------------------------------------------------------
    # Build HTML per residue: <span class="class">AA</span>
    # ---------------------------------------------------------
    html_tokens = []

    for i1, aa in enumerate(scfv_seq, start=1):

        classes = []
        bold = False

        # Vernier
        if i1 in vernier_positions:
            classes.append("blu"); bold = True

        # Linker
        elif off_lk + 1 <= i1 <= off_lk + LK_LEN:
            classes.append("gry"); bold = True

        # VH CDR logic
        elif off_vh + 1 <= i1 <= off_vh + VH_LEN:
            ix = i1 - off_vh - 1
            if vh_nums and 0 <= ix < len(vh_nums):
                n = vh_nums[ix]
                for loop in ("H1", "H2", "H3"):
                    a, b = cdrs.get(loop, (10**9, -10**9))
                    if a <= n <= b:
                        classes.append("mag"); bold = True
                        break

        # VL CDR logic
        elif off_vl + 1 <= i1 <= off_vl + VL_LEN:
            ix = i1 - off_vl - 1
            if vl_nums and 0 <= ix < len(vl_nums):
                n = vl_nums[ix]
                for loop in ("L1", "L2", "L3"):
                    a, b = cdrs.get(loop, (10**9, -10**9))
                    if a <= n <= b:
                        classes.append("grn"); bold = True
                        break

        # Underline VH and VL regions
        if off_vh + 1 <= i1 <= off_vh + VH_LEN:
            classes.append("vh")
        elif off_vl + 1 <= i1 <= off_vl + VL_LEN:
            classes.append("vl")

        class_str = " ".join(classes)
        style = "font-weight:bold" if bold else ""
        style_attr = f' style="{style}"' if style else ""

        html_tokens.append(
            f'<span class="{class_str}"{style_attr}>{html.escape(aa)}</span>'
        )

    scfv_html = "".join(html_tokens)
    scfv_unnumbered = scfv_html  # unnumbered block

    # ---------------------------------------------------------
    # Position sets
    # ---------------------------------------------------------
    loops_positions = set()
    linker_positions_full = set()
    framework_positions_full = set()

    if fixed:
        for key in ("H1","H2","H3","L1","L2","L3"):
            for a,b in fixed.get(key, []):
                loops_positions.update(range(a,b+1))
        for a,b in fixed.get("LINKER", []):
            linker_positions_full.update(range(a,b+1))
        for a,b in fixed.get("FRAMEWORK", []):
            framework_positions_full.update(range(a,b+1))

    vernier_full = set(vernier_positions)
    caps = {1,2,total_len-1,total_len} if total_len>=4 else set()

    fixed_positions_full = loops_positions | linker_positions_full | vernier_full | caps
    design_positions_full = set(range(1,total_len+1)) - fixed_positions_full
    lss_positions_full = loops_positions | vernier_full

    def _line(name,s):
        if not s: return f"{name}: (none)"
        rs = _format_ranges(_collapse_ranges(sorted(s)))
        return f"{name} (n={len(s)}): {rs}"

    pos_lines = "\n".join([
        _line("Loops (H1–H3 + L1–L3)", loops_positions),
        _line("Vernier", vernier_full),
        _line("LSS = Loops ∪ Vernier", lss_positions_full),
        _line("Fixed", fixed_positions_full),
        _line("Design", design_positions_full),
    ])

    # ---------------------------------------------------------
    # CSS
    # ---------------------------------------------------------
    css = """
    body { margin:24px; font-family:Inter,system-ui,-apple-system; background:#fafafa; color:#222; }
    h1 { font-size:22px; margin-bottom:6px; }
    .mono {
        font-family:"Fira Code","Source Code Pro",monospace;
        white-space:pre-wrap;
        word-wrap: break-word;
        line-height:1.45;
        font-size:15px;
        overflow-x:auto;
        max-width:100%;
    }

    .mono span { white-space: pre-wrap; }

    /* Tighter, single-spaced legend */
    .mono.legend-compact {
        white-space: normal;   /* <- stops treating source newlines as extra breaks */
        line-height: 1.2;      /* <- tighten vertical spacing */
        font-size: 13px;       /* optional, makes the legend feel more compact */
    }


    .box { background:#fff; border:1px solid #e5e7eb; padding:14px; border-radius:10px; margin:14px 0; }
    .small { font-size:12px; color:#555; }

    /* Colors */
    .mag { color:#D81B60; }
    .grn { color:#1E88E5; }
    .gry { color:#6b7280; }
    .blu { color:#43A047; }

    .vh { text-decoration:underline; }
    .vl { text-decoration:underline; }
    /* Scrollable FASTA blocks (Top, Bottom, Random, All Designs) */
    .designs-block {
        max-height: 420px;      /* adjust height as desired */
        overflow-y: auto;       /* vertical scrolling */
        overflow-x: auto;       /* horizontal scrolling for long sequences */
        border: 1px solid #ddd;
        background: #fdfdfd;
        padding: 8px;
        margin-top: 8px;
        white-space: pre;       /* preserves FASTA formatting perfectly */
        font-size: 14px;
        line-height: 1.3;
    }
    /* Scrollable ranking tables */
    .ranking-block {
        max-height: 420px;      /* same as FASTA block */
        overflow-y: auto;
        overflow-x: auto;
        border: 1px solid #ddd;
        background: #fafafa;
        padding: 8px;
        margin-top: 8px;
        white-space: pre;       /* keeps alignment intact */
        font-size: 14px;
        line-height: 1.3;
    }
    /* Ensure table respects scroll container */
    .ranking-block table {
        display: block;
        width: max-content;   /* prevents forced stretching */
    }



    """

    def box(title, content):
        return f'<div class="box"><div class="small"><b>{html.escape(title)}</b></div>{content}</div>'

    # Numbered scFv
    scfv_block = make_numbered_block(scfv_html, start=1, wrap=20, colwidth=5)

    # Position sets block
    pos_block = f"<pre class='mono'>{html.escape(pos_lines)}</pre>"

    # ---------------------------------------------------------
    # Build FASTA block from structure input, if provided
    # ---------------------------------------------------------
    pdb_fasta_block = ""
    if pdb_seq:
        pdb_name = os.path.basename(input_path)
        pdb_fasta_block = (
            f"<div class='mono'>>{html.escape(pdb_name)}<br>{html.escape(pdb_seq)}</div>"
)


    # ---------------------------------------------------------
    # HTML Document
    # ---------------------------------------------------------
    html_doc = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>scFv summary — {html.escape(base)}</title>
<style>{css}</style>
</head>
<body>

<h1>scFv summary — {html.escape(base)}</h1>

<div class="box">
  <div class="small"><b>Legend</b>
  </div>
  <div class="mono legend-compact">
      <span class="mag" style="font-weight:bold">H-CDR</span> – Heavy-chain CDR residues<br>
      <span class="grn" style="font-weight:bold">L-CDR</span> – Light-chain CDR residues<br>
      <span class="gry" style="font-weight:bold">Linker</span> – Linker residues<br>
      <span class="blu" style="font-weight:bold">Vernier</span> – Vernier positions<br>
      <span class="vh">VH region</span> – Underlined<br>
      <span class="vl">VL region</span> – Underlined<br>
      <b>MA</b> and <b>RA</b> – flanking tags (not underlined)<br>
  </div>
</div>

{box("FASTA sequence derived from input", pdb_fasta_block)}

<div class="box">
  <div class="small"><b>Full scFv sequence (formatted, no numbers)</b></div>
  <div class="mono">{{SCFV_UNNUMBERED}}</div>
</div>

{box("Assembled scFv (MA–VH–LINKER–VL–RA)", scfv_block)}
{box("Position sets (scFv coordinates)", pos_block)}

</body>
</html>
"""

    # Replace placeholders
    html_doc = html_doc.replace("{SCFV_UNNUMBERED}", scfv_unnumbered)

    out_path = f"{base}_summary.html"
    with open(out_path,"w",encoding="utf-8") as f:
        f.write(html_doc)

    return out_path




# -------------------------
# Main
# -------------------------
def main():
    ap = argparse.ArgumentParser(
        description="ANARCI VH/VL → scFv with deterministic PDB↔FASTA mapping; optional Vernier shell."
    )
    ap.add_argument("input", help="Input FASTA or structure (.pdb/.cif/.mmcif)")
    ap.add_argument("--scheme", default="imgt", help="imgt, kabat, martin, chothia, aho")
    ap.add_argument("--combine", action="store_true", help="Write scFv = MA + VH + linker + VL + RA")
    ap.add_argument("--linker", default="GGGGSGGGGSGGGGS", help="Linker used between VH and VL")
    ap.add_argument("--outfile", default="scfv_output.fasta", help="FASTA output")
    ap.add_argument("--save-anarci", action="store_true", help="Save raw ANARCI text")
    ap.add_argument("--chains", help="Comma-separated chain IDs to include (structure input only)")
    ap.add_argument("--keep-linear-tsv", action="store_true", help="Write pdb_linear_map.tsv")
    ap.add_argument("--dump-mapping", action="store_true", help="Write VH_mapping.tsv / VL_mapping.tsv")
    # Vernier control
    ap.add_argument("--vernier-cutoff", type=float, default=None,
                    help="Å cutoff for Vernier (structure input only)")
    ap.add_argument("--dist", type=float, default=None,
                    help="Alias for --vernier-cutoff (loops_from_sequence.py compatibility)")
    # Simple-output mode (suppress screen printing)
    ap.add_argument(
        "--simple-output",
        choices=["loops", "lss", "framework","vernier","fixed","design"],
        help="Print ONLY space-separated design positions (loops_from_sequence.py compatible)"
    )
    args = ap.parse_args()

    # Backward compatibility
    if args.vernier_cutoff is None and args.dist is not None:
        args.vernier_cutoff = args.dist

    infile = args.input
    ext = os.path.splitext(infile)[1].lower()
    chains_filter = set(args.chains.split(",")) if args.chains else None

    # ----------------------------------------------------
    # Helper to conditionally print (silent in simple mode)
    # ----------------------------------------------------
    def screen_print(*a, **k):
        if not args.simple_output:
            print(*a, **k)

    # --- Input handling ---
    structure = None
    if ext in (".pdb", ".cif", ".mmcif"):
        structure = load_structure(infile)
        pdb_seq, pdb_res = extract_linear_from_structure(structure, chains_filter)
        if not pdb_seq:
            print("No protein residues found.")
            sys.exit(1)

        if args.keep_linear_tsv:
            write_linear_map_tsv("pdb_linear_map.tsv", pdb_seq, pdb_res)

        screen_print(f"# Structure loaded → linear protein length: {len(pdb_seq)}")
        anarci_text = run_anarci_on_seq(pdb_seq, args.scheme)

    else:
        with open(infile) as f:
            raw = f.read()
        seq = "".join(line.strip() for line in raw.splitlines() if not line.startswith(">"))
        pdb_seq, pdb_res = seq, []
        anarci_text = run_anarci_on_fasta(infile, args.scheme)

    if args.save_anarci:
        base = os.path.splitext(os.path.basename(infile))[0]
        with open(f"{base}_anarci.txt", "w") as f:
            f.write(anarci_text)
        screen_print(f"# ANARCI output saved to {base}_anarci.txt")

    # --- Extract VH/VL ---
    (vh, vh_nums), (vl, vl_nums) = extract_domains(anarci_text)
    if not vh and not vl:
        print("No VH/VL domains found.")
        sys.exit(0)

    # --- Write FASTA (always) ---
    with open(args.outfile, "w") as f:
        if vh: f.write(">VH\n" + vh + "\n")
        if vl: f.write(">VL\n" + vl + "\n")
        if args.combine and vh and vl:
            scfv = "MA" + vh + args.linker + vl + "RA"
            f.write(">scFv\n" + scfv + "\n")
    screen_print(f"# FASTA written to {args.outfile}")

    # --- ANARCI→PDB mapping ---
    vh_map = vl_map = {}
    vernier_positions = set()

    if structure:
        seq_left, res_left, seq_right, res_right, pos, L = slice_by_linker(pdb_seq, pdb_res, args.linker)

        if vh:
            try:    vh_map = map_domain_by_substring(vh, seq_left, res_left, "VH")
            except: vh_map = map_domain_by_substring(vh, pdb_seq, pdb_res, "VH")

        if vl:
            try:    vl_map = map_domain_by_substring(vl, seq_right, res_right, "VL")
            except: vl_map = map_domain_by_substring(vl, pdb_seq, pdb_res, "VL")

        if args.dump_mapping:
            if vh_map:
                write_anarci_mapping_tsv("VH_mapping.tsv", "VH", vh, vh_map)
                screen_print("# VH_mapping.tsv written")
            if vl_map:
                write_anarci_mapping_tsv("VL_mapping.tsv", "VL", vl, vl_map)
                screen_print("# VL_mapping.tsv written")

    # --- Compute CDR / Framework / Vernier (if combine) ---
    if args.combine and vh and vl:
        scheme_key = args.scheme.lower()
        if scheme_key not in CDR_MAP:
            print(f"# Unknown scheme '{args.scheme}'")
            return

        cdrs = CDR_MAP[scheme_key]
        MA_LEN, RA_LEN = 2, 2
        vh_len, vl_len, link_len = len(vh), len(vl), len(args.linker)

        vh_off = MA_LEN
        linker_off = vh_off + vh_len
        vl_off = linker_off + link_len
        total_len = MA_LEN + vh_len + link_len + vl_len + RA_LEN

        # --- Build CDR positions ---
        fixed = {}
        for loop in ("H1","H2","H3"):
            a, b = cdrs[loop]
            idxs = [vh_off + i + 1 for i,n in enumerate(vh_nums) if a <= n <= b]
            fixed[loop] = _collapse_ranges(idxs)
        for loop in ("L1","L2","L3"):
            a, b = cdrs[loop]
            idxs = [vl_off + i + 1 for i,n in enumerate(vl_nums) if a <= n <= b]
            fixed[loop] = _collapse_ranges(idxs)

        linker_positions = list(range(linker_off + 1, linker_off + link_len + 1))
        fixed["LINKER"] = _collapse_ranges(linker_positions)

        excluded = set(linker_positions)
        for loop in ("H1","H2","H3","L1","L2","L3"):
            for a,b in fixed[loop]:
                excluded.update(range(a, b+1))
        excluded.update(range(1, MA_LEN+1))
        excluded.update(range(total_len-RA_LEN+1, total_len+1))

        framework_positions = [pos for pos in range(1, total_len+1) if pos not in excluded]
        fixed["FRAMEWORK"] = _collapse_ranges(framework_positions)

        # --- Write fixed positions file (always) ---
        base = os.path.splitext(os.path.basename(args.input))[0]
        fixed_out = f"{base}_fixed_positions.txt"

        # --- Report CDRs, LINKER, FRAMEWORK to both file and screen ---
        lines = []
        lines.append("# scFv fixed positions (1-based)")
        lines.append(f"# Scheme: {args.scheme}")

        for key in ("H1","H2","H3","L1","L2","L3","LINKER","FRAMEWORK"):
            line = f"{key}: {_format_ranges(fixed[key])}"
            lines.append(line)
            screen_print(line)   # <-- print to screen as well

        # Write file
        with open(fixed_out, "w") as f:
            f.write("\n".join(lines) + "\n")

        screen_print(f"# Fixed positions written to {fixed_out}")


        # --- Vernier computation ---
        if args.vernier_cutoff and structure and (vh_map or vl_map):
            cdr_residues = []
            def add_from_map(dom_nums, dom_map, which_loops):
                for loop in which_loops:
                    a,b = cdrs[loop]
                    for i,n in enumerate(dom_nums):
                        if a <= n <= b and i in dom_map:
                            cdr_residues.append(dom_map[i])
            if vh_map:
                add_from_map(vh_nums, vh_map, ("H1","H2","H3"))
            if vl_map:
                add_from_map(vl_nums, vl_map, ("L1","L2","L3"))

            vernier = compute_vernier(structure, cdr_residues, chains_filter, args.vernier_cutoff)
            rev_VH = {res:i for i,res in vh_map.items()} if vh_map else {}
            rev_VL = {res:i for i,res in vl_map.items()} if vl_map else {}

            for res in vernier:
                if res in rev_VH:
                    vernier_positions.add(MA_LEN + rev_VH[res] + 1)
                elif res in rev_VL:
                    vernier_positions.add(MA_LEN + vh_len + link_len + rev_VL[res] + 1)

            vernier_ranges = _collapse_ranges(sorted(vernier_positions))
            vr_line = f"VERNIER({args.vernier_cutoff:.1f}Å): {_format_ranges(vernier_ranges)}"
            with open(fixed_out, "a") as f:
                f.write(vr_line + "\n")
            screen_print(vr_line)

    # --- Final colored preview (screen only) ---
    if not args.simple_output:
        screen_print("\n# Colored previews (final)")
        if vh:
            screen_print(">VH (colored CDRs)")
            screen_print(colorize(vh, vh_nums, args.scheme, "H"), "\n")
        if vl:
            screen_print(">VL (colored CDRs)")
            screen_print(colorize(vl, vl_nums, args.scheme, "L"), "\n")

        if args.combine and vh and vl:
            scheme_key = args.scheme.lower()
            cdrs = CDR_MAP.get(scheme_key, {})
            MA, RA = "MA", "RA"
            scfv_seq = MA + vh + args.linker + vl + RA
            MA_LEN, VH_LEN, LK_LEN, VL_LEN = len(MA), len(vh), len(args.linker), len(vl)
            off_vh, off_lk, off_vl = MA_LEN, MA_LEN + VH_LEN, MA_LEN + VH_LEN + LK_LEN

            painted = []
            for i1,aa in enumerate(scfv_seq, start=1):
                color = ""
                if i1 in vernier_positions:
                    color = C_BLUE
                elif off_lk + 1 <= i1 <= off_lk + LK_LEN:
                    color = C_GRAY
                elif off_vh + 1 <= i1 <= off_vh + VH_LEN:
                    ix = i1 - off_vh - 1
                    if 0 <= ix < len(vh_nums):
                        n = vh_nums[ix]
                        for loop in ("H1","H2","H3"):
                            a,b = cdrs.get(loop, (10**9, -10**9))
                            if a <= n <= b:
                                color = C_MAGENTA
                                break
                elif off_vl + 1 <= i1 <= off_vl + VL_LEN:
                    ix = i1 - off_vl - 1
                    if 0 <= ix < len(vl_nums):
                        n = vl_nums[ix]
                        for loop in ("L1","L2","L3"):
                            a,b = cdrs.get(loop, (10**9, -10**9))
                            if a <= n <= b:
                                color = C_GREEN
                                break
                painted.append(f"{color}{aa}{C_RESET}")

            legend = f"{C_MAGENTA}H-CDRs{C_RESET}, {C_GREEN}L-CDRs{C_RESET}, {C_GRAY}LINKER{C_RESET}, {C_BLUE}VERNIER{C_RESET}"
            screen_print(">scFv (colored)")
            screen_print(f"# Legend: {legend}\n")

            WRAP=100
            line=[]
            for k,token in enumerate(painted, start=1):
                line.append(token)
                if k % WRAP == 0:
                    screen_print("".join(line))
                    line=[]
            if line: screen_print("".join(line))
            screen_print()

    # --- HTML summary (always written) ---
    base = os.path.splitext(os.path.basename(args.input))[0]
    vernier_positions_local = locals().get('vernier_positions', set())
    fixed_local = locals().get('fixed', {})

    write_full_html_summary(
        base=base,
        scheme=args.scheme,
        vh=vh, vh_nums=vh_nums,
        vl=vl, vl_nums=vl_nums,
        linker=args.linker,
        vernier_positions=vernier_positions_local,
        fixed=fixed_local,
        input_path=args.input,
        save_anarci=args.save_anarci,
        pdb_seq=pdb_seq,
    )

    # === SIMPLE OUTPUT MODE (final print only) ===
    if args.simple_output:

        # Flatten helper for ranges like [(a,b),(c,d)]
        def _expand(ranges):
            return [i for a, b in ranges for i in range(a, b + 1)]

        # Core loop definitions
        loops_positions = sorted(
            _expand(fixed["H1"]) +
            _expand(fixed["H2"]) +
            _expand(fixed["H3"]) +
            _expand(fixed["L1"]) +
            _expand(fixed["L2"]) +
            _expand(fixed["L3"])
        )

        # Linker
        linker_positions_full = _expand(fixed["LINKER"])

        # Vernier
        vernier_positions_full = sorted(vernier_positions)

        # Include MA (1–2) and RA (last 2 residues)
        MA_RA_positions = [1, 2, total_len - 1, total_len]

        # Build FIXED = loops + linker + vernier + MA/RA
        fixed_positions_full = sorted(
            set(loops_positions)
            | set(linker_positions_full)
            | set(vernier_positions_full)
            | set(MA_RA_positions)
        )

        # DESIGN = everything NOT FIXED
        design_positions_full = sorted(
            set(range(1, total_len + 1)) - set(fixed_positions_full)
        )

        # --- Output modes ---
        mode = args.simple_output

        if mode == "loops":
            out = loops_positions

        elif mode == "lss":  # loops + secondary shell
            out = sorted(set(loops_positions) | set(vernier_positions_full))

        elif mode == "vernier":  # secondary shell only
            out = vernier_positions_full

        elif mode == "fixed":  # loops + linker + vernier + MA/RA
            out = fixed_positions_full

        elif mode == "design":  # everything not fixed
            out = design_positions_full

        elif mode == "framework":
            out = framework_positions

        print(" ".join(str(i) for i in out))
        return




if __name__ == "__main__":
    main()