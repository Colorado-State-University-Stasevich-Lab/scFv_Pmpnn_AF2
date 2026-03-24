#!/usr/bin/env python3
import sys, re, html
from pathlib import Path

###############################################################
# Appends colorized design sequences + mutation highlighting
# to existing scFv summary HTML.
#
# Mutation logic:
#   - Compare each design to WT
#   - Any difference → wrap in <span class="mutbg"> ... </span>
#
# This keeps original color-coding (mag, grn, blu, vh, vl, etc.)
# but adds a gray background to visually mark mutated sites.
###############################################################


def load_combined_fasta(path):
    entries = []
    header = None
    seq = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    entries.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header:
            entries.append((header, "".join(seq)))

    return entries


def extract_colorized_block(html_text):
    m = re.search(r"<div class=\"mono\">(.+?)</div>", html_text, flags=re.DOTALL)
    if not m:
        print("ERROR: cannot find scFv colored HTML block in summary.")
        sys.exit(1)
    return m.group(1).strip()


def tokenize_html_sequence(html_seq):
    pattern = re.compile(r'(<span.*?>.*?</span>|.)', re.DOTALL)
    tokens = []
    for m in pattern.finditer(html_seq):
        tok = m.group(1)
        if tok.strip():
            tokens.append(tok)
    return tokens


def colorize_with_mutations(template_tokens, new_seq, wt_seq):
    """
    Copy coloring from template for positions within scFv length.
    For extra residues (e.g., linker, epitope), just output plain letters.
    Highlight mutated residues vs WT using <span class="mutbg">...</span>.
    """
    out = []
    template_len = len(template_tokens)
    wt_len = len(wt_seq)

    for i, aa in enumerate(new_seq):

        # Case 1 — within template region → apply original coloring
        if i < template_len and i < wt_len:
            base_tok = template_tokens[i]
            is_mut = (aa != wt_seq[i])

            if base_tok.startswith("<span"):
                colored = re.sub(r">(.)</span>", f">{html.escape(aa)}</span>", base_tok)
            else:
                colored = html.escape(aa)

            if is_mut:
                colored = f'<span class="mutbg">{colored}</span>'

        # Case 2 — beyond template (extra residues)
        else:
            colored = html.escape(aa)

        out.append(colored)

    return "".join(out)



def append_designs_block(summary_path, fasta_path, output_path=None):
    html_text = Path(summary_path).read_text()
    entries = load_combined_fasta(fasta_path)

    # WT = first entry
    wt_name, wt_seq = entries[0]

    # Load original colored scFv template
    original_block = extract_colorized_block(html_text)
    template_tokens = tokenize_html_sequence(original_block)

    # Insert CSS for mutation highlighting AND FASTA-tight formatting
    css_block = """
    <style>
    .mutbg {
        background-color: #e0e0e0;
        border-radius: 3px;
    }
    .mono.designs-block {
        line-height: 1.0;
    }
    </style>
    """
    html_text = html_text.replace("</head>", css_block + "\n</head>")


    # Build unified FASTA-style colored block without ANY block-level wrappers
    lines = []
    for header, seq in entries:
        colored = colorize_with_mutations(template_tokens, seq, wt_seq)
        lines.append(f"&gt;{html.escape(header)}<br>{colored}<br>")


    def header_for(filename):
        if "top" in filename.lower():
            return "Top Designs (FASTA-style, aligned, colored, mutation-highlighted)"
        if "bottom" in filename.lower():
            return "Bottom Designs (FASTA-style, aligned, colored, mutation-highlighted)"
        if "random" in filename.lower():
            return "Random Designs (FASTA-style, aligned, colored, mutation-highlighted)"
        return "All Designs (FASTA-style, aligned, colored, mutation-highlighted)"

    # Determine appropriate header label based on FASTA file name
    fasta_filename = Path(fasta_path).name
    label = header_for(fasta_filename)

    # Build the block with correct heading
    designs_html = (
        "<div class='box'>"
        f"<div class='small'><b>{label}</b></div>"
        "<div class='mono designs-block'>"
        + "".join(lines) +
        "</div></div>"
    )



    # Append before </body>
    new_html = html_text.replace("</body>", designs_html + "\n</body>")

    outpath = output_path or summary_path
    Path(outpath).write_text(new_html)
    print(f"Appended designs + mutation highlighting → {outpath}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: append_designs_to_summary.py summary.html combined.fa [output.html]")
        sys.exit(1)

    summary = sys.argv[1]
    fasta = sys.argv[2]
    out = sys.argv[3] if len(sys.argv) > 3 else None

    append_designs_block(summary, fasta, out)
