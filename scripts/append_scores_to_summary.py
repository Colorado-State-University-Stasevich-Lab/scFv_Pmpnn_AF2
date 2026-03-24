#!/usr/bin/env python3
import pandas as pd
import numpy as np
import html
import argparse
from pathlib import Path

###############################################################################
# COLOR HELPERS
###############################################################################

def color_scale(value, vmin, vmax, invert=False):
    if pd.isna(value):
        return "#ffffff"

    if abs(vmax - vmin) < 1e-9:
        return "rgba(0,255,0,0.35)"  # All same → green

    x = (value - vmin) / (vmax - vmin + 1e-9)
    x = max(0.0, min(1.0, x))
    if invert:
        x = 1.0 - x

    r = int(255 * (1 - x))
    g = int(255 * x)
    return f"rgba({r},{g},0,0.35)"


###############################################################################
# MAIN TABLE FORMATTER
###############################################################################

def dataframe_to_colored_html(df):
    """
    Flexible: If AF2 columns exist → include them.
    Otherwise show seq-only metrics.
    """

    # if "scfvtools_score" in df.columns:
    #     df = df.sort_values("scfvtools_score", ascending=False)

    preferred_order = [
        "Accession",
        "pLDDT",
        "pAE_mean",
        "pTM",
        "ipTM",
        "scfvtools_score",        # moved up
        "scfvtools_blosum_score",
        "scfvtools_blosum_diff_score",
        "Prob. of Solubility",    # moved down
        "log_likelihood_target",
        "log_likelihood"
    ]



    # Keep only columns present
    cols = [c for c in preferred_order if c in df.columns]
    df = df[cols]

    numeric_cols = df.select_dtypes(include=[np.number]).columns
    col_ranges = {
        col: (df[col].min(), df[col].max()) for col in numeric_cols
    }

    # Build HTML
    rows = []

    header = "".join(
        f"<th style='padding:4px 8px; border-bottom:1px solid #999'>{html.escape(col)}</th>"
        for col in df.columns
    )
    rows.append(f"<tr>{header}</tr>")

    for _, row in df.iterrows():
        cells = []
        for col, value in row.items():
            text = html.escape(str(value))
            if col in numeric_cols:
                inv = ("pae" in col.lower())
                vmin, vmax = col_ranges[col]
                bg = color_scale(value, vmin, vmax, invert=inv)
                cells.append(f"<td style='padding:4px 8px; background:{bg}'>{text}</td>")
            else:
                cells.append(f"<td style='padding:4px 8px'>{text}</td>")
        rows.append(f"<tr>{''.join(cells)}</tr>")

    return (
        "<table style='border-collapse:collapse; font-family:monospace; font-size:12px;'>"
        + "".join(rows) +
        "</table>"
    )


###############################################################################
# SMALL TABLE FORMATTER FOR TOP / BOTTOM / RANDOM
###############################################################################

def simple_table(csv_path, title):
    df = pd.read_csv(csv_path)
    return (
        f"<h2>{html.escape(title)}</h2>\n"
        + df.to_html(index=False, escape=True)
    )


###############################################################################
# MAIN APPENDER
###############################################################################

def append_scores_to_summary(main_csv, summary_html, out=None,
                             top_csv=None, bottom_csv=None, random_csv=None):
    df = pd.read_csv(main_csv)
    main_table = dataframe_to_colored_html(df)

    block = (
        "<div class='box'>"
        "<div class='small'><b>Design Ranking Scores</b></div>"
        "<div class='mono ranking-block'>"
        f"{main_table}"
        "</div>"
        "</div>\n"
    )


    # # Add stacked Top / Bottom / Random if present
    # if top_csv and Path(top_csv).exists():
    #     block += simple_table(top_csv, "Top Designs") + "\n"

    # if bottom_csv and Path(bottom_csv).exists():
    #     block += simple_table(bottom_csv, "Bottom Designs") + "\n"

    # if random_csv and Path(random_csv).exists():
    #     block += simple_table(random_csv, "Random Designs") + "\n"

    # Inject block before </body>
    html_text = Path(summary_html).read_text()
    new_html = html_text.replace("</body>", block + "</body>")

    out_path = out or summary_html
    Path(out_path).write_text(new_html)

    print(f"[append_scores_to_summary] Added tables → {out_path}")


###############################################################################
# CLI
###############################################################################

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="Main score table (AF2 merged or pre-AF2)")
    ap.add_argument("--summary", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--top")
    ap.add_argument("--bottom")
    ap.add_argument("--random")

    args = ap.parse_args()

    append_scores_to_summary(
        args.csv, args.summary, args.out,
        top_csv=args.top,
        bottom_csv=args.bottom,
        random_csv=args.random
    )
