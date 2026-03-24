import gradio as gr
import re
import subprocess
import tempfile
import os
from pathlib import Path
import shutil
from typing import Generator, Tuple, Optional

# Path to your pipeline script
PIPELINE_SCRIPT = str((Path(__file__).resolve().parent / "af2_pmpnnV7.sh"))


def _normalize_epitope_chain(value) -> str:
    """Return epitope chain(s) as a comma-separated string or "None".
    Accepts a single chain like "B", a comma/space-separated string like "B,C,F",
    or a list of selections from a multiselect dropdown.
    """
    if value is None:
        return "None"
    # Gradio multiselect returns list[str]
    if isinstance(value, (list, tuple)):
        items = [str(v).strip() for v in value if str(v).strip()]
        items = [v for v in items if v.lower() != "none"]
        if not items:
            return "None"
        return ",".join(items)
    s = str(value).strip()
    if not s or s.lower() == "none":
        return "None"
    # Allow "B,C,F" or "B C F"
    parts = re.split(r"[\s,;]+", s)
    parts = [p.strip() for p in parts if p.strip()]
    parts = [p for p in parts if p.lower() != "none"]
    return ",".join(parts) if parts else "None"


def run_pipeline(
    input_file,
    determine_cdrs: str,
    ss_near_cdrs: int,
    linker_seq: str,
    seqs_per_run: int,
    epitope_chain: str,
    gpu_choice: str,
    skip_folding: bool,
    N_select: int,
    pmpnn_seed: int,     
    run_AF2: bool,
) -> Generator[Tuple[str, Optional[str], Optional[str]], None, None]:
    """
    Generator for Gradio: streams log text as the pipeline runs,
    then returns (log, zip_file_path, html_preview) at the end.

    Outputs:
      - log_text (str)
      - zip_file_path (str or None)
      - html_preview (str or None, HTML)
    """

    # Create isolated workdir for this run
    workdir = Path(tempfile.mkdtemp())
    print(f"[PIPELINE] Temporary run directory: {workdir}", flush=True)
    input_dir = workdir / "input_dir"
    output_dir = workdir / "output_dir"
    input_dir.mkdir(parents=True, exist_ok=True)

    # Copy scripts directory into the working directory
    script_src = Path(__file__).resolve().parent / "scripts"
    script_dst = workdir / "scripts"
    shutil.copytree(script_src, script_dst)

    # Ensure a clean output_dir
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save uploaded file into input_dir
    infile_path = input_dir / Path(input_file.name).name
    shutil.copy(input_file.name, infile_path)

    ext = infile_path.suffix.lower()
    is_fasta = ext in [".fa", ".fasta"]

    log = []

    def emit(
        log_line: str = "",
        zip_path: Optional[Path] = None,
        html_preview: Optional[str] = None,
    ):
        """Helper to join and yield log text + optional outputs."""
        if log_line:
            log.append(log_line)
        full_log = "".join(log)
        return full_log, str(zip_path) if zip_path else None, html_preview

    # Build environment for the bash script
    env = os.environ.copy()
    env["input_dir"] = str(input_dir)
    env["output_dir"] = str(output_dir)

    # Core design / scFv settings
    env["determine_CDRs"] = determine_cdrs
    env["ss_near_CDRs"] = str(ss_near_cdrs)
    env["linker_seq"] = linker_seq
    env["seqs_per_run"] = str(seqs_per_run)
    #env["epitope_chain"] = _normalize_epitope_chain(epitope_chain)
    if is_fasta:
        env["epitope_chain"] = "None"
    else:
        env["epitope_chain"] = _normalize_epitope_chain(epitope_chain)

    # Folding controls
    env["skip_folding"] = "true" if skip_folding else "false"
    env["run_AF2"] = "true" if run_AF2 else "false"

    # ProteinMPNN seed (for reproducible design runs)
    env["pmpnn_seed"] = str(pmpnn_seed)

    # Respect GPU choice for anything CUDA launched by the shell script (ColabFold/JAX)
    if gpu_choice and gpu_choice != "auto":
        env["CUDA_VISIBLE_DEVICES"] = gpu_choice
    else:
        # Make 'auto' mean "no restriction", even if parent process had it set
        env.pop("CUDA_VISIBLE_DEVICES", None)


    # ----- N selection logic: enforce Nt = Nb = Nr ≤ seqs_per_run -----
    requested_N = int(N_select)
    max_allowed = int(seqs_per_run)


    if requested_N > max_allowed:
        requested_N = max_allowed
        log.append(
            f"[WARN] Requested N_select={N_select} but only "
            f"{seqs_per_run} designs per run are generated. "
            f"Clamping to N={requested_N}.\n"
        )
    if requested_N < 1:
        requested_N = 1
        log.append("[WARN] Requested N_select < 1; clamping to N=1.\n")

    # Warn about computational cost if AF2 is enabled and N is large
    if run_AF2 and requested_N > 20:
        log.append(
            "[WARN] You requested N_select="
            f"{requested_N} with run_AF2=True. Folding more than ~20 "
            "designs in each of Top / Bottom / Random can be "
            "computationally expensive. Proceeding anyway.\n"
        )

    env["Nt"] = str(requested_N)
    env["Nb"] = str(requested_N)
    env["Nr"] = str(requested_N)

    log.append(
        f"Using N={requested_N} for Top, Bottom, and Random selections "
        f"(seqs_per_run={seqs_per_run}).\n"
    )

    # Initial log
    log.append(f"Working directory: {workdir}\n")
    log.append(f"Input file: {infile_path}\n")
    log.append("Starting pipeline...\n\n")
    yield emit()
    log.append(
        f"GPU selection: {gpu_choice} (CUDA_VISIBLE_DEVICES={env.get('CUDA_VISIBLE_DEVICES','<unset>')})\n"
    )
    yield emit()


    # Start the pipeline as a subprocess
    process = subprocess.Popen(
        ["bash", PIPELINE_SCRIPT],
        cwd=workdir,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    # Stream stdout line-by-line
    if process.stdout is not None:
        for line in process.stdout:
            log.append(line)
            # Stream incremental updates (no zip or HTML yet)
            yield emit()

    # Wait for completion
    retcode = process.wait()
    log.append(f"\nPipeline finished with exit code {retcode}\n")

    # If failure, just return log and no files
    if retcode != 0:
        log.append("\nERROR: Pipeline failed. Check log above for details.\n")
        yield emit()
        return

    # Zip up the output directory
    zip_base = workdir / "results"
    shutil.make_archive(str(zip_base), "zip", output_dir)
    zip_path = zip_base.with_suffix(".zip")

    # Try to locate a summary HTML to preview
    html_preview: Optional[str] = None
    summary_dir = output_dir / "Summary"
    if summary_dir.exists():
        # Prefer *_summary.html if multiple
        html_files = sorted(summary_dir.glob("*.html"))
        preferred: Optional[Path] = None
        for f in html_files:
            if "summary" in f.name.lower():
                preferred = f
                break
        if preferred is None and html_files:
            preferred = html_files[0]

        if preferred is not None:
            try:
                html_text = preferred.read_text()
                # Simple wrapper so it displays cleanly inside Gradio
                html_preview = (
                    "<div style='max-height:600px; overflow:auto; "
                    "border:1px solid #ccc; padding:8px; background:white;'>"
                    + html_text
                    + "</div>"
                )
                log.append(f"\nSummary HTML preview loaded from: {preferred}\n")
            except Exception as e:
                log.append(f"\nWARNING: Could not read summary HTML: {e}\n")

    # Final yield with log + zip + HTML preview (if any)
    yield emit(zip_path=zip_path, html_preview=html_preview)


# ---------------------- Gradio Interface ---------------------- #

description_md = """
### **How to use**

1. Upload a **PDB/CIF structure** or **FASTA sequence**.

   - **PDB/CIF:** scFv should be chain **A**; specify the epitope chain (e.g. **B**).
   - **FASTA:** provide the scFv sequence, optionally with an epitope:
     ```
     SCFV_SEQUENCE:EPITOPE_SEQUENCE
     ```

2. FASTA inputs are **automatically folded with AlphaFold**.

3. Adjust options and click **Run pipeline**.

4. Monitor the live log. When finished, download **results.zip** or view the **Summary HTML**.
"""

with gr.Blocks(title="scFv Intrabody Design Pipeline") as demo:
    gr.Markdown("# scFv Design Pipeline")
    gr.Markdown(description_md)

    with gr.Row():
        with gr.Column(scale=1):
            input_file = gr.File(
                label="Upload input structure (PDB or model)",
                file_types=[".pdb", ".cif", ".mmcif", ".fasta", ".fa"],
                file_count="single",
            )

            determine_cdrs = gr.Dropdown(
                choices=["martin", "chothia", "kabat", "IMGT"],
                value="martin",
                label="CDR numbering scheme",
                info="Used by make_scfv_vernierV3.py to define CDR loops.",
            )

            ss_near_cdrs = gr.Slider(
                minimum=0,
                maximum=10,
                step=1,
                value=3,
                label="Protect distance near CDRs (Å)",
                info="Positions within this distance of CDRs are NOT mutated.",
            )

            linker_seq = gr.Textbox(
                value="GGGGSGGGGSGGGGS",
                label="Linker sequence",
                info="Peptide linker between VH and VL in the combined scFv.",
            )

            seqs_per_run = gr.Slider(
                minimum=1,
                maximum=1000,
                step=1,
                value=500,
                label="Sequences per run (ProteinMPNN)",
                info="Number of designs ProteinMPNN will generate for the scFv.",
            )

            epitope_chain = gr.Textbox(
                value="B",
                label="Epitope chain ID(s)",
                info=(
                    "One or more chain IDs from the input PDB. "
                    "Examples: 'B' or 'B,C,F'. Use 'None' for scFv fasta inputs, even with epitope sequence present in fasta file."
                    ),
            )

            gpu_choice = gr.Dropdown(
                choices=["auto", "0", "1", "2", "3", "4", "5", "6"],
                value="auto",
                label="GPU",
                info="Choose which GPU ColabFold should use. 'auto' does not restrict visibility. Should be 0 or 1 for sleet.",
            )
            
            skip_folding = gr.Checkbox(
                value=False,  # matches skip_folding="${skip_folding:-true}"
                label="Skip AF2 scFv folding",
                info="Enable if your input structure is already a valid scFv model and you do not want to re-fold it.",
            )

            N_select = gr.Slider(
                minimum=1,
                maximum=50,
                step=1,
                value=20,
                label="N for Top / Bottom / Random",
                info=(
                    "Number of designs in each category (Top, Bottom, Random). "
                    "Must be ≤ sequences per run. When AF2 is enabled, "
                    "N ≲ 20 is strongly recommended due to computational cost."
                ),
            )

            pmpnn_seed = gr.Number(
                value=37,
                precision=0,
                label="ProteinMPNN random seed",
                info="Integer seed for ProteinMPNN; use the same value to reproduce a design run.",
            )

            run_AF2 = gr.Checkbox(
                value=False,
                label="Run AF2 on WT + designs",
                info="If checked, fold WT + designs with AF2; otherwise skip.",
            )

            run_btn = gr.Button("Run pipeline", variant="primary")

        with gr.Column(scale=1):
            log_box = gr.Textbox(
                label="Pipeline log (live)",
                lines=25,
                value="",
                interactive=False,
            )

            results_zip = gr.File(
                label="Download all results (results.zip)",
                interactive=False,
            )

            summary_html = gr.HTML(
                label="Summary HTML preview",
                value="",
            )

    # Wire up the button with streaming output
    run_btn.click(
        fn=run_pipeline,
        inputs=[
            input_file,
            determine_cdrs,
            ss_near_cdrs,
            linker_seq,
            seqs_per_run,
            epitope_chain,
            gpu_choice,
            skip_folding,
            N_select,
            pmpnn_seed,   
            run_AF2,
        ],
        outputs=[log_box, results_zip, summary_html],
    )

if __name__ == "__main__":
    demo.launch(server_name="0.0.0.0", server_port=7860)
