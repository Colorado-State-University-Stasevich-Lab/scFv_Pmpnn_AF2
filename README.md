![scFv Colored by 3 regions](images/F1.large.jpg "Graphical Abstract Intro Code")

# **AI-assisted protein design to rapidly convert antibody sequences to intrabodies targeting diverse peptides and histone modifications**

This repository is based on the original `scFv_Pmpnn_AF2` workflow and extends it with an updated end-to-end shell pipeline, richer scoring and ranking, HTML summaries, and an optional Gradio web interface.

The main entry point in this fork is now:

- `af2_pmpnnV7.sh` — primary shell pipeline
- `web_run_pipeline_full.py` — optional web interface that wraps the shell pipeline

---

## Overview

This pipeline is designed to:

1. take a single input structure or FASTA file
2. identify VH/VL regions and CDRs with ANARCI
3. build an scFv sequence and define non-designable positions
4. optionally run ColabFold to generate a backbone
5. run ProteinMPNN to generate sequence designs
6. score designs with:
   - `scfvtools`
   - `struct-evo`
   - SWI / probability of solubility
7. merge scores, rank designs, and select Top / Bottom / Random sets
8. optionally run AF2 on WT + selected designs
9. generate summary HTML and output tables

---

## Major updates in this fork

Compared with the original repository, this fork emphasizes:

- `af2_pmpnnV7.sh` as the main production pipeline
- `web_run_pipeline_full.py` as an optional Gradio front end
- `scripts/make_scfv_vernierV3.py` for VH/VL extraction, CDR handling, linear sequence↔structure mapping, and Vernier shell detection
- pre-AF2 score merging and ranking
- Top / Bottom / Random design selection
- summary HTML augmentation with ranked tables and colorized sequence blocks

Useful helper scripts added in this fork include:

- `scripts/make_scfv_vernierV3.py`
- `scripts/merge_preAF2_scores.py`
- `scripts/merge_swi_evo_json.py`
- `scripts/select_top_bottom_random.py`
- `scripts/swi_calculator.py`
- `scripts/append_scores_to_summary.py`
- `scripts/append_designs_to_summary.py`

---

## Installation / required software

This repository does **not** set up all environments automatically. You will need to install and manage the required tools yourself.

### Required tools / environments

At minimum, this pipeline expects working installations of:

- **ANARCI** – https://github.com/oxpig/ANARCI  
- **ColabFold / localcolabfold** – https://github.com/sokrypton/ColabFold  
- **ProteinMPNN** – https://github.com/dauparas/ProteinMPNN  
- **scfvtools** – https://github.com/Colorado-State-University-Stasevich-Lab/scfvtools  
- **structural-evolution (struct-evo)** – https://github.com/varun-shanker/structural-evolution  

### Important note about environments

The shell script currently **switches among multiple conda environments during execution**. In its current form, it expects environment activation to work and expects certain tools to exist in those environments.

Each stage of the pipeline depends on a different environment:

| Stage                | Environment     |
|---------------------|----------------|
| Pipeline start      | pmpnn          |
| AF2 folding         | colabfold      |
| ProteinMPNN         | pmpnn          |
| scfvtools scoring   | scfvtools_env  |
| struct-evo scoring  | struct-evo     |

So before using this repository on a new system, review the environment names and paths near the top and middle of `af2_pmpnnV7.sh`.

### Important path assumptions

This fork currently contains machine-specific path assumptions, for example:

- `/home/tagteam/miniforge3/etc/profile.d/conda.sh`
- `/home/tagteam/code/ProteinMPNN`
- `~/Projects/structural-evolution/bin/score_log_likelihoods.py`

You will likely need to edit these for your own machine, or override them where supported.

---

## Main way to run: `af2_pmpnnV7.sh`

This is the main workflow in the repository.

### Basic shell-script workflow

The shell script is written around two directories:

- `input_dir`
- `output_dir`

By default, the script looks for **one input file** inside `input_dir`, processes that file, and writes results into `output_dir`.

### Step 1: prepare an input directory

Create an input directory and place **exactly one input file** inside it.

Supported input file types include:

- `.pdb`
- `.cif`
- `.mmcif`
- `.fa`
- `.fasta`

Examples:

```bash
mkdir -p input_dir
cp my_structure.pdb input_dir/
```

or

```bash
mkdir -p input_dir
cp my_sequence.fasta input_dir/
```

### Step 2: decide whether to use defaults or override variables

At the top of `af2_pmpnnV7.sh`, the pipeline defines user variables like this:

```bash
input_dir="${input_dir:-input_dir}"
output_dir="${output_dir:-output_dir}"
determine_CDRs="${determine_CDRs:-martin}"
...
```

This means:

- if you set an environment variable before launching the script, that value is used
- otherwise the script falls back to the default shown in the file

So you can either:

#### Option A: edit defaults in the script itself

For example, change:

```bash
input_dir="${input_dir:-input_dir}"
output_dir="${output_dir:-output_dir}"
```

to whatever you want.

#### Option B: override on the command line

Example:

```bash
input_dir=my_inputs output_dir=my_outputs determine_CDRs=martin ss_near_CDRs=3 linker_seq=GGGGSGGGGSGGGGS seqs_per_run=100 epitope_chain=B skip_folding=false pmpnn_seed=37 Nt=10 Nb=10 Nr=10 sort_column=scfvtools_blosum_diff_score run_AF2=false bash af2_pmpnnV7.sh
```

This is usually the cleanest way to run the pipeline without editing the script.

---

## User-configurable shell variables

The main shell variables currently include:

- `input_dir` — directory containing the single input file
- `output_dir` — directory where outputs are written
- `determine_CDRs` — CDR numbering scheme, e.g. `martin`, `chothia`, `kabat`, `IMGT`
- `ss_near_CDRs` — protection distance near CDRs in Å
- `linker_seq` — scFv linker sequence
- `seqs_per_run` — number of ProteinMPNN sequences to generate
- `epitope_chain` — epitope chain ID(s) for structure input, e.g. `B` or `B,C,F`
- `skip_folding` — whether to skip initial AF2/ColabFold folding when possible
- `PMPNN` — path to the ProteinMPNN installation
- `pmpnn_seed` — ProteinMPNN seed
- `Nt` — number of top designs to keep
- `Nb` — number of bottom designs to keep
- `Nr` — number of random designs to keep
- `sort_column` — score column used for ranking
- `run_AF2` — whether to run AF2 on WT + selected designs

---

## How input handling works

### Structure input: PDB / CIF / mmCIF

If the input is a structure:

- the script reads the single file in `input_dir`
- CIF / mmCIF can be converted to PDB for compatibility
- the epitope chain(s) can be specified with `epitope_chain`
- the scFv is expected to correspond to chain `A` for downstream design

### FASTA input

If the input is a FASTA file:

- the pipeline detects FASTA automatically
- if the sequence contains `SCFV:EPITOPE`, the script splits this into scFv and epitope parts
- the scFv-only portion is passed into the scFv/Vernier logic
- FASTA inputs are forced through AF2/ColabFold, because ProteinMPNN requires a structure backbone
- for FASTA input, `epitope_chain` is effectively set to `None`

---

## Actual shell-script flow

The current shell pipeline proceeds roughly as follows:

### Step 0 — load one input file

The script finds the first file in `input_dir`. If no file is present, it exits with an error.

### Step 1 — determine CDRs + build scFv FASTA

`scripts/make_scfv_vernierV3.py` is used to:

- identify VH and VL
- define CDR regions
- combine VH + linker + VL into an scFv
- determine fixed / designable positions
- generate mapping files and an HTML summary scaffold

### Step 2 — initial ColabFold / AF2 backbone generation

The script switches to the `colabfold` environment and either:

- skips folding and reuses the input structure, or
- runs `colabfold_batch`

If the input was FASTA, folding is forced.

### Step 3 — ProteinMPNN design

The script switches back to the `pmpnn` environment and:

- parses the chosen structure
- assigns chain A for design
- creates fixed-position JSONL files
- runs `protein_mpnn_run.py`

### Step 4 — build combined FASTA of WT + designs

The script assembles a combined FASTA that includes WT plus generated designs.

### Step 5 — scoring and pre-AF2 ranking

This stage includes:

- `scfvtools` scoring
- switching to the `struct-evo` environment for struct-evo scoring
- SWI calculation
- merging scores with `scripts/merge_preAF2_scores.py`
- selecting Top / Bottom / Random designs with `scripts/select_top_bottom_random.py`

### Step 5E — FASTA subsets + summary augmentation

The pipeline builds FASTA files for:

- top designs
- bottom designs
- random designs

and appends those to the summary HTML.

### Step 6 — optional AF2 on WT + selected designs

If `run_AF2=true`, WT plus selected Top / Bottom / Random designs are folded with AF2 / ColabFold.

### Step 7 — merge post-AF2 metrics

AF2 JSON metrics are merged with prior scores using `scripts/merge_swi_evo_json.py`.

### Step 8 — summary output

HTML summary content is extended with:

- ranked score tables
- colorized design sequences
- mutation highlighting

---

## Key outputs

The exact file set depends on settings, but important outputs typically include:

### In `output_dir/`

- `scfv_output.fasta`
- `combined_multimer.fa`
- `merged_preAF2_scores.csv`
- AF2 / ColabFold result folders
- ProteinMPNN output folders

### In `output_dir/Summary/`

- `*_summary.html`
- `selected_top.csv`
- `selected_bottom.csv`
- `selected_random.csv`
- `selected_for_af2.txt`
- `top_designs.fa`
- `bottom_designs.fa`
- `random_designs.fa`
- mapping / design-position text files

---

## Optional web interface: `web_run_pipeline_full.py`

The web interface is optional and simply wraps the same shell pipeline. This would typically run in a separate environment:

Example setup:

conda create -n scfv_web python=3.10
conda activate scfv_web
pip install gradio fastapi uvicorn

### What it does

The Gradio app:

- creates a temporary working directory
- creates temporary `input_dir` and `output_dir`
- copies the `scripts/` directory into that working directory
- passes user-selected options into the shell script through environment variables
- streams the live shell log
- packages the final outputs into `results.zip`
- renders a summary HTML preview if present

### Run it

```bash
python web_run_pipeline_full.py
```

Then open the local Gradio URL, usually:

```text
http://localhost:7860
```

### Web interface usage summary

- upload a single PDB/CIF/mmCIF/FASTA file
- choose CDR scheme, linker, protection distance, epitope chain(s), GPU, seed, and design counts
- optionally enable AF2 on WT + designs
- click **Run pipeline**
- inspect the live log
- download `results.zip`
- preview the summary HTML in the browser

---

## Suggested usage pattern for new users

A practical first test is:

1. install / verify all required environments
2. edit `PMPNN` and any machine-specific paths in `af2_pmpnnV7.sh`
3. create `input_dir`
4. place one test PDB or FASTA into `input_dir`
5. run:

```bash
bash af2_pmpnnV7.sh
```

or with explicit overrides:

```bash
input_dir=input_dir output_dir=output_dir seqs_per_run=20 run_AF2=false bash af2_pmpnnV7.sh
```

6. inspect the HTML and CSV outputs in `output_dir/Summary`

---

## Acknowledgments

This workflow builds on and depends on several important external tools, including:

- ANARCI
- ProteinMPNN
- ColabFold / AlphaFold2
- scfvtools
- structural-evolution

Please also acknowledge the original upstream repository that this fork extends.
