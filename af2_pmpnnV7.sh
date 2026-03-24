#!/bin/bash
# @tjstasevich based on script from @jderoo
# V5-clean: AF2 multimer + scFv MPNN design (minimal, no redundant logic)

###############################################################
### USER CONFIGURATION (OVERRIDABLE VIA ENV VARS)
###############################################################
# If the variable is set in the environment, use that.
# Otherwise, fall back to the default in the script.

input_dir="${input_dir:-input_dir}"                   # name of input directory
output_dir="${output_dir:-output_dir}"                # name of output directory
determine_CDRs="${determine_CDRs:-martin}"            # CDR numbering scheme to use
ss_near_CDRs="${ss_near_CDRs:-3}"                     # don't change AAs this close to CDRs (angstroms)
linker_seq="${linker_seq:-GGGGSGGGGSGGGGS}"           # linker sequence to use
seqs_per_run="${seqs_per_run:-20}"                     # how many sequences to design with MPNN
epitope_chain="${epitope_chain:-B}"                   # which chain is the epitope in input?
skip_folding="${skip_folding:-true}"                 # if input is scFv of form MA+VH+linker+VL+KA
PMPNN="${PMPNN:-/home/tagteam/code/ProteinMPNN}"
pmpnn_seed="${pmpnn_seed:-37}"                       # default ProteinMPNN seed

# scoring / selection controls
Nt="${Nt:-10}"                       # number of top designs
Nb="${Nb:-10}"                       # number of bottom designs
Nr="${Nr:-10}"                       # number of random designs
sort_column="${sort_column:-scfvtools_blosum_diff_score}"   # which column to sort by
run_AF2="${run_AF2:-false}"           # set to 'false' to skip AF2 folding
###############################################################
### Helpers
###############################################################
# Normalize epitope_chain to an array of chain IDs.
# Accepts: "None", "B", "B,C,F", "B C F"
parse_epitope_chains () {
    local raw="$1"
    raw="${raw//;/,}"
    raw="${raw// /,}"
    # collapse multiple commas
    raw=$(echo "$raw" | tr -s ',')
    IFS=',' read -r -a _chains <<< "$raw"
    epitope_chains=()
    for c in "${_chains[@]}"; do
        c="$(echo "$c" | tr -d '[:space:]')"
        if [ -z "$c" ]; then
            continue
        fi
        if [ "$c" = "None" ] || [ "$c" = "none" ]; then
            continue
        fi
        epitope_chains+=("$c")
    done
}

join_by_colon () {
    local IFS=":"
    echo "$*"
}

###############################################################

# Activate env
if [ "$CONDA_DEFAULT_ENV" != "pmpnn" ]; then
    source /home/tagteam/miniforge3/etc/profile.d/conda.sh
    conda activate pmpnn
fi

mkdir -p "$output_dir"
summary_dir="$output_dir/Summary"
mkdir -p "$summary_dir"

###############################################################
### STEP 0: Load single file
###############################################################
input_file=$(find "$input_dir" -maxdepth 1 -type f | head -n 1)
# --- MINIMAL SAFETY CHECK #1 ---
if [ -z "$input_file" ]; then
    echo "ERROR: No input file found in $input_dir"
    exit 1
fi
# -------------------------------
base_name=$(basename "$input_file")
base_noext="${base_name%.*}"

echo "Input file = $input_file"

###############################################################
### STEP 0A: Detect input type + (optional) CIF → PDB conversion
###############################################################
ext="${input_file##*.}"
ext="${ext,,}"   # lowercase
is_fasta=false
epitope_seq_from_fasta=""   # used only for FASTA scFv:epitope mode

if [[ "$ext" == "fa" || "$ext" == "fasta" ]]; then
    is_fasta=true
    echo "Detected FASTA input."
elif [[ "$ext" == "cif" || "$ext" == "mmcif" ]]; then
    echo "Detected CIF input. Converting to PDB for compatibility..."

    python - << EOF
import gemmi, pathlib
cif_path = pathlib.Path("$input_file")
pdb_path = cif_path.with_suffix(".pdb")
doc = gemmi.cif.read_file(str(cif_path))
structure = gemmi.make_structure_from_block(doc.sole_block())
structure.write_pdb(str(pdb_path))
EOF

    # Update input_file to point to the new PDB
    input_file="${input_file%.*}.pdb"
    base_name=$(basename "$input_file")
    base_noext="${base_name%.*}"
    ext="pdb"

    echo "Converted PDB written to: $input_file"
fi

###############################################################
### STEP 0B: FASTA scFv:epitope splitter (MUST run before ANARCI)
###############################################################
# If FASTA contains:
#   >name
#   SCFV_SEQUENCE:EPITOPE_SEQUENCE
# then ANARCI must see ONLY the scFv part (no ':').
# We also capture the epitope part for later ColabFold multimer folding.
if [ "$is_fasta" = true ]; then
    echo "Checking FASTA for multimer scFv:epitope syntax (pre-ANARCI)..."

    # Concatenate all non-header lines (handles wrapped FASTA)
    full_seq=$(awk 'BEGIN{seq=""} !/^>/ {gsub(/[ \t\r\n]/,""); seq=seq$0} END{print seq}' "$input_file")

    if [ -z "$full_seq" ]; then
        echo "ERROR: FASTA contains no sequence."
        exit 1
    fi

    if echo "$full_seq" | grep -q ":"; then
        scfv_part="${full_seq%%:*}"
        epitope_part="${full_seq#*:}"

        if [ -z "$scfv_part" ] || [ -z "$epitope_part" ]; then
            echo "ERROR: Found ':' in FASTA but could not split into scFv and epitope."
            exit 1
        fi

        echo "Detected scFv:epitope in FASTA."
        echo "  scFv length:    ${#scfv_part}"
        echo "  epitope length: ${#epitope_part}"

        # Write scFv-only FASTA for ANARCI / make_scfv_vernierV3.py
        scfv_only_fasta="$output_dir/scfv_only_input.fasta"
        mkdir -p "$output_dir"
        echo ">scFv" > "$scfv_only_fasta"
        echo "$scfv_part" >> "$scfv_only_fasta"

        # From here on, use scFv-only FASTA for CDR/Vernier extraction
        input_file="$scfv_only_fasta"
        base_name=$(basename "$input_file")
        base_noext="${base_name%.*}"

        # Save epitope for later multimer folding (chain B)
        epitope_seq_from_fasta="$epitope_part"

        # FASTA has no chain IDs to extract; disable structure-based epitope extraction
        epitope_chain="None"
    else
        # FASTA is scFv-only
        epitope_seq_from_fasta=""
        epitope_chain="None"
        echo "No ':' found; treating FASTA as scFv-only."
    fi
fi

###############################################################
### STEP 1: Determine CDRs + scFv FASTA
###############################################################
# Use new scFv/Vernier script
cmd="python scripts/make_scfv_vernierV3.py \
    --scheme $determine_CDRs \
    --combine \
    --linker $linker_seq \
    --save-anarci \
    --dist $ss_near_CDRs \
    --simple-output design \
    $input_file"

# Execute and capture the printed list of positions (fail fast if ANARCI fails)
design_positions="$($cmd)"
rc=$?
if [ $rc -ne 0 ]; then
    echo "ERROR: make_scfv_vernierV3.py failed (exit code $rc)."
    exit 1
fi

# Guard against missing output (prevents cascading mv/awk failures)
if [ ! -f "scfv_output.fasta" ]; then
    echo "ERROR: scfv_output.fasta was not produced by make_scfv_vernierV3.py."
    exit 1
fi

echo "Design positions: $design_positions"

# Write to summary for record-keeping
echo "$design_positions" > "$summary_dir/${base_noext}_design_positions.txt"

# Move summary files (if they exist)
for f in \
    "${base_noext}_summary.html" \
    "${base_noext}_fixed_positions.txt" \
    "${base_noext}_anarci.txt" \
    VH_mapping.tsv \
    VL_mapping.tsv \
    pdb_linear_map.tsv; do
    [ -f "$f" ] && mv "$f" "$summary_dir/"
done

# Move output FASTA produced by make_scfv_vernierV3.py
mv scfv_output.fasta "$output_dir/"
combined_fasta="$output_dir/scfv_output.fasta"
scfv_seq=$(awk '/^>scFv/{getline; print; exit}' "$combined_fasta")

###############################################################
### STEP 2: Run LocalColabFold on scFv:epitope model
###############################################################

echo "DEBUG: Using python $(which python)"

cf_out="$output_dir/colabfold_results"
mkdir -p "$cf_out"
mpnn_fasta="$output_dir/${base_noext}_mpnn_input.fasta"

# Always extract epitope sequence(s) if epitope_chain is defined (PDB/CIF mode)
parse_epitope_chains "$epitope_chain"

epitope_seq=""

# If FASTA provided epitope via scFv:epitope syntax, use it here.
if [ "$is_fasta" = true ] && [ -n "${epitope_seq_from_fasta:-}" ]; then
    epitope_seq="$epitope_seq_from_fasta"
fi

if [ "${#epitope_chains[@]}" -gt 0 ]; then
    if [ "$is_fasta" = true ]; then
        echo "ERROR: epitope_chain was provided ($epitope_chain) but input is FASTA."
        echo "       For FASTA scFv:epitope input, set epitope_chain=None."
        exit 1
    fi

    epitope_seqs=()
    for c in "${epitope_chains[@]}"; do
        seq=$(python scripts/extract_chain_seq.py "$input_file" "$c")
        if [ -z "$seq" ]; then
            echo "ERROR: Could not extract sequence for epitope chain '$c' from $input_file"
            exit 1
        fi
        epitope_seqs+=("$seq")
    done

    # Colon-separated sequence string for LocalColabFold multimer mode
    epitope_seq=$(join_by_colon "${epitope_seqs[@]}")
fi

# Activate ColabFold environment
source /home/tagteam/miniforge3/etc/profile.d/conda.sh
conda deactivate
conda activate colabfold

# Ensure JAX/colabfold can find CUDA + cuDNN libs in this conda env
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

# If user wants to skip folding (e.g., input is already a good scFv model)
# NOTE: If the input is FASTA, we MUST fold to obtain a PDB backbone for MPNN.
if [ "$is_fasta" = true ]; then
    echo "FASTA input → forcing LocalColabFold to generate a PDB backbone (cannot skip folding)."
    skip_folding=false
fi

if [ "$skip_folding" = true ]; then
    echo "Skipping LocalColabFold — using input structure directly."
    rank_file="$input_file"
else
    # Otherwise fold
    if [ -z "$epitope_seq" ]; then
        # Monomer run
        echo -e ">scFv\n$scfv_seq" > "$mpnn_fasta"

        colabfold_batch \
            --templates \
            --num-recycle 6 \
            "$mpnn_fasta" "$cf_out"
    else
        # Multimer: scFv + epitope sequence(s) (epitope becomes chain B)
        echo -e ">scFv_epitope\n${scfv_seq}:${epitope_seq}" > "$mpnn_fasta"

        colabfold_batch \
            --templates \
            --num-recycle 6 \
            --model-type alphafold2_multimer_v2 \
            "$mpnn_fasta" "$cf_out"
    fi

    # Identify the Alphafold result to use
    rank_file=$(find "$cf_out" -maxdepth 1 -type f -name "*rank_001*.pdb" | head -n 1)
fi

echo "Using structure for MPNN: $rank_file"

###############################################################
# STEP 3: Run ProteinMPNN 
###############################################################

# Activate MPNN environment
echo "Activating MPNN environment"
source /home/tagteam/miniforge3/etc/profile.d/conda.sh
conda deactivate
conda activate pmpnn

#export CUDA_VISIBLE_DEVICES=""
#echo "Forcing CPU execution for ProteinMPNN (GB10 workaround)"

dir=$output_dir"/MPNN_input"
mkdir -p "$dir"
cp "$rank_file" "$dir/"

# Old-style JSONL filenames (single files)
parsed_jsonl=$dir"/parsed_pdbs.jsonl"
assigned_jsonl=$dir"/assigned_pdbs.jsonl"
fixed_jsonl=$dir"/fixed_pdbs.jsonl"

chains_to_design="A"
echo "Chains to design = $chains_to_design"
echo "Running ProteinMPNN..."
#printf 'RAW DESIGN POSITIONS = [%q]\n' "$design_positions"

# 1. Parse rank_001 → parsed_pdbs.jsonl
# --- MINIMAL SAFETY CHECK #2 ---
if [ -z "$rank_file" ]; then
    echo "ERROR: rank_file is empty — folding step failed"
    exit 1
fi
# --------------------------------

# Delete ONLY rank_* files EXCEPT rank_001.pdb
find "$dir" -maxdepth 1 -type f -name "rank_*.pdb" ! -name "rank_001.pdb" -delete

# 1. Parse chains → parsed_jsonl.jsonl
python "$PMPNN/helper_scripts/parse_multiple_chains.py" \
    --input_path "$dir" \
    --output_path "$parsed_jsonl"


# 2. Assign chains → assigned_pdbs.jsonl
python "$PMPNN/helper_scripts/assign_fixed_chains.py" \
    --input_path "$parsed_jsonl" \
    --output_path "$assigned_jsonl" \
    --chain_list "$chains_to_design"

# 3. Create fixed positions → fixed_pdbs.jsonl
python "$PMPNN/helper_scripts/make_fixed_positions_dict.py" \
    --input_path "$parsed_jsonl" \
    --output_path "$fixed_jsonl" \
    --chain_list "$chains_to_design" \
    --position_list "$design_positions" \
    --specify_non_fixed

# 4. Run ProteinMPNN (CPU mode)
python "$PMPNN/protein_mpnn_run.py" \
    --jsonl_path "$parsed_jsonl" \
    --chain_id_jsonl "$assigned_jsonl" \
    --fixed_positions_jsonl "$fixed_jsonl" \
    --out_folder "$output_dir" \
    --num_seq_per_target "$seqs_per_run" \
    --sampling_temp "0.1" \
    --seed "$pmpnn_seed" \
    --batch_size 1 \
    --backbone_noise 0.02 \
    --use_soluble_model

mpnn_rc=$?
if [ $mpnn_rc -ne 0 ]; then
    echo "ERROR: ProteinMPNN failed"
    exit 1
fi

# Verify sequences were actually produced
if [ ! -d "$output_dir/seqs" ] || ! ls "$output_dir/seqs"/*.fa >/dev/null 2>&1; then
    echo "ERROR: ProteinMPNN produced no FASTA outputs"
    exit 1
fi

echo "ProteinMPNN completed successfully."


###############################################################
# STEP 4: Build one combined FASTA (WT + designs) for folding
###############################################################

python scripts/make_combined_multimer_fasta.py \
    "$output_dir/seqs" \
    "$epitope_seq" \
    "$output_dir/combined_multimer.fa"


summary_html="$summary_dir/${base_noext}_summary.html"
combined_fa="$output_dir/combined_multimer.fa"

python scripts/append_designs_to_summary.py \
    "$summary_html" \
    "$combined_fa"


###############################################################
# STEP 5: Score designs (SWI, struct-evo, scfvtools) – pre-AF2
###############################################################

# Paths for outputs
clean_fa="$output_dir/combined_multimer_chainA_only.fa"
evo_csv="$output_dir/scores_chainA.csv"
swi_csv="$output_dir/swi_scores.csv"
scfvtools_csv="$output_dir/scfvtools_scores.csv"
merged_preaf2_csv="$output_dir/merged_preAF2_scores.csv"

echo "Cleaning FASTA (struct-evo, SWI, scfvtools only need chain A)..."
awk -v OFS="" '
    /^>/ { print; next }
    {
        # Split on colon and keep only the first sequence (chain A)
        split($0, arr, ":");
        print arr[1]
    }
' "$combined_fa" > "$clean_fa"
echo "Clean FASTA written to: $clean_fa"

########################################
# 5A. SCFVTOOLS scoring
########################################
echo "Running scfvtools scoring..."
source /home/tagteam/miniforge3/etc/profile.d/conda.sh
conda deactivate
conda activate scfvtools_env

python ~/Projects/scfvtools/scripts/run_scfvtools_scoring.py \
    --fasta "$clean_fa" \
    --scheme martin \
    --reference_csv "~/Projects/scfvtools/example_data/diff_H.csv" \
    --reference_name "diff_H" \
    --out_csv "$scfvtools_csv" \
    --summary_html "$summary_html"

echo "SCFVTOOLS scores written to: $scfvtools_csv"


########################################
# 5B. struct-evo scoring
########################################
echo "Switching to struct-evo environment..."
source /home/tagteam/miniforge3/etc/profile.d/conda.sh
conda deactivate
conda activate struct-evo

# Use the WT backbone from the earlier folding (rank_file from STEP 2),
# not from the designs AF2 run (which hasn't happened yet)
if [ "${#epitope_chains[@]}" -eq 0 ]; then
    echo "No epitope chain → scoring scFv only (single-chain mode)"
    struct_evo_args="--chain A"
else
    echo "Epitope chain detected → scoring scFv + epitope (multichain mode)"
    struct_evo_args="--multichain-backbone"
fi

echo "Running struct-evo scoring..."
python ~/Projects/structural-evolution/bin/score_log_likelihoods.py \
    "$rank_file" \
    $struct_evo_args \
    --seqpath "$clean_fa" \
    --outpath "$evo_csv"

echo "struct-evo scores written to: $evo_csv"

########################################
# 5C. SWI scoring
########################################
echo "Running SWI calculator..."
python ./scripts/swi_calculator.py \
    "$clean_fa" \
    -o "$swi_csv"

echo "SWI scores written to: $swi_csv"


###############################################################
# STEP 5D: Merge scores + choose Top / Bottom / Random
###############################################################

merged_preaf2_csv="$output_dir/merged_preAF2_scores.csv"


echo "Merging pre-AF2 scores..."

sort_column="${sort_column:-scfvtools_blosum_diff_score}"

python scripts/merge_preAF2_scores.py \
    --swi "$swi_csv" \
    --evo "$evo_csv" \
    --scfv "$scfvtools_csv" \
    --sort_by "$sort_column" \
    --out "$merged_preaf2_csv"


echo "Selecting Top ($Nt), Bottom ($Nb), Random ($Nr) using sort column '$sort_column'..."

python scripts/select_top_bottom_random.py \
    --input "$merged_preaf2_csv" \
    --sortcol "$sort_column" \
    --Nt "$Nt" \
    --Nb "$Nb" \
    --Nr "$Nr" \
    --outdir "$summary_dir"

###############################################################
##########     STEP 5E: MAKE TOP/BOTTOM/RANDOM FASTAs    #####
##########     AND APPEND THEM TO SUMMARY.html           #####
###############################################################

echo "Building FASTA files for Top / Bottom / Random..."

combined_fa="$output_dir/combined_multimer.fa"

top_fa="$summary_dir/top_designs.fa"
bottom_fa="$summary_dir/bottom_designs.fa"
random_fa="$summary_dir/random_designs.fa"

extract_chainA () {
    local acc_file="$1"
    local out_fa="$2"

    mkdir -p "$output_dir/Summary"
    : > "$out_fa"

    # Add WT (always NR==2)
    awk '
        BEGIN { RS=">"; FS="\n" }
        NR==2 {
            print ">"$1;
            print $2;
        }' "$combined_fa" >> "$out_fa"

    while IFS=, read -r accession seq swi sol seqid ll ll_tgt score; do

        # Skip header row
        if [[ "$accession" =~ ^Accession ]]; then
            continue
        fi

        # Extract design number from the accession
        design_num=$(echo "$accession" | grep -oE 'design_([0-9]+)' | grep -oE '[0-9]+')
        if [ -z "$design_num" ]; then
            echo "WARNING: Could not extract design number from: $accession"
            continue
        fi

        # Extract the matching FASTA entry (header + full sequence)
        awk -v num="$design_num" '
            BEGIN { RS=">"; FS="\n" }
            $1 ~ ("_design_" num "$") {
                print ">"$1;
                print $2;
            }' "$combined_fa" >> "$out_fa"

    done < "$acc_file"
}



# Build FASTAs
extract_chainA "$summary_dir/selected_top.csv"    "$top_fa"
extract_chainA "$summary_dir/selected_bottom.csv" "$bottom_fa"
extract_chainA "$summary_dir/selected_random.csv" "$random_fa"

echo "TOP FASTA:    $top_fa"
echo "BOTTOM FASTA: $bottom_fa"
echo "RANDOM FASTA: $random_fa"

# Append each block separately to the summary
python scripts/append_designs_to_summary.py "$summary_html" "$top_fa"
python scripts/append_designs_to_summary.py "$summary_html" "$bottom_fa"
python scripts/append_designs_to_summary.py "$summary_html" "$random_fa"

echo "Top/Bottom/Random FASTA blocks appended to summary."


###############################################################
# STEP 6: Optional AF2 folding ONLY on selected sequences
###############################################################
# Activate MPNN environment
echo "Activating colabfold environment"
source /home/tagteam/miniforge3/etc/profile.d/conda.sh
conda deactivate
conda activate colabfold

# Ensure JAX/colabfold can find CUDA + cuDNN libs in this conda env
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

merged_csv="$output_dir/merged_scores.csv"
reordered_csv="$output_dir/merged_scores_reordered.csv"
output_struct_dir="$output_dir/colabfold_results"   

if [ "$run_AF2" = "false" ]; then
    echo "Skipping AF2 folding (run_AF2=false)."
    af2_ran=false
else
    echo "Preparing FASTA with WT + top/bottom/random designs..."

    selected_fa="$output_dir/selected_for_af2.fa"
    : > "$selected_fa"

    # Always include WT entry from combined_multimer.fa
    echo "Including WT from combined FASTA..."
    awk 'BEGIN{RS=">"; FS="\n"} NR==2 {print ">"$0}' "$combined_fa" >> "$selected_fa"

    # Include top/bottom/random design FASTAs (already built above)
    echo "Including top/bottom/random FASTA blocks..."
    cat "$top_fa" "$bottom_fa" "$random_fa" >> "$selected_fa"

    # Deduplicate by header in case a design appears in more than one set
    awk '
        /^>/ {
            if (seen[$0]++) {
                skip=1
            } else {
                skip=0
            }
        }
        !skip
    ' "$selected_fa" > "${selected_fa}.tmp" && mv "${selected_fa}.tmp" "$selected_fa"

    # Count and report which designs will be folded
    total=$(grep -c '^>' "$selected_fa")
    echo ""
    echo "============================================================"
    echo "Folding the following $total sequences (including WT):"
    echo "------------------------------------------------------------"
    i=0
    while read -r line; do
        if [[ $line == ">"* ]]; then
            ((i++))
            name="${line#>}"
            echo "  [$i/$total] $name"
        fi
    done < "$selected_fa"
    echo "------------------------------------------------------------"
    echo ""

    mkdir -p "$output_struct_dir"


    echo "Running LocalColabFold on selected FASTA..."
if [ "${#epitope_chains[@]}" -eq 0 ]; then
    # Monomer folding (scFv only)
    colabfold_batch             --templates             --num-recycle 6             "$selected_fa" "$output_struct_dir"
else
    # Multimer folding (scFv + one or more epitope chains)
    colabfold_batch             --templates             --num-recycle 6             --model-type alphafold2_multimer_v2             "$selected_fa" "$output_struct_dir"
fi

af2_ran=true
    echo "Merging SWI + struct-evo + AF2 JSON metrics + scfvtools..."
    python ./scripts/merge_swi_evo_json.py \
        --swi_file "$swi_csv" \
        --evo_file "$evo_csv" \
        --json_dir "$output_struct_dir" \
        --scfvtools_file "$scfvtools_csv" \
        --top "$summary_dir/selected_top.csv" \
        --bottom "$summary_dir/selected_bottom.csv" \
        --random "$summary_dir/selected_random.csv" \
        --outpath "$merged_csv"


    # echo "Reordering merged scores..."
    # python ./scripts/reorder_all_scores.py \
    #     --input "$merged_csv" \
    #     --output "$reordered_csv"

    echo "Creating HTML table"
    python ./scripts/append_scores_to_summary.py \
        --csv "$merged_csv" \
        --summary "$summary_dir/${base_noext}_summary.html" \
        --top "$summary_dir/selected_top.csv" \
        --bottom "$summary_dir/selected_bottom.csv" \
        --random "$summary_dir/selected_random.csv"
    echo "STEP 6 completed."
fi

if [ "$run_AF2" = "false" ]; then
    python ./scripts/append_scores_to_summary.py \
        --csv "$merged_preaf2_csv" \
        --summary "$summary_dir/${base_noext}_summary.html" \
        --top "$summary_dir/selected_top.csv" \
        --bottom "$summary_dir/selected_bottom.csv" \
        --random "$summary_dir/selected_random.csv"

    echo "Summary updated without AF2 structural metrics."
fi



