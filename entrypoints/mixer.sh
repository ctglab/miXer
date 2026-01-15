#!/bin/bash 
set -e

# these are suffixes for subfolders
OUTDIR_EXCA="excavator2_output"
OUTDIR_SVM="svm_processed_output"
OUTDIR_HMM="mixer_windows"
OUTDIR_VCF="mixer_vcfs"

# This is intended to be the entrypoint for the miXer container.
ascii_logo="CgogICAgICAgICAgICAgICAgICDilojilojiloggIOKWiOKWiOKWiOKWiOKWiCDilojilojilojilojiloggICAgICAgICAgICAgICAgICAgCiAgICAgICAgICAgICAgICAg4paR4paR4paRICDilpHilpHilojilojilogg4paR4paR4paI4paI4paIICAgICAgICAgICAgICAgICAgICAKIOKWiOKWiOKWiOKWiOKWiOKWiOKWiOKWiOKWiOKWiOKWiOKWiOKWiCAgIOKWiOKWiOKWiOKWiCAg4paR4paR4paI4paI4paIIOKWiOKWiOKWiCAgICDilojilojilojilojilojiloggIOKWiOKWiOKWiOKWiOKWiOKWiOKWiOKWiCAK4paR4paR4paI4paI4paI4paR4paR4paI4paI4paI4paR4paR4paI4paI4paIIOKWkeKWkeKWiOKWiOKWiCAgIOKWkeKWkeKWiOKWiOKWiOKWiOKWiCAgICDilojilojilojilpHilpHilojilojilojilpHilpHilojilojilojilpHilpHilojilojilogKIOKWkeKWiOKWiOKWiCDilpHilojilojilogg4paR4paI4paI4paIICDilpHilojilojiloggICAg4paI4paI4paI4paR4paI4paI4paIICDilpHilojilojilojilojilojilojiloggIOKWkeKWiOKWiOKWiCDilpHilpHilpEgCiDilpHilojilojilogg4paR4paI4paI4paIIOKWkeKWiOKWiOKWiCAg4paR4paI4paI4paIICAg4paI4paI4paIIOKWkeKWkeKWiOKWiOKWiCDilpHilojilojilojilpHilpHilpEgICDilpHilojilojiloggICAgIAog4paI4paI4paI4paI4paI4paR4paI4paI4paIIOKWiOKWiOKWiOKWiOKWiCDilojilojilojilojilogg4paI4paI4paI4paI4paIIOKWiOKWiOKWiOKWiOKWiOKWkeKWkeKWiOKWiOKWiOKWiOKWiOKWiCAg4paI4paI4paI4paI4paIICAgIArilpHilpHilpHilpHilpEg4paR4paR4paRIOKWkeKWkeKWkeKWkeKWkSDilpHilpHilpHilpHilpEg4paR4paR4paR4paR4paRIOKWkeKWkeKWkeKWkeKWkSAg4paR4paR4paR4paR4paR4paRICDilpHilpHilpHilpHilpEgICAgIAoKCgo="

eval "$(micromamba shell hook --shell bash)"
echo "$ascii_logo" | base64 --decode

# check args 
if [[ "$1" = "preprocessing" ]]; then
    shift 
    JSON_FILE=""
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -j|--json)
                JSON_FILE="$2"; shift 2;;
            *)
                echo "Unknown option: $1"; 
                echo "Usage: mixer.sh preprocessing <config_file>"
                exit 1;;
        esac
    done
    # check required
    if [[ -z "$JSON_FILE" ]]; then
        echo "Missing required argument."
        echo "Usage: mixer.sh preprocessing -j <json_file>"
        exit 1
    fi
    # activate the environment
    micromamba activate excavator2 
    echo "Running preprocessing:"
    python3 -u /app/preprocessingMixer/run_excavator2singularity.py -j "$JSON_FILE"
    
elif [[ "$1" = "inference" ]]; then
        # consume the sub‐command
        shift
        # initialize variables
        JSON_FILE=""
        SAMPLES_FILE=""
        XLR_FILE=""
        SEGDUP_FILE=""
        REFERENCE=""
        BW_ITERATIONS=""
        BW_DELTA=""

        # parse flags
        while [[ $# -gt 0 ]]; do
            case "$1" in
                -j|--json)
                    JSON_FILE="$2"; shift 2;;
                -s|--samples)
                    SAMPLES_FILE="$2"; shift 2;;
                -x|--xlr)
                    XLR_FILE="$2"; shift 2;;
                -sd|--segdup)
                    SEGDUP_FILE="$2"; shift 2;;
                -bw|--baum_welch_iterations)
                    BW_ITERATIONS="$2"; shift 2;;
                -delta|--baum_welch_delta)
                    BW_DELTA="$2"; shift 2;;
                *)
                    echo "Unknown option: $1"; 
                    echo "Usage: mixer.sh inference -j <json> -s <samples> [-x <xlr>] [-sd <segdup>] [-bw <baum_welch_iterations>] [-delta <baum_welch_delta>]"
                    exit 1;;
            esac
        done

        # check required
        if [[ -z "$JSON_FILE" ]]; then
            echo "Missing required argument."
            echo "Usage: mixer.sh inference -j <json> [-s <samples>] [-x <xlr>] [-sd <segdup>] [-bw <baum_welch_iterations>] [-delta <baum_welch_delta>]"
            exit 1
        fi

        micromamba activate mixerPre
        echo "Running inference with:"
        echo "  JSON:    $JSON_FILE"
        [[ -n "$SAMPLES_FILE" ]] && echo "  samples: $SAMPLES_FILE"
        [[ -n "$XLR_FILE"     ]] && echo "  xlr:     $XLR_FILE"
        [[ -n "$SEGDUP_FILE"  ]] && echo "  segdup:  $SEGDUP_FILE"
        [[ -n "$BW_ITERATIONS" ]] && echo "  Baum-Welch iterations: $BW_ITERATIONS"
        [[ -n "$BW_DELTA"     ]] && echo "  Baum-Welch delta: $BW_DELTA"
        [[ -n "$REFERENCE"   ]] && echo "  Reference version: $REFERENCE"

        # Check for enable_intermediate_results_output in JSON (default to false if not present)
        # We do this early to warn the user
        ENABLE_INTERMEDIATE=$(python3 -c "import json, sys; config=json.load(open('$JSON_FILE')); print(str(config.get('enable_intermediate_results_output', False)).lower())")
        
        if [[ "$ENABLE_INTERMEDIATE" != "true" ]]; then
             echo "WARNING: Intermediate files (datasets_testing, *_SVC) will NOT be saved."
             echo "         To save them, set \"enable_intermediate_results_output\": true in your JSON config."
        else
             echo "INFO: Intermediate files will be preserved."
        fi

        python3 -u /app/preprocessingMixer/generate_miXer_datasets.py \
            -j "$JSON_FILE" \
            ${SAMPLES_FILE:+-s "$SAMPLES_FILE"} \
            ${XLR_FILE:+-x "$XLR_FILE"} \
            ${SEGDUP_FILE:+-sd "$SEGDUP_FILE"}

        micromamba activate miXer_ml
        echo "Running SVM inference"
        python3 -u /app/processing/scripts/miXer_inference.py \
            -j "$JSON_FILE" 

        micromamba activate HMMR
        echo "Running HMM filtering"
        Rscript /app/processing/scripts/miXe.R \
            -j "$JSON_FILE" \
            ${BW_ITERATIONS:+-b "$BW_ITERATIONS"} \
            ${BW_DELTA:+-d "$BW_DELTA"} 
        
        micromamba activate miXer_ml
        echo "Writing VCFs"
        python3 -u /app/processing/scripts/Vcf_maker.py \
            -j "$JSON_FILE" \
            ${REFERENCE:+-ref "$REFERENCE"} \

        if [[ "$ENABLE_INTERMEDIATE" != "true" ]]; then
             echo "Cleaning up intermediate results..."
             # Get main outdir from JSON
             MAIN_OUTDIR=$(python3 -c "import json, sys; config=json.load(open('$JSON_FILE')); print(config['main_outdir_host'])")
             EXP_ID=$(python3 -c "import json, sys; config=json.load(open('$JSON_FILE')); print(config['exp_id'])")
             TARGET_DIR="${MAIN_OUTDIR}/${EXP_ID}"
             
             # remove datasets_testing_* directories
             find "$TARGET_DIR" -maxdepth 1 -type d -name "datasets_testing_*" -exec rm -rf {} +
             
             # remove *_SVC directories
             find "$TARGET_DIR" -maxdepth 1 -type d -name "*_SVC" -exec rm -rf {} +
             echo "Intermediate results cleaned up."
        else
             echo "Preserving intermediate results."
        fi

else
    # if no valid sub-command is provided, print usage
    echo "Usage: mixer.sh <sub-command> [options]"
    echo "Sub-commands:"
    echo "  preprocessing -j <json_file>   Run preprocessing with the specified JSON file"
    echo "  inference -j <json> [-s <samples>] [-x <xlr>] [-sd <segdup>]"
    echo "Options:"
    echo "    -j, --json <file>          JSON configuration file for preprocessing"
    echo "    -s, --samples <file>       Samples file for inference"
    echo "    -x, --xlr <file>           XLR file for inference (optional)"
    echo "    -sd, --segdup <file>       SegDup file for inference (optional)"
    echo "    -bw, --baum_welch_iterations <iterations>  Baum-Welch iterations for HMM (optional)"
    echo "    -delta, --baum_welch_delta <delta>  Baum-Welch delta for HMM (optional)"
    echo "    -r, --reference <string>   Reference genome version (optional)"
    exit 1
fi 
