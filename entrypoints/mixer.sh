#!/bin/bash 
set -e

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
        OUTDIR=""
        XLR_FILE=""
        SEGDUP_FILE=""

        # parse flags
        while [[ $# -gt 0 ]]; do
            case "$1" in
                -j|--json)
                    JSON_FILE="$2"; shift 2;;
                -s|--samples)
                    SAMPLES_FILE="$2"; shift 2;;
                -o|--outdir)
                    OUTDIR="$2"; shift 2;;
                -x|--xlr)
                    XLR_FILE="$2"; shift 2;;
                -sd|--segdup)
                    SEGDUP_FILE="$2"; shift 2;;
                *)
                    echo "Unknown option: $1"; 
                    echo "Usage: mixer.sh inference -j <json> -s <samples> -o <outdir> [-x <xlr>] [-sd <segdup>]"
                    exit 1;;
            esac
        done

        # check required
        if [[ -z "$JSON_FILE" || -z "$SAMPLES_FILE" || -z "$OUTDIR" ]]; then
            echo "Missing required argument."
            echo "Usage: mixer.sh inference -j <json> -s <samples> -o <outdir> [-x <xlr>] [-sd <segdup>]"
            exit 1
        fi

        micromamba activate mixerPre
        echo "Running inference with:"
        echo "  JSON:    $JSON_FILE"
        echo "  samples: $SAMPLES_FILE"
        echo "  outdir:  $OUTDIR"
        [[ -n "$XLR_FILE"     ]] && echo "  xlr:     $XLR_FILE"
        [[ -n "$SEGDUP_FILE"  ]] && echo "  segdup:  $SEGDUP_FILE"

        python3 -u /app/preprocessingMixer/generate_miXer_datasets.py \
            -j "$JSON_FILE" \
            -s "$SAMPLES_FILE" \
            -o "$OUTDIR" \
            ${XLR_FILE:+-x "$XLR_FILE"} \
            ${SEGDUP_FILE:+-sd "$SEGDUP_FILE"}
else
    # if no valid sub-command is provided, print usage
    echo "Usage: mixer.sh <sub-command> [options]"
    echo "Sub-commands:"
    echo "  preprocessing -j <json_file>   Run preprocessing with the specified JSON file"
    echo "  inference -j <json> -s <samples> -o <outdir> [-x <xlr>] [-sd <segdup>]"
    echo "Options:"
    echo "    -j, --json <file>          JSON configuration file for preprocessing"
    echo "    -s, --samples <file>       Samples file for inference"
    echo "    -o, --outdir <directory>   Output directory for inference results"
    echo "    -x, --xlr <file>           XLR file for inference (optional)"
    echo "    -sd, --segdup <file>       SegDup file for inference (optional)"
    exit 1
fi 
