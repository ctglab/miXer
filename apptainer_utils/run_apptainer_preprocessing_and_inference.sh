#!/bin/bash

# Check usage
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_path>"
    exit 1
fi

CONFIG_PATH="$1"

# Check if config file exists
if [ ! -f "$CONFIG_PATH" ]; then
    echo "Error: Configuration file not found - $CONFIG_PATH"
    exit 1
fi

# Load variables from JSON config
eval "$(jq -r '@sh "EXP_ID=\(.exp_id) CONFIG=\(.config) TARGET=\(.target) REF=\(.ref) REF37=\(.ref37) THREADS=\(.threads) MAP=\(.map) GAP=\(.gap) CENTRO=\(.centro) CHROM=\(.chrom) PAR=\(.par) PREMADE_CONTROL_RDATA=\(.premade_control_rdata) MAIN_OUTDIR_HOST=\(.main_outdir_host) MIXER_RESOURCES_DIR=\(.mixer_resources_dir) MIXER_SUPPORT_DIR=\(.mixer_support_dir) EXCAVATOR2_SUPPORT_DIR=\(.excavator2_support_dir) FASTA_DIR=\(.fasta_dir) BAM_DIR=\(.bam_dir) SING_DIR=\(.sing_dir) MIXER_APPTAINER_SIF=\(.mixer_apptainer_sif)"' "$CONFIG_PATH")"

# Create per-run temp workspace - scripts will write temporary files here
TEMP_DIR="./temp_${EXP_ID}"
FINAL_OUTPUT_DIR="${MAIN_OUTDIR_HOST}/${EXP_ID}"

mkdir -p "$TEMP_DIR/exca2_output_${EXP_ID}" \
         "$TEMP_DIR/inference_ready_datasets_${EXP_ID}" \
         "$TEMP_DIR/svm_processed_output" \
         "$TEMP_DIR/mixer_windows" \
         "$TEMP_DIR/mixer_vcfs" \
         "$TEMP_DIR/mixer_outputs/${EXP_ID}" \
         "$FINAL_OUTPUT_DIR"

# Define container paths
MIXER_RESOURCES_CONTAINER="/app/mixer_resources/"
MIXER_SUPPORT_CONTAINER="/app/mixer_support/"
EXCAVATOR2_SUPPORT_CONTAINER="/app/excavator2_support/"
FASTA_DIR_CONTAINER="/app/fasta_dir/"
BAM_DIR_CONTAINER="/app/bams/"
SING_DIR_CONTAINER="/app/singfiles/"
MAIN_OUTPUT_DIR_CONTAINER="/app/mixer_outputs/"

# File paths inside the container
MIXER_RESOURCES_CONFIG="${MIXER_RESOURCES_CONTAINER}${CONFIG}"
MIXER_RESOURCES_TARGET="${MIXER_RESOURCES_CONTAINER}${TARGET}"
PAR_PATH="${MIXER_SUPPORT_CONTAINER%/}/$(basename "$PAR")"
MAP_PATH="${EXCAVATOR2_SUPPORT_CONTAINER%/}/$(basename "$MAP")"
GAP_PATH="${EXCAVATOR2_SUPPORT_CONTAINER%/}/$(basename "$GAP")"
CENTRO_PATH="${EXCAVATOR2_SUPPORT_CONTAINER%/}/$(basename "$CENTRO")"
CHROM_PATH="${EXCAVATOR2_SUPPORT_CONTAINER%/}/$(basename "$CHROM")"
REF_PATH="${FASTA_DIR_CONTAINER%/}/$(basename "$REF")"
REF_STRING="${REF37:-$REF}"

# Check if REF37 is empty before joining
# REF_STRING is used for VCF creation
if [ -n "$REF37" ]; then
  REF37_PATH="${FASTA_DIR_CONTAINER%/}/$(basename "$REF37")"
  REF_STRING="$REF37"
else
  REF37_PATH=""
  REF_STRING="$REF"
fi

# Check if PREMADE_CONTROL_RDATA is empty before joining
if [ -n "$PREMADE_CONTROL_RDATA" ]; then
  PREMADE_CONTROL_RDATA_PATH="${MIXER_RESOURCES_CONTAINER%/}/$(basename "$PREMADE_CONTROL_RDATA")"
else
  PREMADE_CONTROL_RDATA_PATH=""
fi

# Ensure output directory exists
if [ -d "$MAIN_OUTDIR_HOST" ]; then
    echo "Output directory exists: $MAIN_OUTDIR_HOST"
else
    echo "Output directory does not exist. Creating: $MAIN_OUTDIR_HOST"
    mkdir -p "$MAIN_OUTDIR_HOST"
    if [ $? -ne 0 ]; then
        echo "[ERROR] Failed to create output directory: $MAIN_OUTDIR_HOST"
        exit 1
    fi
fi

# Run Apptainer container
apptainer exec \
  --bind "$TEMP_DIR/exca2_output_${EXP_ID}":"/app/exca2_output_${EXP_ID}" \
  --bind "$TEMP_DIR/inference_ready_datasets_${EXP_ID}":"/app/inference_ready_datasets_${EXP_ID}" \
  --bind "$TEMP_DIR/svm_processed_output":"/app/svm_processed_output" \
  --bind "$TEMP_DIR/mixer_windows":"/app/mixer_windows" \
  --bind "$TEMP_DIR/mixer_vcfs":"/app/mixer_vcfs" \
  --bind "$TEMP_DIR/mixer_outputs/${EXP_ID}":"/app/mixer_outputs/${EXP_ID}" \
  --bind "$MIXER_RESOURCES_DIR":"$MIXER_RESOURCES_CONTAINER" \
  --bind "$MIXER_SUPPORT_DIR":"$MIXER_SUPPORT_CONTAINER" \
  --bind "$EXCAVATOR2_SUPPORT_DIR":"$EXCAVATOR2_SUPPORT_CONTAINER" \
  --bind "$FASTA_DIR":"$FASTA_DIR_CONTAINER" \
  --bind "$BAM_DIR":"$BAM_DIR_CONTAINER" \
  --bind "$SING_DIR":"$SING_DIR_CONTAINER" \
  --bind "$MAIN_OUTDIR_HOST":"$MAIN_OUTPUT_DIR_CONTAINER" \
  --env EXP_ID="$EXP_ID" \
  --env CONFIG="$MIXER_RESOURCES_CONFIG" \
  --env TARGET="$MIXER_RESOURCES_TARGET" \
  --env REF="$REF_PATH" \
  --env REF37="$REF37_PATH" \
  --env REF_STRING="$REF_STRING" \
  --env THREADS="$THREADS" \
  --env MAP="$MAP_PATH" \
  --env GAP="$GAP_PATH" \
  --env CENTRO="$CENTRO_PATH" \
  --env CHROM="$CHROM_PATH" \
  --env PAR="$PAR_PATH" \
  --env PREMADE_CONTROL_RDATA_PATH="$PREMADE_CONTROL_RDATA_PATH" \
  --env SING_DIR_CONTAINER="$SING_DIR_CONTAINER" \
  --env BAM_DIR_CONTAINER="$BAM_DIR_CONTAINER" \
  --env MAIN_OUTPUT_DIR_CONTAINER="$MAIN_OUTPUT_DIR_CONTAINER" \
  "$MIXER_APPTAINER_SIF" \
  /bin/bash /app/entrypoints/mixer_apptainer_preprocessing_and_inference.sh


echo "Copying results to final output dir: $FINAL_OUTPUT_DIR"

cp -r "$TEMP_DIR/exca2_output_${EXP_ID}" "$FINAL_OUTPUT_DIR/"
cp -r "$TEMP_DIR/mixer_windows" "$FINAL_OUTPUT_DIR/"
cp -r "$TEMP_DIR/mixer_vcfs" "$FINAL_OUTPUT_DIR/"

echo "Removing temporary directory: $TEMP_DIR"
rm -rf "$TEMP_DIR"