#!/bin/bash

# Read variables from JSON configuration file

# Usage:
#   ./run.sh <config_path> [DOCKER_IMAGE]
#   e.g. ./run_docker_preprocessing_and_inference.sh config.json
#        ./run_docker_preprocessing_and_inference.sh config.json myrepo/mixer:dev

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Usage: $0 <config_path> [DOCKER_IMAGE]"
    exit 1
fi

CONFIG_PATH="$1"
IMAGE="${2:-mixer_docker:latest}"

# Check if the configuration file exists
if [ ! -f "$CONFIG_PATH" ]; then
    echo "Error: Configuration file not found - $CONFIG_PATH"
    exit 1
fi

eval "$(jq -r '@sh "  EXP_ID=\(.exp_id)  MIXER_RESOURCES_DIR=\(.mixer_resources_dir)  MIXER_SUPPORT_DIR=\(.support_dir)  EXCAVATOR2_SUPPORT_DIR=\(.support_dir) FASTA_DIR=\(.fasta_dir)  BAM_DIR=\(.bam_dir)  SING_DIR=\(.sing_dir)  MAIN_OUTDIR_HOST=\(.main_outdir_host)  SAMPLELIST=\(.sample_list)  TARGET=\(.target)  REF=\(.ref)  THREADS=\(.threads)  MAP=\(.map)  GAP=\(.gap) CENTRO=\(.centro) CHROM=\(.chrom) PAR=\(.par) PREMADE_CONTROL_RDATA=\(.premade_control_rdata)"' "$CONFIG_PATH")"

MIXER_RESOURCES_CONTAINER="/app/mixer_resources/"
SUPPORT_CONTAINER="/app/support/"
FASTA_DIR_CONTAINER="/app/fasta_dir/"
BAM_DIR_CONTAINER="/app/bams/"
SING_DIR_CONTAINER="/app/singfiles/"
MAIN_OUTPUT_DIR_CONTAINER="/app/mixer_outputs/"

#DEBUG_OUT="/docker_debug_out"
# Creating in-image paths for files
MIXER_RESOURCES_CONFIG="${MIXER_RESOURCES_CONTAINER}${SAMPLELIST}"
MIXER_RESOURCES_TARGET="${MIXER_RESOURCES_CONTAINER}${TARGET}"
PAR_PATH="${SUPPORT_CONTAINER%/}/$(basename "$PAR")"
#XLR_PATH="${SUPPORT_CONTAINER%/}/$(basename "$XLR")"
#SEG_PATH="${SUPPORT_CONTAINER%/}/$(basename "$SEG")"
MAP_PATH="${SUPPORT_CONTAINER%/}/$(basename "$MAP")"
GAP_PATH="${SUPPORT_CONTAINER%/}/$(basename "$GAP")"
CENTRO_PATH="${SUPPORT_CONTAINER%/}/$(basename "$CENTRO")"
CHROM_PATH="${SUPPORT_CONTAINER%/}/$(basename "$CHROM")"
REF_PATH="${FASTA_DIR_CONTAINER%/}/$(basename "$REF")"
REF_STRING="$REF"

# Check if PREMADE_CONTROL_RDATA is empty before joining
if [ -n "$PREMADE_CONTROL_RDATA" ]; then
  PREMADE_CONTROL_RDATA_PATH="${MIXER_RESOURCES_CONTAINER%/}/$(basename "$PREMADE_CONTROL_RDATA")"
else
  PREMADE_CONTROL_RDATA_PATH=""
fi

# Run the docker image
docker run -it --userns=host --privileged \
               -e EXP_ID="$EXP_ID" \
               -v "$MIXER_RESOURCES_DIR":"$MIXER_RESOURCES_CONTAINER" \
               -v "$MIXER_SUPPORT_DIR":"$SUPPORT_CONTAINER" \
               -v "$EXCAVATOR2_SUPPORT_DIR":"$SUPPORT_CONTAINER" \
               -v "$FASTA_DIR":"$FASTA_DIR_CONTAINER" \
               -v "$BAM_DIR":"$BAM_DIR_CONTAINER" \
               -v "$SING_DIR":"$SING_DIR_CONTAINER" \
               -e "PREMADE_CONTROL_RDATA_PATH"="$PREMADE_CONTROL_RDATA_PATH" \
               -e "SING_DIR_CONTAINER=$SING_DIR_CONTAINER" \
               -e "BAM_DIR_CONTAINER=$BAM_DIR_CONTAINER" \
               -e "CONFIG=$MIXER_RESOURCES_CONFIG" \
               -e "TARGET=$MIXER_RESOURCES_TARGET" \
               -e "REF=$REF_PATH" \
               -e "REF_STRING=$REF_STRING" \
               -e "THREADS=$THREADS" \
               -e "MAP=$MAP_PATH" \
               -e "GAP=$GAP_PATH" \
               -e "CENTRO=$CENTRO_PATH" \
               -e "CHROM=$CHROM_PATH" \
               -e "PAR=$PAR_PATH" \
               -e "MAIN_OUTPUT_DIR_CONTAINER=$MAIN_OUTPUT_DIR_CONTAINER" \
               mixer_docker:latest /bin/bash -c "mkdir -p /app/$MAIN_OUTPUT_DIR_CONTAINER && source /app/entrypoints/mixer_preprocessing_and_inference_simplified.sh"

#Debug copies
#docker cp $(docker ps -q -n=1):/app/exca2Out "$DEBUG_OUT"
#docker cp $(docker ps -q -n=1):/app/inference_ready_datasets "$DEBUG_OUT"

# Copy the created folders to the outside
mkdir -p "$MAIN_OUTDIR_HOST"
docker cp $(docker ps -q -n=1):$MAIN_OUTPUT_DIR_CONTAINER "$MAIN_OUTDIR_HOST"

# Stop the running container
docker stop $(docker ps -q -n=1)
# Remove the container
docker container rm -f $(docker ps -q -n=1)



