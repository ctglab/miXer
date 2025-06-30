#!/bin/bash

# Read variables from JSON configuration file

# Check if the configuration file path is provided as a command line argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_path>"
    exit 1
fi

# Set the CONFIG_PATH from the command line argument
CONFIG_PATH="$1"

# Check if the configuration file exists
if [ ! -f "$CONFIG_PATH" ]; then
    echo "Error: Configuration file not found - $CONFIG_PATH"
    exit 1
fi

eval "$(jq -r '@sh "KIT=\(.kit) THREADS=\(.threads) MIXER_SUPPORT_DIR=\(.mixer_support_dir) PREPARED_DATASET_DIR=\(.prepared_dataset_dir) MIXER_RESOURCES_CONTAINER=\(.mixer_resources_container) MIXER_SUPPORT_CONTAINER=\(.mixer_support_container) PREPARED_DATASET_CONTAINER=\(.prepared_dataset_container) PAR=\(.par) REF=\(.ref) REF37=\(.ref37) OUTDIR_SVM=\(.outdir_svm) OUTDIR_HMM=\(.outdir_hmm) OUTDIR_VCF=\(.outdir_vcf) OUTDIR_SVM_HOST=\(.outdir_svm_host) OUTDIR_HMM_HOST=\(.outdir_hmm_host) OUTDIR_VCF_HOST=\(.outdir_vcf_host)"' "$CONFIG_PATH")"


# Creating in-image paths for files
PAR_PATH="${MIXER_SUPPORT_CONTAINER%/}/$(basename "$PAR")"

# Check if REF37 is an empty string
if [ -n "$REF37" ]; then
    REF_STRING="$REF37"
else
    REF_STRING="$REF"
fi

# Run the docker image
docker run -it \
               -e KIT="$KIT" \
               -v "$MIXER_SUPPORT_DIR":"$MIXER_SUPPORT_CONTAINER" \
               -v "$PREPARED_DATASET_DIR":"$PREPARED_DATASET_CONTAINER" \
               -e "PREPARED_DATASET_CONTAINER"="$PREPARED_DATASET_CONTAINER"\
               -e "THREADS=$THREADS" \
               -e "PAR=$PAR_PATH" \
               -e "REF_STRING=$REF_STRING" \
               -e "OUTDIR_SVM=$OUTDIR_SVM" \
               -e "OUTDIR_HMM=$OUTDIR_HMM" \
               -e "OUTDIR_VCF=$OUTDIR_VCF" \
               mixer_docker:latest /bin/bash -c "mkdir -p /app/$OUTDIR_SVM /app/$OUTDIR_HMM /app/$OUTDIR_VCF && source /app/entrypoints/mixer_full_inference.sh"


# Copy the created folders to the outside
mkdir -p "$OUTDIR_SVM_HOST"
docker cp $(docker ps -q -n=1):$OUTDIR_SVM "$OUTDIR_SVM_HOST"
mkdir -p "$OUTDIR_HMM_HOST"
docker cp $(docker ps -q -n=1):$OUTDIR_HMM "$OUTDIR_HMM_HOST"
mkdir -p "$OUTDIR_VCF_HOST"
docker cp $(docker ps -q -n=1):$OUTDIR_VCF "$OUTDIR_VCF_HOST"
# Stop the running container
docker stop $(docker ps -q -n=1)


