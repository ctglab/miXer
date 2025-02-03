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

eval "$(jq -r '@sh "KIT=\(.kit) CONFIG=\(.config) TARGET=\(.target) REF=\(.ref) REF37=\(.ref37) THREADS=\(.threads) MAP=\(.map) GAP=\(.gap) CENTRO=\(.centro) CHROM=\(.chrom) PAR=\(.par) XLR=\(.xlr) SEG=\(.seg) PREMADE_CONTROL_RDATA=\(.premade_control_rdata) OUTDIR_SVM=\(.outdir_svm) OUTDIR_HMM=\(.outdir_hmm) OUTDIR_VCF=\(.outdir_vcf) OUTDIR_SVM_HOST=\(.outdir_svm_host) OUTDIR_HMM_HOST=\(.outdir_hmm_host) OUTDIR_VCF_HOST=\(.outdir_vcf_host) MIXER_RESOURCES_DIR=\(.mixer_resources_dir) MIXER_SUPPORT_DIR=\(.mixer_support_dir) EXCAVATOR2_SUPPORT_DIR=\(.excavator2_support_dir) FASTA_DIR=\(.fasta_dir) BAM_DIR=\(.bam_dir) SING_DIR=\(.sing_dir) MIXER_RESOURCES_CONTAINER=\(.mixer_resources_container) MIXER_SUPPORT_CONTAINER=\(.mixer_support_container) EXCAVATOR2_SUPPORT_CONTAINER=\(.excavator2_support_container) FASTA_DIR_CONTAINER=\(.fasta_dir_container) BAM_DIR_CONTAINER=\(.bam_dir_container) SING_DIR_CONTAINER=\(.sing_dir_container)"' "$CONFIG_PATH")"

#DEBUG_OUT="/docker_debug_out"
# Creating in-image paths for files
MIXER_RESOURCES_CONFIG="${MIXER_RESOURCES_CONTAINER}${CONFIG}"
MIXER_RESOURCES_TARGET="${MIXER_RESOURCES_CONTAINER}${TARGET}"
PAR_PATH="${MIXER_SUPPORT_CONTAINER%/}/$(basename "$PAR")"
XLR_PATH="${MIXER_SUPPORT_CONTAINER%/}/$(basename "$XLR")"
SEG_PATH="${MIXER_SUPPORT_CONTAINER%/}/$(basename "$SEG")"
MAP_PATH="${EXCAVATOR2_SUPPORT_CONTAINER%/}/$(basename "$MAP")"
GAP_PATH="${EXCAVATOR2_SUPPORT_CONTAINER%/}/$(basename "$GAP")"
CENTRO_PATH="${EXCAVATOR2_SUPPORT_CONTAINER%/}/$(basename "$CENTRO")"
CHROM_PATH="${EXCAVATOR2_SUPPORT_CONTAINER%/}/$(basename "$CHROM")"
REF_PATH="${FASTA_DIR_CONTAINER%/}/$(basename "$REF")"


: '
echo "DEBUG_OUT: $DEBUG_OUT"
echo "MIXER_RESOURCES_CONFIG: $MIXER_RESOURCES_CONFIG"
echo "MIXER_RESOURCES_TARGET: $MIXER_RESOURCES_TARGET"
echo "PAR_PATH: $PAR_PATH"
echo "XLR_PATH: $XLR_PATH"
echo "SEG_PATH: $SEG_PATH"
echo "MAP_PATH: $MAP_PATH"
echo "GAP_PATH: $GAP_PATH"
echo "CENTRO_PATH: $CENTRO_PATH"
echo "CHROM_PATH: $CHROM_PATH"
echo "REF_PATH: $REF_PATH"
'

# Debugging inside the container script
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

# Run the docker image
docker run -it --userns=host --privileged \
               -e KIT="$KIT" \
               -v "$MIXER_RESOURCES_DIR":"$MIXER_RESOURCES_CONTAINER" \
               -v "$MIXER_SUPPORT_DIR":"$MIXER_SUPPORT_CONTAINER" \
               -v "$EXCAVATOR2_SUPPORT_DIR":"$EXCAVATOR2_SUPPORT_CONTAINER" \
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
               -e "REF37=$REF37_PATH" \
               -e "THREADS=$THREADS" \
               -e "MAP=$MAP_PATH" \
               -e "GAP=$GAP_PATH" \
               -e "CENTRO=$CENTRO_PATH" \
               -e "CHROM=$CHROM_PATH" \
               -e "PAR=$PAR_PATH" \
               -e "XLR=$XLR_PATH" \
               -e "SEG=$SEG_PATH" \
               -e "OUTDIR_SVM=$OUTDIR_SVM" \
               -e "OUTDIR_HMM=$OUTDIR_HMM" \
               -e "OUTDIR_VCF=$OUTDIR_VCF" \
               mixer_docker:latest /bin/bash -c "mkdir -p /app/$OUTDIR_SVM /app/$OUTDIR_HMM /app/$OUTDIR_VCF && source /app/entrypoints/mixer_preprocessing_and_inference.sh"

#Debug copies
#docker cp $(docker ps -q -n=1):/app/exca2Out "$DEBUG_OUT"
#docker cp $(docker ps -q -n=1):/app/inference_ready_datasets "$DEBUG_OUT"

# Copy the created folders to the outside
mkdir -p "$OUTDIR_SVM_HOST"
docker cp $(docker ps -q -n=1):/app/$OUTDIR_SVM "$OUTDIR_SVM_HOST"
mkdir -p "$OUTDIR_HMM_HOST"
docker cp $(docker ps -q -n=1):$OUTDIR_HMM "$OUTDIR_HMM_HOST"
mkdir -p "$OUTDIR_VCF_HOST"
docker cp $(docker ps -q -n=1):$OUTDIR_VCF "$OUTDIR_VCF_HOST"
# Stop the running container
docker stop $(docker ps -q -n=1)
# Remove the container
docker container rm -f $(docker ps -q -n=1)



