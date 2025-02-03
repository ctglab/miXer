#!/bin/bash
# Get Conda home directory
CONDA_HOME=$(conda info --base)
INF_SCRIPT_DIR="/app/processing/scripts"
HMM_RESOURCE_DIR="/app/processing/scripts/HMM_resources"
HMM_PATH="${OUTDIR_SVM}/${KIT}_SVC"
# Get exca2 config file basename
EXCA2CONFIG=$(basename "${CONFIG}")
source "${CONDA_HOME}/etc/profile.d/conda.sh"  # Load conda functions

###1. run SVM inference 
echo "Running SVM inference"
conda activate ${CONDA_HOME}/envs/miXer_ml
python ${INF_SCRIPT_DIR}/miXer_inference.py -tst ${PREPARED_DATASET_CONTAINER} -tst_exp_tags ${KIT} -output ${OUTDIR_SVM} -ntr ${THREADS}

###2. run HMM filtering
echo "Running HMM filtering"
conda activate ${CONDA_HOME}/envs/HMMR
Rscript ${INF_SCRIPT_DIR}/miXe.R -w ${INF_SCRIPT_DIR} -D ${HMM_PATH} -r ${PAR} -o ${OUTDIR_HMM} -t ${THREADS} -z ${HMM_RESOURCE_DIR}

###3. run VCF writer
echo "Writing VCFs"
python ${INF_SCRIPT_DIR}/Vcf_maker.py -bd ${OUTDIR_HMM} -od ${OUTDIR_VCF} -r ${REF_STRING}

