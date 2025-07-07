#!/bin/bash


# Get Conda home directory
CONDA_HOME=$(conda info --base)
PREPR_SCRIPT_DIR="/app/preprocessingMixer"
INF_SCRIPT_DIR="/app/processing/scripts"
PREPARED_DATASET_DIR="/app/inference_ready_datasets"
OUTDIR_EXCA="/app/exca2_output"
OUTDIR_SVM="/app/svm_processed_output"
OUTDIR_HMM="/app/mixer_windows"
OUTDIR_VCF="/app/mixer_vcfs"
HMM_PATH="${OUTDIR_SVM}/${EXP_ID}_SVC"
HMM_RESOURCE_DIR="/app/processing/scripts/HMM_resources"
# Get exca2 config file basename
EXCA2CONFIG=$(basename "${CONFIG}")
source "${CONDA_HOME}/etc/profile.d/conda.sh"  # Load conda functions

###0. Create internal config file for EXCAVATOR2 and check correctness
echo "Checking config"
conda activate ${CONDA_HOME}/envs/miXer_ml
python -u ${PREPR_SCRIPT_DIR}/check_config.py -cf ${CONFIG} -o ${OUTDIR_EXCA}_${EXP_ID} -bd ${BAM_DIR_CONTAINER}


###1. run EXCAVATOR2:
echo "Running EXCAVATOR2"
conda activate ${CONDA_HOME}/envs/mixersing

python -u ${PREPR_SCRIPT_DIR}/run_excavator2singularity.py -o ${OUTDIR_EXCA}_${EXP_ID} -b ${BAM_DIR_CONTAINER} -e ${EXP_ID} -t ${TARGET} -cf ${OUTDIR_EXCA}_${EXP_ID}/${EXCA2CONFIG}_forExca2 -r ${REF} -m ${MAP} -g ${GAP} -cm ${CENTRO} -ch ${CHROM} -th ${THREADS} -ed ${SING_DIR_CONTAINER}

###2. create miXer datasets:
echo "Creating miXer datasets"
conda activate ${CONDA_HOME}/envs/mixerPre
if [ -z "$PREMADE_CONTROL_RDATA_PATH" ]; then
	python -u ${PREPR_SCRIPT_DIR}/generate_miXer_datasets.py -o ${PREPARED_DATASET_DIR}_${EXP_ID} -t ${TARGET} -r ${REF} -m ${MAP} -p ${PAR} -g ${GAP} -cm ${CENTRO} -f ${OUTDIR_EXCA}_${EXP_ID}/${EXP_ID}_excavator2_output/DataAnalysis_w50k/Control/RCNorm/Control.NRC.RData -e ${EXP_ID} -cf ${OUTDIR_EXCA}_${EXP_ID}/${EXCA2CONFIG}_forExca2 -n ${OUTDIR_EXCA}_${EXP_ID}/${EXP_ID}_excavator2_output/DataPrepare_w50k/*/RCNorm/*RData
else
	echo "PREMADE_CONTROL_RDATA_PATH is not empty, is it correct?: $PREMADE_CONTROL_RDATA_PATH"
	python -u ${PREPR_SCRIPT_DIR}/generate_miXer_datasets.py -o ${PREPARED_DATASET_DIR}_${EXP_ID} -t ${TARGET} -r ${REF} -m ${MAP} -p ${PAR} -g ${GAP} -cm ${CENTRO} -f ${PREMADE_CONTROL_RDATA_PATH} -e ${EXP_ID} -cf ${OUTDIR_EXCA}_${EXP_ID}/${EXCA2CONFIG}_forExca2 -n ${OUTDIR_EXCA}_${EXP_ID}/${EXP_ID}_excavator2_output/DataPrepare_w50k/*/RCNorm/*RData
fi

PREPARED_SVM_DIR=$(ls -td "${PREPARED_DATASET_DIR}_${EXP_ID}/${EXP_ID}_datasets_testing_"* | head -n 1)

# Check if a folder was found
if [ -n "$PREPARED_SVM_DIR" ]; then
    echo "Most recent folder for ${kit}: $PREPARED_SVM_DIR"
else
    echo "No folders found for ${kit}"
fi

###3. run SVM inference 
echo "Running SVM inference"
conda activate ${CONDA_HOME}/envs/miXer_ml
python ${INF_SCRIPT_DIR}/miXer_inference.py -tst ${PREPARED_SVM_DIR} -tst_exp_tags ${EXP_ID} -output ${OUTDIR_SVM} -ntr ${THREADS}

###4. run HMM filtering
echo "Running HMM filtering"
conda activate ${CONDA_HOME}/envs/HMMR
Rscript ${INF_SCRIPT_DIR}/miXe.R -w ${INF_SCRIPT_DIR} -D ${HMM_PATH} -r ${PAR} -o ${OUTDIR_HMM} -t ${THREADS} -z ${HMM_RESOURCE_DIR}

###5. run WCF writer
echo "Writing VCFs"
python ${INF_SCRIPT_DIR}/Vcf_maker.py -td ${OUTDIR_HMM} -od ${OUTDIR_VCF} -r ${REF_STRING}

echo "Rounding up output files"
OUTDIR="${MAIN_OUTPUT_DIR_CONTAINER}${EXP_ID}"
mkdir -p "${OUTDIR}"

#Copying EXCAVATOR2 outputs in main output dir
cp -r "${OUTDIR_EXCA}_${EXP_ID}" "${OUTDIR}"
#Copying SVM-ready datasets in main output dir 
#cp -r "${PREPARED_SVM_DIR}" "${OUTDIR}"
#Copying SVM-processed datasets in main output dir
#cp -r "${HMM_PATH}" "${OUTDIR}"
#Copying HMM-ready datasets in main ouput dir
cp -r "${OUTDIR_HMM}" "${OUTDIR}"
#Copying VCF output in main output dir
cp -r "${OUTDIR_VCF}" "${OUTDIR}"