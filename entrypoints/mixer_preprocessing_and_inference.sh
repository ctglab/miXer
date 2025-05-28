#!/bin/bash
set -euo pipefail
mkdir -p /app/logs

LOGFILE="/app/logs/pipeline_$(date +'%Y%m%d_%H%M%S').log"
exec > >(tee -a "$LOGFILE") 2>&1

# Get Conda home directory
CONDA_HOME=$(conda info --base)
PREPR_SCRIPT_DIR="/app/preprocessingMixer"
INF_SCRIPT_DIR="/app/processing/scripts"
OUTDIR_EXCA="/app/exca2Out"
PREPARED_DATASET_DIR="/app/inference_ready_datasets"
HMM_PATH="${OUTDIR_SVM}/${KIT}_SVC"
HMM_RESOURCE_DIR="/app/processing/scripts/HMM_resources"
# Get exca2 config file basename
EXCA2CONFIG=$(basename "${CONFIG}")
source "${CONDA_HOME}/etc/profile.d/conda.sh"  # Load conda functions

###1. coverage:
echo "Calculating coverage"
conda activate ${CONDA_HOME}/envs/mixersing
### Allowing namespace creation for singularity file
python -u ${PREPR_SCRIPT_DIR}/computeAutosomalCvg.py -o ${OUTDIR_EXCA}/${KIT} -t ${TARGET} -cf ${CONFIG} -e ${KIT} -md ${SING_DIR_CONTAINER} -bd ${BAM_DIR_CONTAINER}

###2. merge samples to create XXXY:
echo "Merging samples if needed"
conda activate ${CONDA_HOME}/envs/mergedBam
if [ -e "${OUTDIR_EXCA}/${KIT}/${KIT}_autosomal_cvg.txt" ]; then
   python -u ${PREPR_SCRIPT_DIR}/mergeMixerTrainsamples.py -o ${OUTDIR_EXCA}/${KIT} -c ${OUTDIR_EXCA}/${KIT}/${KIT}_autosomal_cvg.txt -cf ${CONFIG} -th ${THREADS} -e ${KIT}
else
   python -u ${PREPR_SCRIPT_DIR}/mergeMixerTrainsamples.py -o ${OUTDIR_EXCA}/${KIT} -cf ${CONFIG} -th ${THREADS} -e ${KIT}
fi

###3. run EXCAVATOR2:
echo "Running EXCAVATOR2"
conda activate ${CONDA_HOME}/envs/mixersing
# Check for REF37 if it is empty 

if [ ! -z "$REF37" ]; then
	echo "REF37 is not empty, is it correct?: $REF37"
	python -u ${PREPR_SCRIPT_DIR}/run_excavator2singularity.py -o ${OUTDIR_EXCA}/${KIT} -b ${BAM_DIR_CONTAINER} -e ${KIT} -t ${TARGET} -cf ${OUTDIR_EXCA}/${KIT}/${EXCA2CONFIG}_forExca2 -r37 ${REF37} -m ${MAP} -g ${GAP} -cm ${CENTRO} -ch ${CHROM} -th ${THREADS} -ed ${SING_DIR_CONTAINER}
else
	python -u ${PREPR_SCRIPT_DIR}/run_excavator2singularity.py -o ${OUTDIR_EXCA}/${KIT} -b ${BAM_DIR_CONTAINER} -e ${KIT} -t ${TARGET} -cf ${OUTDIR_EXCA}/${KIT}/${EXCA2CONFIG}_forExca2 -r ${REF} -m ${MAP} -g ${GAP} -cm ${CENTRO} -ch ${CHROM} -th ${THREADS} -ed ${SING_DIR_CONTAINER}
fi


###4. create miXer datasets:
echo "Creating miXer datasets"
conda activate ${CONDA_HOME}/envs/mixerPre
if [ -z "$PREMADE_CONTROL_RDATA_PATH" ]; then
	python -u ${PREPR_SCRIPT_DIR}/generate_miXer_datasets.py -o ${PREPARED_DATASET_DIR}/${KIT} -t ${TARGET} -r ${REF} -m ${MAP} -p ${PAR} -g ${GAP} -cm ${CENTRO} -x ${XLR} -s ${SEG} -f ${OUTDIR_EXCA}/${KIT}/${KIT}_excavator2_output/DataAnalysis_w50k/Control/RCNorm/Control.NRC.RData -e ${KIT} -cf ${OUTDIR_EXCA}/${KIT}/${EXCA2CONFIG}_forExca2 -n ${OUTDIR_EXCA}/${KIT}/${KIT}_excavator2_output/DataPrepare_w50k/*/RCNorm/*RData
else
	echo "PREMADE_CONTROL_RDATA_PATH is not empty, is it correct?: $PREMADE_CONTROL_RDATA_PATH"
	python -u ${PREPR_SCRIPT_DIR}/generate_miXer_datasets.py -o ${PREPARED_DATASET_DIR}/${KIT} -t ${TARGET} -r ${REF} -m ${MAP} -p ${PAR} -g ${GAP} -cm ${CENTRO} -x ${XLR} -s ${SEG} -f ${PREMADE_CONTROL_RDATA_PATH} -e ${KIT} -cf ${OUTDIR_EXCA}/${KIT}/${EXCA2CONFIG}_forExca2 -n ${OUTDIR_EXCA}/${KIT}/${KIT}_excavator2_output/DataPrepare_w50k/*/RCNorm/*RData
fi

PREPARED_SVM_DIR=$(ls -td "${PREPARED_DATASET_DIR}/${KIT}/${KIT}_datasets_testing_"* | head -n 1)

# Check if a folder was found
if [ -n "$PREPARED_SVM_DIR" ]; then
    echo "Most recent folder for ${kit}: $PREPARED_SVM_DIR"
else
    echo "No folders found for ${kit}"
fi

###5. run SVM inference 
echo "Running SVM inference"
conda activate ${CONDA_HOME}/envs/miXer_ml
python ${INF_SCRIPT_DIR}/miXer_inference.py -tst ${PREPARED_SVM_DIR} -tst_exp_tags ${KIT} -output ${OUTDIR_SVM} -ntr ${THREADS}

###6. run HMM filtering
echo "Running HMM filtering"
conda activate ${CONDA_HOME}/envs/HMMR
Rscript ${INF_SCRIPT_DIR}/miXe.R -w ${INF_SCRIPT_DIR} -D ${HMM_PATH} -r ${PAR} -o ${OUTDIR_HMM} -t ${THREADS} -z ${HMM_RESOURCE_DIR}

###7. run WCF writer
echo "Writing VCFs"
python ${INF_SCRIPT_DIR}/Vcf_maker.py -bd ${OUTDIR_HMM} -od ${OUTDIR_VCF} -r ${REF_STRING}


