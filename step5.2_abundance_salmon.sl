#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: Salmon gene abundance estimation
#===========================================================================

#SBATCH --account=envi.envi
#SBATCH --job-name=salmon
#SBATCH --partition=smp
#SBATCH --time=48:00:00
#SBATCH --qos=48h
#SBATCH --mem=220G
#SBATCH --cpus-per-task=64
#SBATCH --output=salmon_%A_%a.out

# Set the array according to the number of samples.
# Example: #SBATCH --array=1-33%6

# Add your e-mail settings if SLURM notifications are needed.
# Example:
# #SBATCH --mail-type=END,FAIL
# #SBATCH --mail-user=your.email@institute.de


#===========================================================================
# VARIABLES
#===========================================================================

WORK=${PWD}

OUTDIR="output"

OUT_TADPOLE="out.tadpole"

OUT_SALMON="out.salmon"


SALMON_INDEX="non_redundant_pCDS_catalog.index"

END_MERGED="_tadpole_ecc_fastp_merged_R2.fq.gz"


CPU=${SLURM_CPUS_PER_TASK}


mkdir -p ${OUTDIR}/${OUT_SALMON}


#===========================================================================
# PREPARE SAMPLE
#===========================================================================

cd ${WORK}/${OUTDIR}/${OUT_TADPOLE}

FILE_MERGED=$(ls *${END_MERGED} | sed -n ${SLURM_ARRAY_TASK_ID}p)

ID=${FILE_MERGED%${END_MERGED}}

cd ${WORK}


echo "============================================================"
echo "Sample: ${ID}"
echo "============================================================"


#===========================================================================
# SALMON QUANTIFICATION
#===========================================================================

module load salmon/1.10.2


srun salmon quant \
    -i ${OUTDIR}/${OUT_SALMON}/${SALMON_INDEX} \
    --meta \
    --libType A \
    -1 ${WORK}/${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc_fastp_R1.fq.gz \
    -2 ${WORK}/${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc_fastp_R2.fq.gz \
    --validateMappings \
    -o ${OUTDIR}/${OUT_SALMON}/${ID}_paired \
    --minScoreFraction 0.95


srun salmon quant \
    -i ${OUTDIR}/${OUT_SALMON}/${SALMON_INDEX} \
    --meta \
    --libType A \
    -r ${WORK}/${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc_fastp_merged_R2.fq.gz \
    --validateMappings \
    -o ${OUTDIR}/${OUT_SALMON}/${ID}_merged \
    --minScoreFraction 0.95


module unload salmon/1.10.2


echo "============================================================"
echo "Salmon quantification finished."
echo "============================================================"
