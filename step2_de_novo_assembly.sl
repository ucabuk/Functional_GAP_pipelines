#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: Error correction with Tadpole followed by MEGAHIT assembly
#===========================================================================

#SBATCH --account=envi.envi
#SBATCH --job-name=assembly
#SBATCH --partition=fat
#SBATCH --time=48:00:00
#SBATCH --qos=48h
#SBATCH --cpus-per-task=64
#SBATCH --mem=450G
#SBATCH --output=megahit_%A_%a.out

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

OUT_FASTP="out.fastp"

OUT_TADPOLE="out.tadpole"

OUT_MEGAHIT="out.megahit"


END_R1="_fastp_R1.fq.gz"

END_R2="_fastp_R2.fq.gz"

END_MERGED="_fastp_merged_R2.fq.gz"


CPU=${SLURM_CPUS_PER_TASK}


#===========================================================================
# PREPARE SAMPLE
#===========================================================================

mkdir -p ${OUTDIR}/${OUT_TADPOLE}

mkdir -p ${OUTDIR}/${OUT_MEGAHIT}


cd ${OUTDIR}/${OUT_FASTP}

FILE_MERGED=$(ls *${END_MERGED} | sed -n ${SLURM_ARRAY_TASK_ID}p)

ID=${FILE_MERGED%${END_MERGED}}

cd ${WORK}


TMP="tmp_tadpole/${ID}"

mkdir -p ${TMP}


echo "============================================================"
echo "Sample: ${ID}"
echo "============================================================"


#===========================================================================
# TADPOLE ERROR CORRECTION
#===========================================================================

module load bbmap/39.01


srun tadpole.sh \
    in=${OUTDIR}/${OUT_FASTP}/${FILE_MERGED} \
    out=${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc${END_MERGED} \
    mode=correct \
    k=50


srun tadpole.sh \
    in=${OUTDIR}/${OUT_FASTP}/${ID}${END_R1} \
    in2=${OUTDIR}/${OUT_FASTP}/${ID}${END_R2} \
    out=${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc${END_R1} \
    out2=${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc${END_R2} \
    mode=correct \
    k=50


module unload bbmap/39.01


echo "============================================================"
echo "Tadpole correction finished."
echo "============================================================"


#===========================================================================
# MEGAHIT ASSEMBLY
#===========================================================================

module load megahit/1.2.9


echo "Assembly is performed with MEGAHIT."


srun megahit \
    --presets meta-large \
    -m 1 \
    -t ${CPU} \
    --tmp-dir ${TMP} \
    --min-contig-len 300 \
    -1 ${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc${END_R1} \
    -2 ${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc${END_R2} \
    -r ${OUTDIR}/${OUT_TADPOLE}/${ID}_tadpole_ecc${END_MERGED} \
    -o ${OUTDIR}/${OUT_MEGAHIT}/${ID}


module unload megahit/1.2.9


echo "============================================================"
echo "MEGAHIT assembly finished."
echo "============================================================"
