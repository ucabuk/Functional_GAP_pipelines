#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: Prodigal gene prediction from MEGAHIT assemblies
#===========================================================================

#SBATCH --account=envi.envi
#SBATCH --job-name=prodigal
#SBATCH --partition=smp
#SBATCH --time=10:00:00
#SBATCH --qos=12h
#SBATCH --cpus-per-task=1
#SBATCH --output=prodigal_%A_%a.out

# Set the array according to the number of samples.
# Example: #SBATCH --array=1-33%3

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

OUT_MEGAHIT="out.megahit"

OUT_PRODIGAL="out.prodigal"


END_MERGED="_tadpole_ecc_fastp_merged_R2.fq.gz"


#===========================================================================
# PREPARE SAMPLE
#===========================================================================

cd ${OUTDIR}/${OUT_TADPOLE}

FILE_MERGED=$(ls *${END_MERGED} | sed -n ${SLURM_ARRAY_TASK_ID}p)

ID=${FILE_MERGED%${END_MERGED}}

cd ${WORK}


mkdir -p ${OUTDIR}/${OUT_PRODIGAL}


echo "============================================================"
echo "Sample: ${ID}"
echo "============================================================"


#===========================================================================
# PRODIGAL GENE PREDICTION
#===========================================================================

module load prodigal/2.6.3


srun prodigal \
    -g 11 \
    -p meta \
    -i ${OUTDIR}/${OUT_MEGAHIT}/${ID}/final.contigs.fa \
    -a ${OUTDIR}/${OUT_PRODIGAL}/${ID}.faa \
    -d ${OUTDIR}/${OUT_PRODIGAL}/${ID}.fna \
    -f gff \
    -o ${OUTDIR}/${OUT_PRODIGAL}/${ID}.gff


module unload prodigal/2.6.3


echo "============================================================"
echo "Prodigal gene prediction finished."
echo "============================================================"
