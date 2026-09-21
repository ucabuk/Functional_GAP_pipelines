#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: Merge paired and merged Salmon abundance results
#===========================================================================

#SBATCH --account=envi.envi
#SBATCH --job-name=salmon_merge
#SBATCH --partition=smp
#SBATCH --time=08:00:00
#SBATCH --qos=12h
#SBATCH --mem=60G
#SBATCH --cpus-per-task=8
#SBATCH --output=salmon_merge_%j.out

# Add your e-mail settings if SLURM notifications are needed.
# Example:
# #SBATCH --mail-type=END,FAIL
# #SBATCH --mail-user=your.email@institute.de


#===========================================================================
# VARIABLES
#===========================================================================

WORK=${PWD}

OUTDIR="output"

OUT_SALMON="out.salmon"

OUT_SALMON_MERGED="out.salmon_merged_paired"

PYTHON_SCRIPT="${WORK}/Python_script/sum_up_qc_merged_paired.py"


mkdir -p ${OUTDIR}/${OUT_SALMON_MERGED}


#===========================================================================
# REMOVE OLD MERGED OUTPUTS
#===========================================================================

rm -f \
    ${OUTDIR}/${OUT_SALMON_MERGED}/${OUTDIR}_all_raw_quant.sf \
    ${OUTDIR}/${OUT_SALMON_MERGED}/${OUTDIR}_gene_quant.raw.count.len \
    ${OUTDIR}/${OUT_SALMON_MERGED}/${OUTDIR}_all_tpm_quant.sf


#===========================================================================
# COMBINE PAIRED AND MERGED SALMON RESULTS
#===========================================================================

srun python3 \
    ${PYTHON_SCRIPT} \
    ${OUTDIR}/${OUT_SALMON} \
    ${OUTDIR}/${OUT_SALMON_MERGED}


#===========================================================================
# PREPARE QUANT DIRECTORIES
#===========================================================================

for file in ${OUTDIR}/${OUT_SALMON_MERGED}/*.sf
do

    filename=$(basename "$file")

    ID="${filename%_merged_paired.quant.sf}"

    mkdir -p ${OUTDIR}/${OUT_SALMON_MERGED}/${ID}

    mv "$file" \
        ${OUTDIR}/${OUT_SALMON_MERGED}/${ID}/quant.sf

done


#===========================================================================
# MERGE SALMON OUTPUTS
#===========================================================================

module load salmon/1.10.2


srun salmon quantmerge \
    --quants ${OUTDIR}/${OUT_SALMON_MERGED}/* \
    --output ${OUTDIR}/${OUTDIR}_all_tpm_quant.sf


srun salmon quantmerge \
    --quants ${OUTDIR}/${OUT_SALMON_MERGED}/* \
    --column numreads \
    --output ${OUTDIR}/${OUTDIR}_all_raw_quant.sf


srun salmon quantmerge \
    --quants ${OUTDIR}/${OUT_SALMON_MERGED}/* \
    --column len \
    --output ${OUTDIR}/${OUTDIR}_gene_quant.raw.count.len


mv ${OUTDIR}/${OUTDIR}_all_tpm_quant.sf \
    ${OUTDIR}/${OUT_SALMON_MERGED}/.


mv ${OUTDIR}/${OUTDIR}_all_raw_quant.sf \
    ${OUTDIR}/${OUT_SALMON_MERGED}/.


mv ${OUTDIR}/${OUTDIR}_gene_quant.raw.count.len \
    ${OUTDIR}/${OUT_SALMON_MERGED}/.


module unload salmon/1.10.2


echo "============================================================"
echo "Salmon results merged."
echo "============================================================"
