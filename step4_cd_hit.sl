#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: CD-HIT clustering to construct a non-redundant gene catalog
#===========================================================================

#SBATCH --account=envi.envi
#SBATCH --job-name=cdhit
#SBATCH --partition=smp
#SBATCH --time=48:00:00
#SBATCH --qos=48h
#SBATCH --cpus-per-task=64
#SBATCH --output=cdhit_%j.out

# Add your e-mail settings if SLURM notifications are needed.
# Example:
# #SBATCH --mail-type=END,FAIL
# #SBATCH --mail-user=your.email@institute.de


#===========================================================================
# VARIABLES
#===========================================================================

WORK=${PWD}

OUTDIR="output"

OUT_PRODIGAL="out.prodigal"


IN="redundant_prot_catalog.faa"

OUT="non_redundant_prot_catalog.faa"


CPU=${SLURM_CPUS_PER_TASK}


#===========================================================================
# REMOVE OLD CATALOG FILES
#===========================================================================

rm -f \
    ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/redundant_prot_catalog.faa \
    ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/redundant_pCDS_catalog.fna \
    ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/non_redundant_prot_catalog.faa \
    ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/non_redundant_prot_catalog.faa.clstr \
    ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/non_redundant_prot_catalog.id \
    ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/non_redundant_pCDS_catalog.fna


#===========================================================================
# CREATE REDUNDANT CATALOG
#===========================================================================

cat ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/*.faa \
    > ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/redundant_prot_catalog.faa


cat ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/*.fna \
    > ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/redundant_pCDS_catalog.fna


sed -i 's/ //g' \
    ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/redundant_prot_catalog.faa


sed -i 's/ //g' \
    ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/redundant_pCDS_catalog.fna


#===========================================================================
# CD-HIT
#===========================================================================

module load cd-hit/4.8.1


srun cd-hit \
    -i ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/${IN} \
    -o ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/${OUT} \
    -n 5 \
    -c 0.95 \
    -G 0 \
    -M 0 \
    -d 0 \
    -aS 0.85 \
    -T ${CPU}


module unload cd-hit/4.8.1


#===========================================================================
# GET NON-REDUNDANT GENE IDS
#===========================================================================

grep '>' \
    ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/${OUT} \
    | cut -c 2- \
    > ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/non_redundant_prot_catalog.id


#===========================================================================
# CREATE NON-REDUNDANT PCDS CATALOG
#===========================================================================

module load seqtk/1.4


srun seqtk subseq \
    ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/redundant_pCDS_catalog.fna \
    ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/non_redundant_prot_catalog.id \
    > ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/non_redundant_pCDS_catalog.fna


module unload seqtk/1.4


echo "============================================================"
echo "Non-redundant gene catalog finished."
echo "============================================================"
