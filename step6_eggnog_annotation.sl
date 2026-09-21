#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: Functional annotation of the non-redundant protein catalog with eggNOG-mapper
#===========================================================================

#SBATCH --account=envi.envi
#SBATCH --job-name=eggnog
#SBATCH --partition=smp
#SBATCH --time=48:00:00
#SBATCH --qos=48h
#SBATCH --cpus-per-task=32
#SBATCH --output=eggnog_%j.out

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

OUT_EGGNOG="out.eggnog"

OUT_TMP="tmp_eggnog"


IN="non_redundant_prot_catalog.faa"

OUT="non_redundant_prot_catalog_eggNOG"


CPU=${SLURM_CPUS_PER_TASK}


# Change this path according to the eggNOG database location on your system.
EGGNOG_DB_DIR="/albedo/work/projects/p_cbd/emapperdb-5.0.2"

export EGGNOG_DATA_DIR=${EGGNOG_DB_DIR}


EGGNOG_MODE="diamond"


#===========================================================================
# PREPARE OUTPUT
#===========================================================================

mkdir -p ${WORK}/${OUTDIR}/${OUT_EGGNOG}

mkdir -p ${WORK}/${OUTDIR}/${OUT_EGGNOG}/${OUT_TMP}


#===========================================================================
# EGGNOG-MAPPER
#===========================================================================

module load eggnog-mapper/2.1.12


srun emapper.py \
    -m ${EGGNOG_MODE} \
    -i ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/${IN} \
    -o ${OUT} \
    --output_dir ${WORK}/${OUTDIR}/${OUT_EGGNOG} \
    --temp_dir ${WORK}/${OUTDIR}/${OUT_EGGNOG}/${OUT_TMP} \
    --cpu ${CPU} 


module unload eggnog-mapper/2.1.12


echo "============================================================"
echo "eggNOG-mapper annotation finished."
echo "============================================================"
