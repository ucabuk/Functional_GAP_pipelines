#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: Build Salmon index from the non-redundant pCDS catalog
#===========================================================================

#SBATCH --account=envi.envi
#SBATCH --job-name=salmon_index
#SBATCH --partition=fat
#SBATCH --time=48:00:00
#SBATCH --qos=48h
#SBATCH --mem=450G
#SBATCH --cpus-per-task=64
#SBATCH --output=salmon_index_%j.out

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

OUT_SALMON="out.salmon"


OUT_NUCL="non_redundant_pCDS_catalog.fna"

SALMON_INDEX="non_redundant_pCDS_catalog.index"


CPU=${SLURM_CPUS_PER_TASK}


mkdir -p ${OUTDIR}/${OUT_SALMON}


#===========================================================================
# SALMON INDEX
#===========================================================================

module load salmon/1.10.2


if [[ -f ${OUTDIR}/${OUT_SALMON}/${SALMON_INDEX}/versionInfo.json ]]; then

    echo "Salmon index already exists:"
    echo "${OUTDIR}/${OUT_SALMON}/${SALMON_INDEX}"

else

    echo "Building Salmon index..."

    srun salmon index \
        -t ${OUTDIR}/${OUT_PRODIGAL}/${OUT_NUCL} \
        -i ${OUTDIR}/${OUT_SALMON}/${SALMON_INDEX} \
        -p ${CPU}

fi


module unload salmon/1.10.2


echo "============================================================"
echo "Salmon index finished."
echo "============================================================"
