#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: Merge eggNOG, MMseqs2 taxonomy and Salmon abundance results
#===========================================================================

#SBATCH --account=envi.envi
#SBATCH --job-name=annotation_merge
#SBATCH --partition=smp
#SBATCH --time=12:00:00
#SBATCH --qos=12h
#SBATCH --mem=240G
#SBATCH --cpus-per-task=16
#SBATCH --output=annotation_merge_%j.out

# Add your e-mail settings if SLURM notifications are needed.
# Example:
# #SBATCH --mail-type=END,FAIL
# #SBATCH --mail-user=your.email@institute.de


#===========================================================================
# VARIABLES
#===========================================================================

WORK=${PWD}

OUTDIR="output"

OUT_EGGNOG="out.eggnog"

OUT_MMSEQS="out.mmseqs2"

OUT_SALMON_MERGED="out.salmon_merged_paired"

OUT_ANNOTATION="out.annotation_merged"


R_SCRIPT="${WORK}/R_script/annotation_merge.R"


EGGNOG_FILE="${WORK}/${OUTDIR}/${OUT_EGGNOG}/non_redundant_prot_catalog_eggNOG.emapper.annotations"

ABUNDANCE_FILE="${WORK}/${OUTDIR}/${OUT_SALMON_MERGED}/output_all_raw_quant.sf"

TAXONOMY_FILE="${WORK}/${OUTDIR}/${OUT_MMSEQS}/non_redundant_prot_catalog.faa.result.tsv"


# Change this according to your metadata file.
# Set to NONE if metadata will not be used.
METADATA_FILE="${WORK}/metadata.tsv"


OUTPUT_DIR="${WORK}/${OUTDIR}/${OUT_ANNOTATION}"


mkdir -p ${OUTPUT_DIR}


#===========================================================================
# CHECK INPUT FILES
#===========================================================================

for FILE in \
    "${EGGNOG_FILE}" \
    "${ABUNDANCE_FILE}" \
    "${TAXONOMY_FILE}"
do

    if [[ ! -f "${FILE}" ]]; then
        echo "Input file not found: ${FILE}"
        exit 1
    fi

done


if [[ "${METADATA_FILE}" != "NONE" && ! -f "${METADATA_FILE}" ]]; then
    echo "Metadata file not found: ${METADATA_FILE}"
    exit 1
fi


#===========================================================================
# R ENVIRONMENT
#===========================================================================

# Change these according to your system.

CONDA_SH="${HOME}/miniforge3/etc/profile.d/conda.sh"

CONDA_ENV="functional_gap_r"


source ${CONDA_SH}

conda activate ${CONDA_ENV}


#===========================================================================
# ANNOTATION MERGE
#===========================================================================

Rscript ${R_SCRIPT} \
    ${EGGNOG_FILE} \
    ${ABUNDANCE_FILE} \
    ${TAXONOMY_FILE} \
    ${METADATA_FILE} \
    ${OUTPUT_DIR}


conda deactivate


echo "============================================================"
echo "Annotation merge finished."
echo "Output: ${OUTPUT_DIR}"
echo "============================================================"
