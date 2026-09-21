#!/bin/bash
#SBATCH --account=envi.envi
#SBATCH --job-name=pydamage_merge
#SBATCH --partition=smp
#SBATCH --time=12:00:00
#SBATCH --qos=12h
#SBATCH --cpus-per-task=16
#SBATCH --mem=240G
#SBATCH --output=pydamage_merge_%j.out

# Add your e-mail settings if SLURM notifications are needed.
# Example:
# #SBATCH --mail-type=END,FAIL
# #SBATCH --mail-user=your.email@institute.de

set -euo pipefail

WORK=${PWD}

OUTDIR="output"
OUT_PYDAMAGE="out.pydamage"
OUT_KRAKEN="out.kraken"
OUT_FINAL="out.pydamage_annotated"

PYDAMAGE_DIR="${WORK}/${OUTDIR}/${OUT_PYDAMAGE}"
KRAKEN_DIR="${WORK}/${OUTDIR}/${OUT_KRAKEN}"
FINAL_DIR="${WORK}/${OUTDIR}/${OUT_FINAL}"

R_SCRIPT="${WORK}/R_script/pydamage_annotation_merge.R"

# Optional metadata.
# Use NONE when no metadata file should be added.
METADATA_FILE="NONE"
# Example:
# METADATA_FILE="${WORK}/metadata.tsv"

# NCBI taxonomy resources are downloaded automatically on first use.
NCBI_TAX_DIR="${WORK}/resources/ncbi_taxonomy"
NCBI_ARCHIVE="${NCBI_TAX_DIR}/new_taxdump.tar.gz"
NCBI_MD5="${NCBI_TAX_DIR}/new_taxdump.tar.gz.md5"
NCBI_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
NCBI_MD5_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz.md5"

RANKED_LINEAGE="${NCBI_TAX_DIR}/rankedlineage.dmp"
MERGED_DMP="${NCBI_TAX_DIR}/merged.dmp"

CONDA_SH="${HOME}/miniforge3/etc/profile.d/conda.sh"
CONDA_ENV="functional_gap_r"

# Input checks
#===================================================================

if [[ ! -d "${PYDAMAGE_DIR}" ]]; then
    echo "ERROR: pyDamage output directory not found:"
    echo "${PYDAMAGE_DIR}"
    exit 1
fi

if [[ ! -d "${KRAKEN_DIR}" ]]; then
    echo "ERROR: Kraken2 output directory not found:"
    echo "${KRAKEN_DIR}"
    exit 1
fi

if [[ ! -f "${R_SCRIPT}" ]]; then
    echo "ERROR: R script not found:"
    echo "${R_SCRIPT}"
    exit 1
fi

mkdir -p "${FINAL_DIR}"
mkdir -p "${NCBI_TAX_DIR}"

# NCBI taxonomy
#===================================================================

if [[ ! -s "${RANKED_LINEAGE}" || ! -s "${MERGED_DMP}" ]]; then

    echo "NCBI taxonomy files were not found."
    echo "Downloading the latest NCBI new_taxdump archive..."

    rm -f "${NCBI_ARCHIVE}" "${NCBI_MD5}"

    if command -v curl >/dev/null 2>&1; then
        curl -fL --retry 3 "${NCBI_URL}" -o "${NCBI_ARCHIVE}"
        curl -fL --retry 3 "${NCBI_MD5_URL}" -o "${NCBI_MD5}"

    elif command -v wget >/dev/null 2>&1; then
        wget --tries=3 "${NCBI_URL}" -O "${NCBI_ARCHIVE}"
        wget --tries=3 "${NCBI_MD5_URL}" -O "${NCBI_MD5}"

    else
        echo "ERROR: curl or wget is required to download NCBI taxonomy."
        exit 1
    fi

    if command -v md5sum >/dev/null 2>&1 && [[ -s "${NCBI_MD5}" ]]; then
        echo "Checking NCBI archive checksum..."
        (
            cd "${NCBI_TAX_DIR}"
            md5sum -c "$(basename "${NCBI_MD5}")"
        )
    else
        echo "WARNING: md5sum is unavailable. Checksum validation was skipped."
    fi

    echo "Extracting rankedlineage.dmp and merged.dmp..."

    tar -xzf "${NCBI_ARCHIVE}" \
        -C "${NCBI_TAX_DIR}" \
        rankedlineage.dmp merged.dmp

    rm -f "${NCBI_ARCHIVE}" "${NCBI_MD5}"

else
    echo "Using existing NCBI taxonomy files:"
    echo "${RANKED_LINEAGE}"
    echo "${MERGED_DMP}"
fi

if [[ ! -s "${RANKED_LINEAGE}" || ! -s "${MERGED_DMP}" ]]; then
    echo "ERROR: NCBI taxonomy extraction failed."
    exit 1
fi

# R environment
#===================================================================

if [[ ! -f "${CONDA_SH}" ]]; then
    echo "ERROR: Conda initialization script not found:"
    echo "${CONDA_SH}"
    exit 1
fi

source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

# Merge pyDamage, Kraken2, NCBI taxonomy and optional metadata
#===================================================================

Rscript "${R_SCRIPT}" \
    "${PYDAMAGE_DIR}" \
    "${KRAKEN_DIR}" \
    "${RANKED_LINEAGE}" \
    "${MERGED_DMP}" \
    "${METADATA_FILE}" \
    "${FINAL_DIR}"

conda deactivate

echo "Finished."
echo "Final output directory:"
echo "${FINAL_DIR}"
