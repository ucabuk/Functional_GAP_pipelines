#!/bin/bash
#SBATCH --account=envi.envi
#SBATCH --job-name=kraken2
#SBATCH --partition=fat
#SBATCH --time=48:00:00
#SBATCH --qos=48h
#SBATCH --cpus-per-task=64
#SBATCH --mem=600G
#SBATCH --output=kraken2_%A_%a.out

# Set the array according to the number of samples.
# Example: #SBATCH --array=1-38%6
#
# Add your e-mail settings if SLURM notifications are needed.
# Example:
# #SBATCH --mail-type=END,FAIL
# #SBATCH --mail-user=your.email@institute.de

set -euo pipefail

WORK=${PWD}

OUTDIR="output"
OUT_FASTP="out.fastp"
OUT_MEGAHIT="out.megahit"
OUT_KRAKEN="out.kraken"

END_MERGED="_fastp_merged_R2.fq.gz"

CPU=${SLURM_CPUS_PER_TASK}

# Kraken2 database
# Change this path if a different Kraken2 database is used on your system.
KRAKEN_DB="/albedo/work/projects/p_biodiv_dbs/nt_2022_10_db"

# Please do not change the confidence level for the standard workflow.
CONFIDENCE="0"

# Identify sample
#===================================================================

mapfile -t MERGED_FILES < <(find "${WORK}/${OUTDIR}/${OUT_FASTP}" -maxdepth 1 -type f -name "*${END_MERGED}" | sort)

if [[ ${#MERGED_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No merged FASTQ files found in ${WORK}/${OUTDIR}/${OUT_FASTP}"
    exit 1
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set."
    echo "Submit this script as a SLURM array."
    exit 1
fi

if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > ${#MERGED_FILES[@]} )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} is outside the sample range 1-${#MERGED_FILES[@]}."
    exit 1
fi

MERGED_FILE=${MERGED_FILES[$((SLURM_ARRAY_TASK_ID - 1))]}
MERGED_NAME=$(basename "${MERGED_FILE}")
SAMPLE_ID=${MERGED_NAME%${END_MERGED}}

CONTIGS="${WORK}/${OUTDIR}/${OUT_MEGAHIT}/${SAMPLE_ID}/final.contigs.fa"

KRAKEN_DIR="${WORK}/${OUTDIR}/${OUT_KRAKEN}"
RAW_OUTPUT="${KRAKEN_DIR}/${SAMPLE_ID}_conf${CONFIDENCE}_contigs.kraken"
REPORT_OUTPUT="${KRAKEN_DIR}/${SAMPLE_ID}_conf${CONFIDENCE}_contigs.kraken.report"
TABLE_OUTPUT="${KRAKEN_DIR}/${SAMPLE_ID}_conf${CONFIDENCE}_contigs.kraken.tsv"

echo "Sample : ${SAMPLE_ID}"
echo "Contigs: ${CONTIGS}"

# Input checks
#===================================================================

if [[ ! -s "${CONTIGS}" ]]; then
    echo "ERROR: MEGAHIT contig file not found or empty:"
    echo "${CONTIGS}"
    exit 1
fi

if [[ ! -d "${KRAKEN_DB}" ]]; then
    echo "ERROR: Kraken2 database not found:"
    echo "${KRAKEN_DB}"
    exit 1
fi

mkdir -p "${KRAKEN_DIR}"

# Kraken2
#===================================================================

module load kraken2/2.1.3

srun kraken2 \
    --db "${KRAKEN_DB}" \
    --threads "${CPU}" \
    --confidence "${CONFIDENCE}" \
    --output "${RAW_OUTPUT}" \
    --report "${REPORT_OUTPUT}" \
    "${CONTIGS}"

# Create a normalized table for downstream merging.
# Kraken2 output columns:
# status, sequence ID, taxid, sequence length, LCA hit list
awk -v sample="${SAMPLE_ID}" '
BEGIN {
    FS=OFS="\t"
    print "sample_id","kraken_status","contig_id","kraken_taxid","sequence_length","kraken_hitlist"
}
{
    print sample,$1,$2,$3,$4,$5
}
' "${RAW_OUTPUT}" > "${TABLE_OUTPUT}"

echo "Kraken2 completed:"
echo "${TABLE_OUTPUT}"

module unload kraken2/2.1.3
