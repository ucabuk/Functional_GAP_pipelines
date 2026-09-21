#!/bin/bash
#SBATCH --account=envi.envi
#SBATCH --job-name=pydamage
#SBATCH --partition=smp
#SBATCH --time=48:00:00
#SBATCH --qos=48h
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --output=pydamage_%A_%a.out

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
OUT_BOWTIE="out.bowtie"
OUT_PYDAMAGE="out.pydamage"

END_MERGED="_fastp_merged_R2.fq.gz"

CPU=${SLURM_CPUS_PER_TASK}

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

BAM="${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}.merge_paired.sorted.bam"
PYDAMAGE_SAMPLE_DIR="${WORK}/${OUTDIR}/${OUT_PYDAMAGE}/${SAMPLE_ID}"

echo "Sample: ${SAMPLE_ID}"
echo "BAM   : ${BAM}"

# Input checks
#===================================================================

if [[ ! -s "${BAM}" ]]; then
    echo "ERROR: BAM file not found or empty:"
    echo "${BAM}"
    exit 1
fi

mkdir -p "${PYDAMAGE_SAMPLE_DIR}"

# pyDamage
#===================================================================

module load samtools/1.20
module load pydamage/1.0

if [[ ! -f "${BAM}.bai" ]]; then
    echo "BAM index not found. Creating index..."
    srun samtools index -@ "${CPU}" "${BAM}"
fi

echo "Running pyDamage 1.0..."
echo "C->T and G->A transitions are retained."
echo "The original pyDamage output schema is preserved."

srun pydamage \
    --outdir "${PYDAMAGE_SAMPLE_DIR}" \
    analyze "${BAM}" \
    -p "${CPU}" \
    --force

RESULT_FILE="${PYDAMAGE_SAMPLE_DIR}/pydamage_results.csv"

if [[ ! -s "${RESULT_FILE}" ]]; then
    echo "ERROR: pyDamage did not produce:"
    echo "${RESULT_FILE}"
    exit 1
fi

echo "pyDamage completed:"
echo "${RESULT_FILE}"

module unload pydamage/1.0
module unload samtools/1.20
