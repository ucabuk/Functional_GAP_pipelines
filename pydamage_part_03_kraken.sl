#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: Part 3 - Kraken2 classification of assemblies
#===========================================================================

#SBATCH --account=envi.envi
#SBATCH --job-name=kraken2_assembly
#SBATCH --partition=fat
#SBATCH --time=48:00:00
#SBATCH --mem=600G
#SBATCH --qos=48h
#SBATCH --cpus-per-task=64
#SBATCH --output=kraken2_%A_%a.out

# Set the array according to the number of samples.
# Example: #SBATCH --array=1-38%6

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

OUT_MEGAHIT="out.megahit"


END_MERGED="_fastp_merged_R2.fq.gz"

cd ${WORK}/${OUTDIR}/${OUT_FASTP}

FILE_MERGED=$(ls *${END_MERGED} | sed -n ${SLURM_ARRAY_TASK_ID}p)

SAMPLE_ID=${FILE_MERGED%${END_MERGED}}

cd ${WORK}


echo "============================================================"
echo "Sample: ${SAMPLE_ID}"
echo "============================================================"


CPU=${SLURM_CPUS_PER_TASK}


module load kraken2/2.1.3


#===========================================================================
# KRAKEN2 FOR ASSEMBLIES
#===========================================================================

DB="/albedo/work/projects/p_biodiv_dbs/nt_2022_10_db"

# Please do not change the confidence level. High conf. level gives nothing on assemblies...
CONFIDENCE="0"


if [[ -f ${WORK}/${OUTDIR}/${OUT_MEGAHIT}/${SAMPLE_ID}_conf${CONFIDENCE}_contig.kraken ]]; then

    echo "Removing existing Kraken file for ${SAMPLE_ID}..."

    rm -f \
        ${WORK}/${OUTDIR}/${OUT_MEGAHIT}/${SAMPLE_ID}_conf${CONFIDENCE}_contig.kraken \
        ${WORK}/${OUTDIR}/${OUT_MEGAHIT}/${SAMPLE_ID}_conf${CONFIDENCE}_contig.kraken.report

fi


srun kraken2 \
    --confidence ${CONFIDENCE} \
    --db ${DB} \
    ${WORK}/${OUTDIR}/${OUT_MEGAHIT}/${SAMPLE_ID}/final.contigs.fa \
    --threads ${CPU} \
    --output ${WORK}/${OUTDIR}/${OUT_MEGAHIT}/${SAMPLE_ID}_conf${CONFIDENCE}_contig.kraken \
    --report ${WORK}/${OUTDIR}/${OUT_MEGAHIT}/${SAMPLE_ID}_conf${CONFIDENCE}_contig.kraken.report


awk \
    -v prefix="$SAMPLE_ID" \
    '{print prefix, $0}' \
    OFS="\t" \
    ${WORK}/${OUTDIR}/${OUT_MEGAHIT}/${SAMPLE_ID}_conf${CONFIDENCE}_contig.kraken \
    | cut -f1-4 \
    > ${WORK}/${OUTDIR}/${OUT_MEGAHIT}/${SAMPLE_ID}_conf${CONFIDENCE}_contig_added_sample_name.kraken


echo "Kraken2 analysis completed for ${SAMPLE_ID}."


module unload kraken2/2.1.3
