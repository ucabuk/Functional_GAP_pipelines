#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: map reads to the assemblies
#===========================================================================

#SBATCH --account=envi.envi
#SBATCH --job-name=pydamage_bt2
#SBATCH --partition=smp
#SBATCH --time=48:00:00
#SBATCH --mem=240G
#SBATCH --qos=48h
#SBATCH --cpus-per-task=64
#SBATCH --output=pydamage_bt2_%A_%a.out

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

OUT_BOWTIE="out.bowtie"

OUT_PYDAMAGE="out.pydamage"


mkdir -p ${WORK}/${OUTDIR}/${OUT_BOWTIE}
mkdir -p ${WORK}/${OUTDIR}/${OUT_PYDAMAGE}


END_MERGED="_fastp_merged_R2.fq.gz"

cd ${WORK}/${OUTDIR}/${OUT_FASTP}

FILE_MERGED=$(ls *${END_MERGED} | sed -n ${SLURM_ARRAY_TASK_ID}p)

SAMPLE_ID=${FILE_MERGED%${END_MERGED}}

cd ${WORK}


echo "============================================================"
echo "Sample: ${SAMPLE_ID}"
echo "============================================================"


CPU=${SLURM_CPUS_PER_TASK}


# Check the available versions on the server if needed.
module load bowtie2
module load samtools/1.20
module load bamtools/2.5.2


REF=${WORK}/${OUTDIR}/${OUT_MEGAHIT}/${SAMPLE_ID}/final.contigs.fa

BT2_INDEX=${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_contigs_index


#===========================================================================
# INPUT READS
#===========================================================================

MERGED_READS=${WORK}/${OUTDIR}/${OUT_FASTP}/${SAMPLE_ID}_fastp_merged_R2.fq.gz

PAIRED_R1=${WORK}/${OUTDIR}/${OUT_FASTP}/${SAMPLE_ID}_fastp_R1.fq.gz

PAIRED_R2=${WORK}/${OUTDIR}/${OUT_FASTP}/${SAMPLE_ID}_fastp_R2.fq.gz


#===========================================================================
# BOWTIE2 INDEX
#===========================================================================

echo "============================================================"
echo "Running Bowtie2"
echo "Reference: ${REF}"
echo "============================================================"


srun bowtie2-build \
    ${REF} \
    ${BT2_INDEX}


#===========================================================================
# MAP MERGED READS
#===========================================================================

srun bowtie2 \
    --end-to-end \
    --very-sensitive \
    -N 1 \
    -p ${CPU} \
    -x ${BT2_INDEX} \
    -U ${MERGED_READS} \
    -S ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_merged.sam


#===========================================================================
# MAP PAIRED-END READS
#===========================================================================

srun bowtie2 \
    --end-to-end \
    --very-sensitive \
    -N 1 \
    -p ${CPU} \
    -x ${BT2_INDEX} \
    -1 ${PAIRED_R1} \
    -2 ${PAIRED_R2} \
    -S ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_paired.sam


#===========================================================================
# SORT SAM -> BAM
#===========================================================================

srun samtools sort \
    -@ ${CPU} \
    -o ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_merged.sorted.bam \
    ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_merged.sam


srun samtools sort \
    -@ ${CPU} \
    -o ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_paired.sorted.bam \
    ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_paired.sam


rm -f \
    ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_merged.sam

rm -f \
    ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_paired.sam


#===========================================================================
# MERGE BAM FILES
#===========================================================================

srun bamtools merge \
    -in ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_merged.sorted.bam \
    -in ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_paired.sorted.bam \
    -out ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}.merge_paired.bam


rm -f \
    ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_merged.sorted.bam

rm -f \
    ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_paired.sorted.bam


#===========================================================================
# SORT AND INDEX FINAL BAM
#===========================================================================

srun samtools sort \
    -@ ${CPU} \
    -o ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}.merge_paired.sorted.bam \
    ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}.merge_paired.bam


srun samtools index \
    -@ ${CPU} \
    ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}.merge_paired.sorted.bam


rm -f \
    ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}.merge_paired.bam


echo "============================================================"
echo "Bowtie2 mapping finished."
echo "============================================================"
