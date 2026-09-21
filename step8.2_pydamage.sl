#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: Part 2 - run pyDamage on Bowtie2 bam files.
#===========================================================================

#SBATCH --account=envi.envi
#SBATCH --job-name=pydamage
#SBATCH --partition=smp
#SBATCH --time=48:00:00
#SBATCH --mem=240G
#SBATCH --qos=48h
#SBATCH --cpus-per-task=64
#SBATCH --output=pydamage_%A_%a.out

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


module load samtools/1.20
module load pydamage/1.0


BAM=${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}.merge_paired.sorted.bam


#===========================================================================
# CHECK BAM INDEX
#===========================================================================

if [[ ! -f ${BAM}.bai ]]; then

    echo "BAM index not found. Creating index..."

    srun samtools index \
        -@ ${CPU} \
        ${BAM}

fi


#===========================================================================
# PYDAMAGE
#===========================================================================

srun pydamage \
    --outdir ${WORK}/${OUTDIR}/${OUT_PYDAMAGE}/${SAMPLE_ID} \
    analyze \
    ${BAM} \
    -p ${CPU} \
    --force


#===========================================================================
# ADD SAMPLE NAME TO PYDAMAGE OUTPUT
#===========================================================================

awk \
    -F, \
    -v OFS=, \
    -v prefix="${SAMPLE_ID}" \
    'NR==1 {
        print "sample_name," $0;
        next
    }
    {
        print prefix, $0
    }' \
    "${WORK}/${OUTDIR}/${OUT_PYDAMAGE}/${SAMPLE_ID}/pydamage_results.csv" \
    > \
    "${WORK}/${OUTDIR}/${OUT_PYDAMAGE}/${SAMPLE_ID}_name_added_pydamage_result.csv"


module unload pydamage/1.0


echo "============================================================"
echo "pyDamage finished."
echo "============================================================"
