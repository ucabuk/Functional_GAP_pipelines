#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: pyDamage workflow - FASTP uncorrected reads + Bowtie2
#===========================================================================

#SBATCH --account=envi.envi
#SBATCH --job-name=pydamage_bt2
#SBATCH --partition=smp
#SBATCH --time=48:00:00
#SBATCH --mem=240G
#SBATCH --qos=48h
#SBATCH --array=1-38%6
#SBATCH --cpus-per-task=64
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ugur.cabuk@awi.de
#SBATCH --output=pydamage_bt2_%A_%a.out


#===========================================================================
# VARIABLES
#===========================================================================

WORK=${PWD}

OUTDIR="output"

OUT_FASTP="out.fastp"

OUT_MEGAHIT="out.megahit"

OUT_BOWTIE="out.bowtie"

OUT_PYDAMAGE="out.pydamage"


# Which step to run? 
RUN_BOWTIE="YES"
RUN_PYDAMAGE="YES"


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


# check always new version in the server
module load bowtie2
module load samtools/1.20
module load bamtools/2.5.2


REF=${WORK}/${OUTDIR}/${OUT_MEGAHIT}/${SAMPLE_ID}/final.contigs.fasta

# Bowtie2 index prefix
BT2_INDEX=${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_contigs_index


#===========================================================================
# INPUT READS
#===========================================================================

MERGED_READS=${WORK}/${OUTDIR}/${OUT_FASTP}/${SAMPLE_ID}_fastp_merged_R2.fq.gz

PAIRED_R1=${WORK}/${OUTDIR}/${OUT_FASTP}/${SAMPLE_ID}_fastp_R1.fq.gz

PAIRED_R2=${WORK}/${OUTDIR}/${OUT_FASTP}/${SAMPLE_ID}_fastp_R2.fq.gz


#===========================================================================
# BOWTIE2 MAPPING
#===========================================================================

if [ "${RUN_BOWTIE}" = "YES" ]; then

    echo "============================================================"
    echo "Running Bowtie2"
    echo "Reference: ${REF}"
    echo "============================================================"


    #=======================================================================
    # BUILD BOWTIE2 INDEX
    #=======================================================================

    srun bowtie2-build \
        ${REF} \
        ${BT2_INDEX}


    #=======================================================================
    # MAP MERGED READS
    #
    # --end-to-end:
    #     prevents terminal soft clipping.
    #
    # --very-sensitive:
    #     increases mapping sensitivity.
    #
    # -N 1:
    #     allows one mismatch in the seed.
    #
    #
    #=======================================================================

    srun bowtie2 \
        --end-to-end \
        --very-sensitive \
        -N 1 \
        -p ${CPU} \
        -x ${BT2_INDEX} \
        -U ${MERGED_READS} \
        -S ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_merged.sam


    #=======================================================================
    # MAP PAIRED-END UNMERGED READS
    #=======================================================================

    srun bowtie2 \
        --end-to-end \
        --very-sensitive \
        -N 1 \
        -p ${CPU} \
        -x ${BT2_INDEX} \
        -1 ${PAIRED_R1} \
        -2 ${PAIRED_R2} \
        -S ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_paired.sam


    #=======================================================================
    # SORT SAM -> BAM
    #=======================================================================

    srun samtools sort \
        -@ ${CPU} \
        -o ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_merged.sorted.bam \
        ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_merged.sam


    srun samtools sort \
        -@ ${CPU} \
        -o ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_paired.sorted.bam \
        ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_paired.sam


    #=======================================================================
    # REMOVE TEMPORARY SAM FILES
    #=======================================================================

    rm -f \
        ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_merged.sam

    rm -f \
        ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_paired.sam


    #=======================================================================
    # MERGE:
    # merged reads BAM
    # +
    # paired-end BAM
    #=======================================================================

    srun bamtools merge \
        -in ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_merged.sorted.bam \
        -in ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_paired.sorted.bam \
        -out ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}.merge_paired.bam


    #=======================================================================
    # REMOVE INTERMEDIATE BAM FILES
    #=======================================================================

    rm -f \
        ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_merged.sorted.bam

    rm -f \
        ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}_out_paired.sorted.bam


    #=======================================================================
    # SORT FINAL MERGED BAM
    #=======================================================================

    srun samtools sort \
        -@ ${CPU} \
        -o ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}.merge_paired.sorted.bam \
        ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}.merge_paired.bam


    #=======================================================================
    # INDEX FINAL BAM
    #=======================================================================

    srun samtools index \
        -@ ${CPU} \
        ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}.merge_paired.sorted.bam


    #=======================================================================
    # REMOVE UNSORTED MERGED BAM
    #=======================================================================

    rm -f \
        ${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}.merge_paired.bam


    echo "============================================================"
    echo "Bowtie2 mapping finished."
    echo "============================================================"

else

    echo "Skipping Bowtie2 analysis."

fi


#===========================================================================
# PYDAMAGE
#===========================================================================

if [ "${RUN_PYDAMAGE}" = "YES" ]; then

    module load pydamage/1.0


    BAM=${WORK}/${OUTDIR}/${OUT_BOWTIE}/${SAMPLE_ID}.merge_paired.sorted.bam


    #=======================================================================
    # CHECK BAM INDEX
    #=======================================================================

    if [[ ! -f ${BAM}.bai ]]; then

        echo "BAM index not found. Creating index..."

        srun samtools index \
            -@ ${CPU} \
            ${BAM}

    fi


    #=======================================================================
    # RUN PYDAMAGE
    #=======================================================================

    srun pydamage \
        --outdir ${WORK}/${OUTDIR}/${OUT_PYDAMAGE}/${SAMPLE_ID} \
        analyze \
        ${BAM} \
        -p ${CPU} \
        --force


    #=======================================================================
    # ADD SAMPLE NAME TO pyDamage CSV
    #=======================================================================

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


    module unload pydamage/0.72


    echo "============================================================"
    echo "pyDamage finished."
    echo "============================================================"
    conda deactivate
else

    echo "Skipping pyDamage."

fi
