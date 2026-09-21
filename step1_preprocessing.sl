#!/bin/bash

#===========================================================================
# slurm batch script to run the shotgun pipeline step 1
# on several sample using arrays
# Version 0.4
#
# by Lars Harms
#
# contact: lars.harms@awi.de
#
# slurm options and variables under >set required variables<
# have to be modified by the user
#=============================================================================


INDIR="path-to-raw-sequencing-data"

R1_ENDING="_R1.fastq.gz"
R2_ENDING="_R2.fastq.gz"
DEDUP="TRUE"
FILTER="--low_complexity_filter"

# given variables (please do not change)
#===================================================================
WORK=${PWD}
OUTDIR="output"
OUT_FQC="out.fastqc"
PRE="pre"
POST="post"
OUT_DEDUP="out.dedup"
OUT_FASTP="out.fastp"

CPU=${SLURM_CPUS_PER_TASK}

# prepare environment
#===================================================================
mkdir -p ${OUTDIR}/${OUT_FQC}/${PRE}
mkdir -p ${OUTDIR}/${OUT_FQC}/${POST}
mkdir -p ${OUTDIR}/${OUT_DEDUP}
mkdir -p ${OUTDIR}/${OUT_FASTP}


cd ${INDIR}
FILE_R1=$(ls *${R1_ENDING} | sed -n ${SLURM_ARRAY_TASK_ID}p)
FILE_R2=$(ls *${R2_ENDING} | sed -n ${SLURM_ARRAY_TASK_ID}p)

FILEBASE=${FILE_R1%${R1_ENDING}}

OUT_R1_CL="${FILEBASE}_dedup_R1.fq.gz"
OUT_R2_CL="${FILEBASE}_dedup_R2.fq.gz"

OUT_R1="${FILEBASE}_fastp_R1.fq.gz"
OUT_R2="${FILEBASE}_fastp_R2.fq.gz"

OUT_MERGED="${FILEBASE}_fastp_merged_R2.fq.gz"

cd ${WORK}

# tasks to be performed
#===================================================================

# FASTQC PRE
#----------

fastqc -q -o ${OUTDIR}/${OUT_FQC}/${PRE} -t 2 ${INDIR}/${FILE_R1} ${INDIR}/${FILE_R2}

# CLUMPIFY
#----------
if [ ${DEDUP} == "FALSE" ]; then
        echo "Removal of read duplications is turned off."
else
        srun clumpify.sh in=${INDIR}/${FILE_R1} in2=${INDIR}/${FILE_R2} out=${OUTDIR}/${OUT_DEDUP}/${OUT_R1_CL} out2=${OUTDIR}/${OUT_DEDUP}/${OUT_R2_CL} dedupe=t
fi

# FASTP
#----------
if [ ${DEDUP} == "FALSE" ]; then       
        fastp --in1 ${INDIR}/${FILE_R1} --in2 ${INDIR}/${FILE_R2} --out1 ${OUTDIR}/${OUT_FASTP}/${OUT_R1} --out2 ${OUTDIR}/${OUT_FASTP}/${OUT_R2} -m --merged_out ${OUTDIR}/${OUT_FASTP}/${OUT_MERGED} ${FILTER} -w ${CPU} --verbose --json=${OUTDIR}/${OUT_FASTP}/${FILEBASE}.json --html=${OUTDIR}/${OUT_FASTP}/${FILEBASE}.html
else
        srun fastp --in1 ${OUTDIR}/${OUT_DEDUP}/${OUT_R1_CL} --in2 ${OUTDIR}/${OUT_DEDUP}/${OUT_R2_CL} --out1 ${OUTDIR}/${OUT_FASTP}/${OUT_R1} --out2 ${OUTDIR}/${OUT_FASTP}/${OUT_R2} -m --merged_out ${OUTDIR}/${OUT_FASTP}/${OUT_MERGED} ${FILTER} -w ${CPU} --verbose --json=${OUTDIR}/${OUT_FASTP}/${FILEBASE}.json --html=${OUTDIR}/${OUT_FASTP}/${FILEBASE}.html
fi

# FASTQC POST
#----------
srun fastqc -q -o ${OUTDIR}/${OUT_FQC}/${POST} -t 3 ${OUTDIR}/${OUT_FASTP}/${OUT_R1} ${OUTDIR}/${OUT_FASTP}/${OUT_R2} ${OUTDIR}/${OUT_FASTP}/${OUT_MERGED}
