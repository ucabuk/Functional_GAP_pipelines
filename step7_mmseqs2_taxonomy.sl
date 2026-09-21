#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# Contact: ugur.cabuk@awi.de
# Desc: Taxonomic assignment of the non-redundant protein catalog with MMseqs2
#===========================================================================

#SBATCH --account=envi.envi
#SBATCH --job-name=mmseqs2
#SBATCH --partition=fat
#SBATCH --time=48:00:00
#SBATCH --qos=48h
#SBATCH --mem=600G
#SBATCH --cpus-per-task=64
#SBATCH --output=mmseqs2_%j.out

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

OUT_MMSEQS="out.mmseqs2"


INPUT="non_redundant_prot_catalog.faa"


CPU=${SLURM_CPUS_PER_TASK}


# Change this path according to the MMseqs2 NR database on your system.
NR_DB="/albedo/work/projects/p_biodiv_shotgun/05_GAP_data/nr_mmseq2_db/NR_ncbi"


LCA_RANK="species,genus,family,order,class,phylum,kingdom,superkingdom"


mkdir -p ${OUTDIR}/${OUT_MMSEQS}


#===========================================================================
# MMSEQS2 DATABASE
#===========================================================================

module load mmseqs2/15.6f452


srun mmseqs createdb \
    ${OUTDIR}/${OUT_PRODIGAL}/${INPUT} \
    ${OUTDIR}/${OUT_MMSEQS}/${INPUT}.db2


#===========================================================================
# TAXONOMIC ASSIGNMENT
#===========================================================================

srun mmseqs taxonomy \
    --split-memory-limit 500G \
    --threads ${CPU} \
    -e 0.00001 \
    --tax-lineage 1 \
    --lca-ranks ${LCA_RANK} \
    ${OUTDIR}/${OUT_MMSEQS}/${INPUT}.db2 \
    ${NR_DB} \
    ${OUTDIR}/${OUT_MMSEQS}/${INPUT}.result \
    ${OUTDIR}/${OUT_MMSEQS}/tmp_${INPUT} \
    --lca-mode 3 \
    -s 5


#===========================================================================
# CREATE TSV OUTPUT
#===========================================================================

srun mmseqs createtsv \
    ${OUTDIR}/${OUT_MMSEQS}/${INPUT}.db2 \
    ${NR_DB} \
    ${OUTDIR}/${OUT_MMSEQS}/${INPUT}.result \
    ${OUTDIR}/${OUT_MMSEQS}/${INPUT}.result.tsv


module unload mmseqs2/15.6f452


echo "============================================================"
echo "MMseqs2 taxonomic assignment finished."
echo "============================================================"
