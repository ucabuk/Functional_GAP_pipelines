#!/bin/bash
#===========================================================================
# Author: Ugur Cabuk
# contact: ugur.cabuk@awi.de
# Desc: eggNOG Functional Annotation
#===========================================================================

WORK=${PWD}
OUTDIR="output"
OUT_PRODIGAL="out.prodigal"
OUT_METAEUK="out.metaeuk"
OUT_TIARA="out.tiara"
OUT_EGGNOG="out.eggnog"

# For prokGAP example.
IN="non_redundant_PROKGAP_protein.faa"
OUT="non_redundant_PROKGAP_protein"

EGGNOG_DB_DIR="/path-to-database/eggnog"
export EGGNOG_DATA_DIR=${EGGNOG_DB_DIR}

mkdir -p tmp

EGGNOG_MODE="diamond"

emapper.py -m ${EGGNOG_MODE} -i ${WORK}/${OUTDIR}/${OUT_PRODIGAL}/${IN} -o ${WORK}/${OUTDIR}/${OUT_EGGNOG}/${OUT} --temp_dir tmp
