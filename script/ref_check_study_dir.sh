#!/bin/bash 

WDIR=$1
TAG=$2

cd $WDIR/study_strains/$TAG

echo "==> Checking study strain $TAG:"

if [[ -s $TAG.gene.gff && -s $TAG.chrom.sizes && -s $TAG.genome.fa && -s $TAG.genome.fa.fai && \
-s $TAG.rRNA.bed && -s $TAG.prophage.bed && "$(ls -A ${TAG}.STAR)" ]]
then
  echo "- all necessary annotation files and directories for $TAG exist and are not empty."
else 
  >&2 echo "ERROR: some of the necessary annotation files for $TAG do not exist or are empty!"
  exit 1
fi
