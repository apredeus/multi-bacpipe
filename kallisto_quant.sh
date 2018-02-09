#!/bin/bash 

## PIPELINE VERSION

TAG=$1
WDIR=$2
REF=$3
STRAND=$4
cd $WDIR/cleaned_fastqs 

SINGLE=""
FLAG=""
READS=""
WDIR=`pwd`

if [[ $TAG == "" || $REF == "" || $STRAND == "" ]]
then 
  echo "ERROR: Please provide 1) output name (tag) assuming <tag>.fastq.gz or <tag>.R1.fastq.gz/<tag>.R2.fastq.gz input; 2) species/assembly alias, e.g. genprime_v23"
  exit 1
fi 

if [[ $STRAND == "RF" ]]
then
  FLAG="--rf-stranded"
elif [[ $STRAND == "FR" ]]
then
  FLAG="--fr-stranded"
fi 

if [[ -e $TAG.fastq.gz ]]
then 
  echo "kallisto: processing alignment of $TAG as single-end (l 200, s 50)."
  READS="$TAG.fastq.gz"
  SINGLE="--single -l 200 -s 50"
elif [[ -e $TAG.R1.fastq.gz && -e $TAG.R2.fastq.gz ]]
then
  echo "kallisto: processing alignment of $TAG as paired-end."
  READS="$TAG.R1.fastq.gz $TAG.R2.fastq.gz"
else
  echo "ERROR: The reqiured fastq.gz files were not found!" 
  exit 1
fi

kallisto quant -i $REF $FLAG $SINGLE --plaintext -o ${TAG}_kallisto $READS &> $TAG.kallisto.log
mv ${TAG}_kallisto/abundance.tsv $TAG.tsv
rm -rf ${TAG}_kallisto
 
