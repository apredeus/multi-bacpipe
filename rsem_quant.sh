#!/bin/bash 

## PIPELINE VERSION

TAG=$1
REFDIR=$2
SPECIES=$3
STRAND=$4

READS=""
REF=$REFDIR/RSEM/${SPECIES}_rsem
EXTREF=$REFDIR/RSEM/${SPECIES}.extended_rsem

if [[ $TAG == "" || $REFDIR == "" || $SPECIES == "" || $STRAND == "" ]]
then 
  echo "ERROR: Please provide 1) output name (tag) assuming <tag>.fastq.gz or <tag>.R1.fastq.gz/<tag>.R2.fastq.gz input; 2) species/assembly alias, e.g. genprime_v23; 3) strandedness as NONE/FR/RF"
  exit 1
fi 

if [[ $STRAND == "NONE" ]]
then
  FLAG="--strandedness none"
  echo "Proceeding using stranded setting $STRAND"
elif [[ $STRAND == "FR" ]] 
then
  FLAG="--strandedness forward"
  echo "Proceeding using stranded setting $STRAND"
elif [[ $STRAND == "RF" ]]
then
  FLAG="--strandedness reverse"
  echo "Proceeding using stranded setting $STRAND"
else 
  echo "ERROR: you must set strand variable to either NONE, FR, or RF"
  exit
fi

if [[ -e $TAG.fastq.gz ]]
then
  echo "RSEM: processing quantification as single-end, using bowtie2 index $REF."
  READS="$TAG.fastq.gz"
elif [[ -e $TAG.R1.fastq.gz && -e $TAG.R2.fastq.gz ]]
then
  echo "RSEM: processing quantification as paired-end, using bowtie2 index $REF."
  READS="--paired-end $TAG.R1.fastq.gz $TAG.R2.fastq.gz"
else
  echo "ERROR: The reqiured fastq.gz files were not found!" 
  exit 1
fi

if [[ ! -e $REF.grp || ! -e $EXTREF.grp ]]
then 
  echo "ERROR: rsem index $REF or $EXTREF does not exist!" 
  exit 1
fi

rsem-calculate-expression --bowtie2 --estimate-rspd $READS $EXTREF ${TAG}.expanded.rsem &> $TAG.expanded.rsem.log 
rsem-calculate-expression --bowtie2 --estimate-rspd $READS $REF    ${TAG}.rsem &> $TAG.rsem.log 
