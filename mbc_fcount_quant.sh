#!/bin/bash 

## PIPELINE VERSION

TAG=$1
WDIR=$2
REFDIR=$3
SPECIES=$4
STRAND=$5

MFLAG=""
SFLAG=""
PAIRED="" 

cd $WDIR/fastqs 

if [[ -e $TAG.R1.fastq.gz && -e $TAG.R2.fastq.gz ]]
then
  PAIRED="-p" 
fi

MFLAG="-M --fraction" 

cd $WDIR/featureCounts

if [[ $TAG == "" || $SPECIES == "" || $REFDIR == "" ]]
then 
  echo "ERROR: Please provide 1) output name (tag) assuming <tag>.fastq.gz or <tag>.R1.fastq.gz/<tag>.R2.fastq.gz input; 2) species/assembly alias, e.g. genprime_v23; 3) strandedness as NONE/FR/RF"
  exit 1
fi
 
if [[ $STRAND == "NONE" ]]
then
  SFLAG="0"
  echo "featureCounts: processing sample $TAG, strandedness $STRAND (-s $SFLAG), PE options: $PAIRED"
elif [[ $STRAND == "FR" ]]
then
  SFLAG="1"
  echo "featureCounts: processing sample $TAG, strandedness $STRAND (-s $SFLAG), PE options: $PAIRED"
elif [[ $STRAND == "RF" ]]
then
  SFLAG="2"
  echo "featureCounts: processing sample $TAG, strandedness $STRAND (-s $SFLAG), PE options: $PAIRED"
else
  echo "ERROR: you must set strand variable to either NONE, FR, or RF"
  exit
fi

cp ../strand/$TAG.fc.s$SFLAG.tsv $TAG.fc.tsv 
cp ../strand/logs/$TAG.fc.s$SFLAG.log $TAG.fc.log
