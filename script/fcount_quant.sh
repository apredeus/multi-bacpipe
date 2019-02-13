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

cd $WDIR/featureCounts

if [[ $STRAND == "NONE" ]]
then
  SFLAG="0"
  echo "featureCounts: processing sample $TAG (strain $SPECIES), strandedness $STRAND (-s $SFLAG), PE options: $PAIRED"
elif [[ $STRAND == "FR" ]]
then
  SFLAG="1"
  echo "featureCounts: processing sample $TAG (strain $SPECIES), strandedness $STRAND (-s $SFLAG), PE options: $PAIRED"
elif [[ $STRAND == "RF" ]]
then
  SFLAG="2"
  echo "featureCounts: processing sample $TAG (strain $SPECIES), strandedness $STRAND (-s $SFLAG), PE options: $PAIRED"
else
  echo "ERROR: you must set strand variable to either NONE, FR, or RF"
  exit
fi

cp ../strand/$TAG.fc.s$SFLAG.tsv $TAG.fc.tsv 
cp ../strand/logs/$TAG.fc.s$SFLAG.log $TAG.fc.log
