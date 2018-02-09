#!/bin/bash 

## PIPELINE VERSION

TAG=$1
WDIR=$2
REFDIR=$3
SPECIES=$4
REF=$REFDIR/Assemblies/$SPECIES.gene.gff
cd $WDIR/bams 

if [[ $TAG == "" || $SPECIES == "" || $REFDIR == "" ]]
then 
  echo "ERROR: Please provide 1) output name (tag); 2) reference directory (see manual); 3) species tag (see manual)"
  exit 1
fi
 
featureCounts -t gene -g ID -O -s 0 -a $REF -o $TAG.fc.s0.tsv $TAG.bam &> $TAG.fc.s0.log
featureCounts -t gene -g ID -O -s 1 -a $REF -o $TAG.fc.s1.tsv $TAG.bam &> $TAG.fc.s1.log
featureCounts -t gene -g ID -O -s 2 -a $REF -o $TAG.fc.s2.tsv $TAG.bam &> $TAG.fc.s2.log
