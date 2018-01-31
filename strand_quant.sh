#!/bin/bash 

## PIPELINE VERSION

TAG=$1
REFDIR=$2
SPECIES=$3
REF=$REFDIR/Assemblies/$SPECIES.CDS.gff

if [[ $TAG == "" || $SPECIES == "" || $REFDIR == "" ]]
then 
  echo "ERROR: Please provide 1) output name (tag) assuming <tag>.fastq.gz or <tag>.R1.fastq.gz/<tag>.R2.fastq.gz input; 2) species/assembly alias, e.g. genprime_v23; 3) strandedness as NONE/FR/RF"
  exit 1
fi
 
featureCounts -t CDS -g ID -O -s 0 -a $REF -o $TAG.fc.s0.tsv $TAG.bam &> $TAG.fc.s0.log
featureCounts -t CDS -g ID -O -s 1 -a $REF -o $TAG.fc.s1.tsv $TAG.bam &> $TAG.fc.s1.log
featureCounts -t CDS -g ID -O -s 2 -a $REF -o $TAG.fc.s2.tsv $TAG.bam &> $TAG.fc.s2.log
