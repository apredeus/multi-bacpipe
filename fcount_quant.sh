#!/bin/bash 

## PIPELINE VERSION

TAG=$1
REFDIR=$2
SPECIES=$3
STRAND=$4
FLAG=""

REF=$REFDIR/Assemblies/$SPECIES.genes.gff
#EXTREF=$REFDIR/Assemblies/$SPECIES.extended.gff

if [[ $TAG == "" || $SPECIES == "" || $REFDIR == "" ]]
then 
  echo "ERROR: Please provide 1) output name (tag) assuming <tag>.fastq.gz or <tag>.R1.fastq.gz/<tag>.R2.fastq.gz input; 2) species/assembly alias, e.g. genprime_v23; 3) strandedness as NONE/FR/RF"
  exit 1
fi
 
if [[ $STRAND == "NONE" ]]
then
  FLAG="0"
  echo "featureCounts: processing sample $TAG using stranded setting $STRAND"
elif [[ $STRAND == "FR" ]]
then
  FLAG="1"
  echo "featureCounts: processing sample $TAG using stranded setting $STRAND"
elif [[ $STRAND == "RF" ]]
then
  FLAG="2"
  echo "featureCounts: processing sample $TAG using stranded setting $STRAND"
else
  echo "ERROR: you must set strand variable to either NONE, FR, or RF"
  exit
fi


#featureCounts -t CDS   -g ID -O -s $FLAG -a $EXTREF -o $TAG.fc.extended.CDS.tsv   $TAG.bam &> $TAG.fc.extended.CDS.log
#featureCounts -t ncRNA -g ID -O -s $FLAG -a $EXTREF -o $TAG.fc.extended.ncRNA.tsv $TAG.bam &> $TAG.fc.extended.ncRNA.log
cp ../strand/$TAG.fc.s$FLAG.tsv $TAG.fc.tsv 
cp ../strand/logs/$TAG.fc.s$FLAG.log $TAG.fc.log
