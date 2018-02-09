#!/bin/bash 

## PIPELINE VERSION

WDIR=$1
REFDIR=$2
SPECIES=$3
CPUS=$4
STRAND=$5

cd $WDIR/bams

for i in *.bam 
do 
  TAG=${i%%.bam}
  while [ $(jobs | wc -l) -ge $CPUS ] ; do sleep 5; done
  fcount_quant.sh $TAG $WDIR $REFDIR $SPECIES $STRAND & 
done

wait

cd ../featureCounts
mkdir logs
mv ../bams/*.fc.*log logs 
mv ../bams/*.fc.*tsv . 

echo "ALL FEATURECOUNTS QUANTIFICATION IS DONE!"
