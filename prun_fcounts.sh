#!/bin/bash 

## PIPELINE VERSION

REFDIR=$1
SPECIES=$2
NJOB=$3

cd bams

for i in *.bam 
do 
  TAG=${i%%.bam}
  echo "featureCounts: processing sample $TAG, file $i.."
  while [ $(jobs | wc -l) -ge $NJOB ] ; do sleep 5; done
  ../fcount_quant.sh $TAG $REFDIR $SPECIES & 
done

wait 

echo "ALL FEATURECOUNTS QUANTIFICATION IS DONE!"
echo
echo
