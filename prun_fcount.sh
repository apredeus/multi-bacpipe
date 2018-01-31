#!/bin/bash 

## PIPELINE VERSION

REFDIR=$1
SPECIES=$2
CPUS=$3
STRAND=$4

cd bams

for i in *.bam 
do 
  TAG=${i%%.bam}
  echo "featureCounts: processing sample $TAG, file $i.."
  while [ $(jobs | wc -l) -ge $CPUS ] ; do sleep 5; done
  ../fcount_quant.sh $TAG $REFDIR $SPECIES $STRAND & 
done

wait

rm *summary
cd ../featureCounts
mkdir logs
mv ../bams/*.fc.*log logs 
mv ../bams/*.fc.*tsv . 

echo "ALL FEATURECOUNTS QUANTIFICATION IS DONE!"
echo
echo
