#!/bin/bash 

## PIPELINE VERSION

WDIR=$1
REFDIR=$2
SPECIES=$3
CPUS=$4

for i in *.bam 
do 
  TAG=${i%%.bam}
  echo "featureCounts: processing sample $TAG, file $i.."
  while [ $(jobs | wc -l) -ge $CPUS ] ; do sleep 5; done
  strand_quant.sh $TAG $WDIR $REFDIR $SPECIES & 
done

wait

rm *summary
cd ../strand
mkdir logs
mv ../bams/*.fc.*log logs 
mv ../bams/*.fc.*tsv . 

echo "ALL STRANDEDNESS EVALUATION IS DONE!"
