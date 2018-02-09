#!/bin/bash 

## PIPELINE VERSION

WDIR=$1
REFDIR=$2
SPECIES=$3
CPUS=$4

cd $WDIR/bams

for i in *bam
do
  TAG=${i%%.bam}
  while [ $(jobs | wc -l) -ge $CPUS ] ; do sleep 5; done
  calculate_coverage.sh $TAG $WDIR $REFDIR $SPECIES & 
done
wait

cd ../tdfs_and_bws
mkdir logs
mv ../bams/*tdf .
mv ../bams/*bw .
mv ../bams/*.log logs 

echo "ALL COVERAGE CALCULATIONS ARE DONE!"
