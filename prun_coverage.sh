#!/bin/bash 

## PIPELINE VERSION

cd bams

REFDIR=$1
SPECIES=$2
CPUS=$3 ## how many parallel jobs to allow

for i in *bam
do
  TAG=${i%%.bam}
  while [ $(jobs | wc -l) -ge $CPUS ] ; do sleep 5; done
  echo "Making bigWig and TDF files for tag $i.." 
  ../calculate_coverage.sh $TAG $REFDIR $SPECIES &> $TAG.coverage.log & 
done
wait

cd ../tdfs_and_bws
mkdir logs
mv ../bams/*tdf .
mv ../bams/*bw .
mv ../bams/*coverage.log logs 

echo "ALL COVERAGE CALCULATIONS ARE DONE!"
echo
echo
