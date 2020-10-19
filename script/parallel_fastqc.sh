#!/bin/bash 

## PIPELINE VERSION

WDIR=$1
CPUS=$2

for i in *.fastq.gz
do
  while [ $(jobs | wc -l) -ge $CPUS ] ; do sleep 5; done
  echo "fastqc: Gathering sequencing metrics for sample $i"
  echo "command: fastqc -q $i"
  fastqc -q $i & 
done 
wait 

rm *zip
mv *html ../FastQC

echo "ALL FASTQC PROCESSING IS DONE!"
