#!/bin/bash 

## PIPELINE VERSION

CPUS=$1

cd fastqs

for i in *.fastq.gz
do
  while [ $(jobs | wc -l) -ge $CPUS ] ; do sleep 5; done
  echo "fastqc: Gathering sequencing metrics for sample $i"
  fastqc -q $i & 
done 
wait 

rm *zip
mv *html ../FastQC

echo "ALL FASTQC PROCESSING IS DONE!"
