#!/bin/bash 

## PIPELINE VERSION

cd fastqs

for i in *.fastq.gz
do
  echo "fastqc: Gathering sequencing metrics for sample $i"
  fastqc -q $i & 
done 
wait 

rm *zip
mv *html ../FastQC

echo "ALL FASTQC PROCESSING IS DONE!"
echo
echo

