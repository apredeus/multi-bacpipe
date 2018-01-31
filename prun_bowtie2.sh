#!/bin/bash 

## PIPELINE VERSION

cd fastqs

REFDIR=$1
SPECIES=$2
CPUS=$3 ## how many cores for individual jobs
NJOB=$((CPUS/4))

KK=`for i in *fastq.gz
do 
  TAG1=${i%%.fastq.gz}
  TAG2=${TAG1%%.R?}
  echo $TAG2
done | sort | uniq`

for i in $KK
do
  while [ $(jobs | wc -l) -ge $NJOB ] ; do sleep 5; done
  ../bowtie2_align.sh $i $REFDIR $SPECIES 4 &
done
wait

mv *_rrna.log *_trna.log ../rtRNA
cd ../bams 
mkdir logs 
mv ../fastqs/*bam .
mv ../fastqs/*bam.bai . 
mv ../fastqs/*log logs 

echo "ALL BOWTIE2 ALIGMENT IS DONE!"
echo
echo
