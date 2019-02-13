#!/bin/bash 

## PIPELINE VERSION

WDIR=$1
REFDIR=$2
CONFIG=$3
CPUS=$4
NJOB=$((CPUS/4))

N=`grep -v "^Reference" $CONFIG | wc -l`
KK=`grep -v "^Reference" $CONFIG | cut -f 1`
PP=`grep -v "^Reference" $CONFIG | cut -f 2`
a=( $KK ) 
b=( $PP ) 

for i in `seq 0 $((N-1))`
do
  while [ $(jobs | wc -l) -ge $NJOB ] ; do sleep 5; done
  star_align.sh ${a[$i]} $WDIR $REFDIR ${b[$i]} 4 &
done
wait

cd ../bams 
mkdir logs
mv ../fastqs/*_star/*bam .
mv ../fastqs/*_star/*bam.bai . 
mv ../fastqs/*_star/*.star.log logs 

rm -rf ../fastqs/*_star/

echo "ALL STAR ALIGMENTS ARE DONE!"
