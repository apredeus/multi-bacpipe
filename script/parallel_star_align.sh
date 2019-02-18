#!/bin/bash 

## PIPELINE VERSION

SDIR=$1
WDIR=$2
REFDIR=$3
CONFIG=$4
CPUS=$5
## TODO: adjustable number of cores dep on RAM/CPUs
NJOB=$((CPUS/4))

N=`grep -v "^Reference" $CONFIG | wc -l`
KK=`grep -v "^Reference" $CONFIG | cut -f 1`
PP=`grep -v "^Reference" $CONFIG | cut -f 2`
a=( $KK ) 
b=( $PP ) 

for i in `seq 0 $((N-1))`
do
  while [ $(jobs | wc -l) -ge $NJOB ] ; do sleep 5; done
  $SDIR/script/star_align.sh ${a[$i]} $WDIR $REFDIR ${b[$i]} 4 &
done
wait

cd ../bams 
mkdir -p logs
mv ../fastqs/*_star/*bam .
mv ../fastqs/*_star/*bam.bai . 
mv ../fastqs/*_star/*.star.log logs 

rm -rf ../fastqs/*_star/

echo "ALL STAR ALIGMENTS ARE DONE!"
