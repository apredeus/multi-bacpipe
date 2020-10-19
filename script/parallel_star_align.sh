#!/bin/bash 

## PIPELINE VERSION

SDIR=$1
WDIR=$2
REFDIR=$3
CONFIG=$4
CPUS=$5
NJOB=""  ## number of samples to be aligned at the same time
NRUN=""  ## number of cores used per alignment; not more than 8.

N=`grep -v "^Reference" $CONFIG | wc -l`
KK=`grep -v "^Reference" $CONFIG | cut -f 1`
PP=`grep -v "^Reference" $CONFIG | cut -f 2`
a=( $KK ) 
b=( $PP ) 

if (( $CPUS <= $N ))
then 
  NRUN=1
  NJOB=$CPUS
elif (( $CPUS > $N && $CPUS < $N*4)) 
then
  NRUN=2
  NJOB=$((CPUS/2))
elif (( $CPUS >= $N*4 && $CPUS < $N*8)) 
then 
  NRUN=4
  NJOB=$((CPUS/4))
elif (( $CPUS >= $N*8 ))
then
  NRUN=8
  NJOB=$((CPUS/8))
fi

for i in `seq 0 $((N-1))`
do
  while [ $(jobs | wc -l) -ge $NJOB ] ; do sleep 5; done
  $SDIR/script/star_align.sh ${a[$i]} $WDIR $REFDIR ${b[$i]} $NRUN &
done
wait

cd ../bams 
mkdir -p logs
mv ../fastqs/*_star/*bam .
mv ../fastqs/*_star/*bam.bai . 
mv ../fastqs/*_star/*.star.log logs 

rm -rf ../fastqs/*_star/

echo "ALL STAR ALIGMENTS ARE DONE!"
