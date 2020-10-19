#!/bin/bash 

## PIPELINE VERSION

SDIR=$1
WDIR=$2
REFDIR=$3
CONFIG=$4
CPUS=$5
NJOB=""  ## number of samples to be processed at the same time 
NRUN=""  ## number of cores used per process; not more than 8.

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


cd $WDIR/bams 

for i in `seq 0 $((N-1))`
do
  while [ $(jobs | wc -l) -ge $NJOB ] ; do sleep 5; done
  $SDIR/script/calculate_coverage.sh ${a[$i]} $WDIR $REFDIR ${b[$i]} $NRUN &
done
wait

cd ../tdfs_and_bws
mkdir -p logs
mv ../bams/*tdf .
mv ../bams/*bw .
mv ../bams/*.log logs 

echo "ALL COVERAGE CALCULATIONS ARE DONE!"
