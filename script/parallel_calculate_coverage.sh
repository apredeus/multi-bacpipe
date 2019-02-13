#!/bin/bash 

## PIPELINE VERSION

WDIR=$1
REFDIR=$2
CONFIG=$3
CPUS=$4

N=`grep -v "^Reference" $CONFIG | wc -l`
KK=`grep -v "^Reference" $CONFIG | cut -f 1`
PP=`grep -v "^Reference" $CONFIG | cut -f 2`
a=( $KK )
b=( $PP )

cd $WDIR/bams 

for i in `seq 0 $((N-1))`
do
  while [ $(jobs | wc -l) -ge $CPUS ] ; do sleep 5; done
  calculate_coverage.sh ${a[$i]} $WDIR $REFDIR ${b[$i]} &
done
wait

cd ../tdfs_and_bws
mkdir logs
mv ../bams/*tdf .
mv ../bams/*bw .
mv ../bams/*.log logs 

echo "ALL COVERAGE CALCULATIONS ARE DONE!"
