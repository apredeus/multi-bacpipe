#!/bin/bash 

## PIPELINE VERSION

SDIR=$1
WDIR=$2
REFDIR=$3
CONFIG=$4
CPUS=$5
## TODO: adjust number of cores used by calculate coverage scrip

N=`grep -v "^Reference" $CONFIG | wc -l`
KK=`grep -v "^Reference" $CONFIG | cut -f 1`
PP=`grep -v "^Reference" $CONFIG | cut -f 2`
a=( $KK )
b=( $PP )

cd $WDIR/bams 

for i in `seq 0 $((N-1))`
do
  while [ $(jobs | wc -l) -ge $CPUS ] ; do sleep 5; done
  $SDIR/script/calculate_coverage.sh ${a[$i]} $WDIR $REFDIR ${b[$i]} &
done
wait

cd ../tdfs_and_bws
mkdir -p logs
mv ../bams/*tdf .
mv ../bams/*bw .
mv ../bams/*.log logs 

echo "ALL COVERAGE CALCULATIONS ARE DONE!"
