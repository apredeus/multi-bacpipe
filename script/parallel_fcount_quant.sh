#!/bin/bash 

## PIPELINE VERSION

SDIR=$1
WDIR=$2
REFDIR=$3
CONFIG=$4
CPUS=$5
STRAND=$6


N=`grep -v "^Reference" $CONFIG | wc -l`
KK=`grep -v "^Reference" $CONFIG | cut -f 1`
PP=`grep -v "^Reference" $CONFIG | cut -f 2`
a=( $KK ) 
b=( $PP ) 

for i in `seq 0 $((N-1))`
do
  while [ $(jobs | wc -l) -ge $CPUS ] ; do sleep 5; done
  $SDIR/script/fcount_quant.sh ${a[$i]} $WDIR $REFDIR ${b[$i]} $STRAND &
done
wait

cd ../featureCounts
mkdir -p logs
mv *log logs 

echo "ALL FEATURECOUNTS QUANTIFICATION IS DONE!"
