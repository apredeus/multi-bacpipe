#!/bin/bash 

## PIPELINE VERSION

WDIR=$1
REFDIR=$2
CONFIG=$3
CPUS=$4
STRAND=$5


N=`grep -v "^Reference" $CONFIG | wc -l`
KK=`grep -v "^Reference" $CONFIG | cut -f 1`
PP=`grep -v "^Reference" $CONFIG | cut -f 2`
a=( $KK ) 
b=( $PP ) 

for i in `seq 0 $((N-1))`
do
  while [ $(jobs | wc -l) -ge $CPUS ] ; do sleep 5; done
  mbc_fcount_quant.sh ${a[$i]} $WDIR $REFDIR ${b[$i]} $STRAND &
done
wait

rm *summary
cd ../featureCounts
mkdir logs
mv ../bams/*fc.*log logs 

echo "ALL FEATURECOUNTS QUANTIFICATION IS DONE!"
