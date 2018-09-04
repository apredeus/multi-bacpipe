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

for i in `seq 0 $((N-1))`
do
  while [ $(jobs | wc -l) -ge $CPUS ] ; do sleep 5; done
  mbc_strand_quant.sh ${a[$i]} $WDIR $REFDIR ${b[$i]} &
done
wait

rm *summary
cd ../strand
mkdir logs
mv ../bams/*.fc.*log logs 
mv ../bams/*.fc.*tsv . 

echo "ALL STRANDEDNESS EVALUATION IS DONE!"
