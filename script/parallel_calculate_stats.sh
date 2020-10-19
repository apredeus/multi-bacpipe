#!/bin/bash 

SDIR=$1
WDIR=$2
CONFIG=$3
CPUS=$4

N=`grep -v "^Reference" $CONFIG | wc -l`
KK=`grep -v "^Reference" $CONFIG | cut -f 1`
PP=`grep -v "^Reference" $CONFIG | cut -f 2`
a=( $KK ) 
b=( $PP ) 

cd $WDIR/fastqs

for i in `seq 0 $((N-1))`
do
  while [ $(jobs | wc -l) -ge $CPUS ] ; do sleep 5; done
  $SDIR/script/calculate_stats.sh ${a[$i]} $WDIR ${b[$i]} &  
done
wait

cd $WDIR/stats 
echo 
echo "Read number statistics for the processed dataset:"
echo
echo -e "Sample\tStrain\tAll\trRNA\tunmapped\tmultimap\t1_loc\tassigned" > $$.counts.tmp
echo -e "===================\t========\t=======\t=======\t=======\t=======\t=======\t=======" >> $$.counts.tmp
cat *.counts >> $$.counts.tmp
cat $$.counts.tmp | column -t 
rm $$.counts.tmp 

echo 
echo "Read percentage statistics for the processed dataset:"
echo
echo -e "Sample\tStrain\trRNA\tunmapped\tmultimap\t1_loc\tassigned" > $$.percent.tmp
echo -e "===================\t========\t=======\t=======\t=======\t=======\t=======" >> $$.percent.tmp
cat *.percent >> $$.percent.tmp
cat $$.percent.tmp | column -t 
rm $$.percent.tmp 

echo
echo "ALL RNA-STAT CALCULATIONS ARE DONE!"
