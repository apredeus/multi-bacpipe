#!/bin/bash 

## PIPELINE VERSION

WDIR=$1

cd $WDIR/fastqs 
KK=`for i in *fastq.gz
do 
  TAG1=${i%%.fastq.gz}
  TAG2=${TAG1%%.R?}
  echo $TAG2
done | sort | uniq`

cd ../stats

for i in $KK
do
  calculate_stats.sh $i $WDIR &  
done
wait

echo 
echo "Read number statistics for the processed dataset:"
echo -e "Sample\tAll\trRNA\ttRNA\tunmapped\tmultimap\t1_loc\tassigned" > $$.counts.tmp
echo -e "===================\t=======\t=======\t=======\t=======\t=======\t=======\t=======" >> $$.counts.tmp
cat *.counts >> $$.counts.tmp
cat $$.counts.tmp | column -t 
rm $$.counts.tmp 

echo 
echo "Read percentage statistics for the processed dataset:"
echo -e "Sample\trRNA\ttRNA\tunmapped\tmultimap\t1_loc\tassigned" > $$.percent.tmp
echo -e "===================\t=======\t=======\t=======\t=======\t=======\t=======" >> $$.percent.tmp
cat *.percent >> $$.percent.tmp
cat $$.percent.tmp | column -t 
rm $$.percent.tmp 

echo
echo "ALL RNA-STAT CALCULATIONS ARE DONE!"
