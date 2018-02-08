#!/bin/bash 

REFDIR=$1
SPECIES=$2
CPUS=$3
STRAND=$4

REF=$REFDIR/kallisto/${SPECIES}_kallisto
if [[ -e $REF ]]
then
  echo "kallisto: using reference $REF"
else 
  echo "ERROR: kallisto reference $REF not found!" 
fi 

cd cleaned_fastqs 

KK=`for i in *fastq.gz
do 
  TAG1=${i%%.fastq.gz}
  TAG2=${TAG1%%.R?}
  echo $TAG2
done | sort | uniq`

for i in $KK
do 
  while [ $(jobs | wc -l) -ge $CPUS ] ; do sleep 5; done
  ../kallisto_quant.sh $i $REF $STRAND & 
done
wait

cd ../kallisto
mkdir logs 
mv ../cleaned_fastqs/*log logs 
mv ../cleaned_fastqs/*tsv . 

echo "ALL KALLISTO QUANTIFICATIONS ARE DONE!"
echo
echo
