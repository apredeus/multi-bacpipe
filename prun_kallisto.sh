#!/bin/bash 

WDIR=$1
REFDIR=$2
SPECIES=$3
CPUS=$4
STRAND=$5
cd $WDIR/cleaned_fastqs 

REF=$REFDIR/kallisto/${SPECIES}_kallisto
if [[ -e $REF ]]
then
  echo "kallisto: using reference $REF"
else 
  echo "ERROR: kallisto reference $REF not found!" 
fi 


KK=`for i in *fastq.gz
do 
  TAG1=${i%%.fastq.gz}
  TAG2=${TAG1%%.R?}
  echo $TAG2
done | sort | uniq`

for i in $KK
do 
  while [ $(jobs | wc -l) -ge $CPUS ] ; do sleep 5; done
  kallisto_quant.sh $i $WDIR $REF $STRAND & 
done
wait

cd ../kallisto
mkdir logs 
mv ../cleaned_fastqs/*log logs 
mv ../cleaned_fastqs/*tsv . 

echo "ALL KALLISTO QUANTIFICATIONS ARE DONE!"
echo
echo
