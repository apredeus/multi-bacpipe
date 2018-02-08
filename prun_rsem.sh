#!/bin/bash 

REFDIR=$1
SPECIES=$2
CPUS=$3
STRAND=$4

REF=$REFDIR/RSEM/${SPECIES}_rsem
if [[ -e $REF.ti ]]
then
  echo "RSEM: using reference $REF"
else 
  echo "ERROR: rsem reference $REF not found!" 
fi 

source activate rsem
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
  ../rsem_quant.sh $i $REF $STRAND & 
done
wait

rm -rf *rsem.transcript.bam *rsem.stat *rsem.isoforms.results
cd ../RSEM
mkdir logs 
mv ../cleaned_fastqs/*log logs 
mv ../cleaned_fastqs/*results . 

echo "ALL RSEM QUANTIFICATION IS DONE!"
