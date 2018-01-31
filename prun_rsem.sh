#!/bin/bash 

## PIPELINE VERSION

REFDIR=$1
SPECIES=$2
CPUS=$3
STRAND=$4

#ANN=$REFDIR/Assemblies/${SPECIES}/${SPECIES}.3col

source activate rsem

cd fastqs 
KK=`for i in *fastq.gz
do 
  TAG1=${i%%.fastq.gz}
  TAG2=${TAG1%%.R?}
  echo $TAG2
done | sort | uniq`

for i in $KK
do 
  echo "RSEM: processing sample $i.."
  while [ $(jobs | wc -l) -ge $CPUS ] ; do sleep 5; done
  ../rsem_quant.sh $i $REFDIR $SPECIES $STRAND & 
done

wait

cd ../RSEM
mkdir logs 
mv ../fastqs/*log logs 
mv ../fastqs/*results . 

echo "ALL RSEM QUANTIFICATION IS DONE!"
echo
echo
