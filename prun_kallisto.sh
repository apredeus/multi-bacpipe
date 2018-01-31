#!/bin/bash 

## PIPELINE VERSION

REFDIR=$1
SPECIES=$2
STRAND=$3
CPUS=$4

ANN=$REFDIR/Assemblies/${SPECIES}/${SPECIES}.4col

cd fastqs 

KK=`for i in *fastq.gz
do 
  TAG1=${i%%.fastq.gz}
  TAG2=${TAG1%%.R?}
  echo $TAG2
done | sort | uniq`

for i in $KK
do 
  echo "kallisto: processing sample $i.."
  ../kallisto_quant.sh $i $REFDIR $SPECIES $STRAND $CPUS
done

cd ../kallisto
mkdir logs 
mv ../fastqs/*log logs 
mv ../fastqs/*tsv . 

echo -e "Transcript_id\tGene_id\tSymbol\tGene_type" > names
cat $ANN | grep -v -i symbol >> names

for i in *.tsv
do
  TAG=${i%%.tsv}
  echo $TAG > $TAG.tmp
  ## discard the header, use counts 
  awk '{if (NR>1) print $1,$4}' $i | sort -k1,1 | awk '{print $2}' >> $TAG.tmp
done

paste names *.tmp > Kallisto_transcript_counts.tsv
rm *.tmp names 


echo "ALL KALLISTO QUANTIFICATIONS ARE DONE!"
echo
echo
