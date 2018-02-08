#!/bin/bash

REFDIR=$1
SPECIES=$2
ANN=$REFDIR/Assemblies/${SPECIES}.3col


## make expression table of counts/TPM for featureCounts
cd featureCounts
for i in *.fc.tsv
do
  TAG=${i%%.fc.tsv}
  ../fcount_tpm.pl $i > $TAG.fc_tpm.tsv
done 

echo -e "Gene_id\tSymbol\tGene_type" > $$.names
cat $ANN | sort -k1,1 >> $$.names

KK1=`ls *.fc_tpm.tsv`

for i in $KK1
do
  TAG=${i%%.fc_tpm.tsv}
  echo $TAG > $TAG.counts.tmp
  echo $TAG > $TAG.TPM.tmp
  awk '{if (NR>2) print}' $i | sort -k1,1 | cut -f 7 >> $TAG.counts.tmp
  awk '{if (NR>2) print}' $i | sort -k1,1 | cut -f 8 >> $TAG.TPM.tmp 
done

paste $$.names *.counts.tmp  > ../exp_tables/featureCounts.counts.tsv  
paste $$.names *.TPM.tmp     > ../exp_tables/featureCounts.TPM.tsv  
rm *.counts.tmp *.TPM.tmp $$.names

## make expression table of counts/TPM for kallisto

cd ../RSEM

echo -e "Gene_id\tSymbol\tGene_type" > $$.names
cat $ANN | sort -k1,1 >> $$.names

KK2=`ls *.genes.results`

for i in $KK2
do
  TAG=${i%%.rsem.genes.results}
  echo $TAG > $TAG.counts.tmp
  echo $TAG > $TAG.TPM.tmp
  awk '{if (NR>1) print}' $i | sort -k1,1 | cut -f 5 >> $TAG.counts.tmp
  awk '{if (NR>1) print}' $i | sort -k1,1 | cut -f 6 >> $TAG.TPM.tmp 
done

paste $$.names *.counts.tmp  > ../exp_tables/rsem.counts.tsv  
paste $$.names *.TPM.tmp     > ../exp_tables/rsem.TPM.tsv  
rm *.counts.tmp *.TPM.tmp $$.names

## make expression table of counts/TPM for RSEM
cd ../kallisto

echo -e "Gene_id\tSymbol\tGene_type" > $$.names
cat $ANN | sort -k1,1 >> $$.names

KK3=`ls *.tsv`

for i in $KK3
do
  TAG=${i%%.tsv}
  echo $TAG > $TAG.counts.tmp
  echo $TAG > $TAG.TPM.tmp
  awk '{if (NR>1) print}' $i | sort -k1,1 | cut -f 4 >> $TAG.counts.tmp
  awk '{if (NR>1) print}' $i | sort -k1,1 | cut -f 5 >> $TAG.TPM.tmp 
done

paste $$.names *.counts.tmp  > ../exp_tables/kallisto.counts.tsv  
paste $$.names *.TPM.tmp     > ../exp_tables/kallisto.TPM.tsv  
rm *.counts.tmp *.TPM.tmp $$.names

echo "ALL EXPRESSION TABLE PROCESSING IS DONE!" 
