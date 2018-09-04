#!/bin/bash

WDIR=$1
REFDIR=$2
SPECIES=$3
ANN=$REFDIR/$SPECIES/${SPECIES}.3col

## This would need to be significantly re-written for multi-strain


## make expression table of counts/TPM for featureCounts
cd $WDIR/featureCounts
for i in *fc.tsv
do
  TAG=${i%%.tsv}
  bc_fcount_tpm.pl $i > $TAG.tpm.tsv
done 

echo -e "Gene_id\tSymbol\tGene_type" > $$.names
cat $ANN | sort -k1,1 >> $$.names

FC=`ls *.fc.tpm.tsv`

for i in $FC
do
  TAG=${i%%.fc.tpm.tsv}
  echo $TAG > $TAG.counts.tmp
  echo $TAG > $TAG.TPM.tmp
  awk '{if (NR>2) print}' $i | sort -k1,1 | cut -f 7 >> $TAG.counts.tmp
  awk '{if (NR>2) print}' $i | sort -k1,1 | cut -f 8 >> $TAG.TPM.tmp 
done

paste $$.names *.counts.tmp      > ../exp_tables/featureCounts.counts.tsv  
paste $$.names *.TPM.tmp         > ../exp_tables/featureCounts.TPM.tsv  
rm *counts.tmp *TPM.tmp $$.names

echo "ALL EXPRESSION TABLE PROCESSING IS DONE!" 
