#!/bin/bash

WDIR=$1
CONFIG=$2
REFSTR=$3

N1=`grep -v "^Reference" $CONFIG | wc -l`
NR=`grep    "^Reference" $CONFIG | wc -l`
NS=`grep -v "^Reference" $CONFIG | cut -f 2 | sort | uniq | wc -l` 
KK=`grep -v "^Reference" $CONFIG | cut -f 2 | sort | uniq`

cd $WDIR/featureCounts

## make expression table of counts/TPM for featureCounts
for i in *fc.tsv
do
  TAG=${i%%.tsv}
  mbc_fcount_tpm.pl $i > $TAG.tpm.tsv
done 

FC=`ls *.fc.tpm.tsv`
for i in $FC
do
  TAG=${i%%.fc.tpm.tsv}
  echo $TAG > $TAG.counts.tmp
  echo $TAG > $TAG.TPM.tmp
  awk '{if (NR>2) print}' $i | sort -k1,1 | cut -f 7 >> $TAG.counts.tmp
  awk '{if (NR>2) print}' $i | sort -k1,1 | cut -f 8 >> $TAG.TPM.tmp 
done

for i in $KK 
do
  echo "Making raw count and TPM tables for strain $i."
  PP=`grep -P "\t$i" $CONFIG | cut -f 1`
  FIRST=`echo $PP | awk '{print $1}'`
  echo -e $i > $i.names
  awk '{if (NR>2) print $1}' $FIRST.fc.tpm.tsv | sort -k1,1 >> $i.names
  CNT_LIST=`echo $PP | perl -ne '@t=split/\s+/; foreach $i (@t) {print "$i.counts.tmp "}; print "\n"'`
  TPM_LIST=`echo $PP | perl -ne '@t=split/\s+/; foreach $i (@t) {print "$i.TPM.tmp "}; print "\n"'`
  paste $i.names $CNT_LIST > ../exp_tables/$i.counts.tsv
  paste $i.names $TPM_LIST > ../exp_tables/$i.TPM.tsv
done

rm *.names *.counts.tmp *TPM.tmp

cd $WDIR/exp_tables

## run a final script 
## both CDS_orthologs.tsv and ncRNA_orthologs.tsv should be annotated with chr/plasmid/prophage

unify_exp_tables.pl $WDIR $CONFIG $WDIR/refstr/annotated_CDS.tsv $WDIR/refstr/annotated_ncRNA.tsv $REFSTR

echo "ALL EXPRESSION TABLE PROCESSING IS DONE!" 
