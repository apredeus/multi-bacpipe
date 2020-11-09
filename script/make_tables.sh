#!/bin/bash

SDIR=$1
WDIR=$2
CONFIG=$3
KEYSTR=$4

REFNUM=`grep    "^Reference" $CONFIG | wc -l`
STUDY=`grep -v "^Reference" $CONFIG | cut -f 2 | sort | uniq`

cd $WDIR/featureCounts

if [[ $REFNUM != 0 ]]
then
  ## make expression table of counts/TPM for featureCounts
  ## multi-strain workflow
  for i in *fc.tsv
  do
    SAMPLE=${i%%.tsv}
    $SDIR/script/fcount_tpm.pl $i > $SAMPLE.tpm.tsv
  done 
  
  FC=`ls *.fc.tpm.tsv`
  for i in $FC
  do
    SAMPLE=${i%%.fc.tpm.tsv}
    echo $SAMPLE > $SAMPLE.counts.tmp
    echo $SAMPLE > $SAMPLE.TPM.tmp
    awk '{if (NR>2) print}' $i | sort -k1,1 | cut -f 7 >> $SAMPLE.counts.tmp
    awk '{if (NR>2) print}' $i | sort -k1,1 | cut -f 8 >> $SAMPLE.TPM.tmp 
  done
  
  for i in $STUDY 
  do
    echo "Making raw count and TPM tables for strain $i."
    PP=`grep -P "\t$i$" $CONFIG | cut -f 1`
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
  ## TODO: this would need to be rewritten for a new unified processing   
 
  $SDIR/script/unify_exp_tables.pl $WDIR $CONFIG $WDIR/orthologs.tsv $KEYSTR
  
  echo "ALL EXPRESSION TABLE PROCESSING IS DONE!"
else
  ## single-strain workflow, a lot simpler processing of the final tables 
  ## just make tables, annotate with name/type/prophage (latter is optional)  
  TAG=`cut -f 2 $CONFIG | uniq`
 
  for i in *fc.tsv
  do
    SAMPLE=${i%%.tsv}
    $SDIR/script/fcount_tpm.pl $i > $SAMPLE.tpm.tsv
  done 
  
  for i in *.fc.tpm.tsv
  do
    SAMPLE=${i%%.fc.tpm.tsv}
    echo $SAMPLE > $SAMPLE.counts.tmp
    echo $SAMPLE > $SAMPLE.TPM.tmp
    awk '{if (NR>2) print}' $i | sort -k1,1 | cut -f 7 >> $SAMPLE.counts.tmp
    awk '{if (NR>2) print}' $i | sort -k1,1 | cut -f 8 >> $SAMPLE.TPM.tmp 
  done
  
  echo "Making raw count and TPM tables for strain $TAG."
  PP=`grep -P "\t$TAG$" $CONFIG | cut -f 1`
  FIRST=`echo $PP | awk '{print $1}'`
  echo -e $TAG > $TAG.names
  awk '{if (NR>2) print $1}' $FIRST.fc.tpm.tsv | sort -k1,1 >> $TAG.names
  CNT_LIST=`echo $PP | perl -ne '@t=split/\s+/; foreach $i (@t) {print "$i.counts.tmp "}; print "\n"'`
  TPM_LIST=`echo $PP | perl -ne '@t=split/\s+/; foreach $i (@t) {print "$i.TPM.tmp "}; print "\n"'`
  paste $TAG.names $CNT_LIST > ../exp_tables/$TAG.counts.tsv
  paste $TAG.names $TPM_LIST > ../exp_tables/$TAG.TPM.tsv
 
  rm *.names *.counts.tmp *TPM.tmp
  
  ## annotate tables with gene names, gene types, and plasmid/chromosome/prophage (optional)
  $SDIR/script/simple_exp_tables.pl $WDIR $TAG
  
  echo "ALL EXPRESSION TABLE PROCESSING IS DONE!"
fi
