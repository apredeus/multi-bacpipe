#!/bin/bash

WDIR=$1
TAG=$2

cd $WDIR/ref_strains/$TAG

if [[ -s $TAG.clean.gff ]]
then 
  echo "==> Checking reference strain $TAG: non-empty file $TAG.clean.gff found."
else 
  >&2 echo "ERROR: file $TAG.clean.gff is not found in /ref_strains/$TAG directory!"
  exit 1
fi

N1=`grep -c -P "\t" $TAG.clean.gff`
N2=`perl -ne 'if (m/ID=(.*?);/) {print "$1\n"}' $TAG.clean.gff | sort | uniq | wc -l`

if [[ $N1 == $N2 ]]
then 
  echo "==> $N1 unique locus tags found in $TAG.clean.gff"
else
  >&2 echo "ERROR: it seems like there are some problems with $TAG.clean.gff!"
  >&2 echo "Please make sure that each feature has a unique locus tag in the ID= field."
  exit 1
fi
