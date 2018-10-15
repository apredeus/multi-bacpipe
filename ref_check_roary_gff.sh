#!/bin/bash

WDIR=$1
TAG=$2

cd $WDIR/refstr

if [[ -s $TAG.roary.gff ]]
then 
  echo "==> Checking reference strain $TAG: non-empty file $TAG.roary.gff found."
else 
  >&2 echo "ERROR: file $TAG.roary.gff is not found in refstr directory!"
  exit 1
fi

FASTA=`grep "##FASTA" $TAG.roary.gff`
N1=`grep -c -P "\t" $TAG.roary.gff`
N2=`grep -c -P "\tCDS\t" $TAG.roary.gff`
N3=`perl -ne '@t=split/\t+/; if ($t[8]=~m/ID=(.*?);/) {print "$1\n"}' $TAG.roary.gff | sort | uniq | wc -l`
NAMES=`perl -ne '@t=split/\t+/; if ($t[8]=~m/ID=(.*?);.*Name=(.*?);/) {print "$1\t$2\n"}' $TAG.roary.gff | awk '$1==$2'`

if [[ $FASTA != "" && $N1 == $N2 && $N1 == $N3 && $NAMES == "" ]]
then 
  echo "==> Checking file $TAG.roary.gff for pipeline compatibility: $N1 unique locus tags found."
else
  >&2 echo "ERROR: it seems like there are some problems with $TAG.roary.gff!"
  >&2 echo "Please make sure that: 1) all CDS IDs are unique;"
  >&2 echo "2) file is in Roary format, i.e. has GTF file followed by ##FASTA delimiter, and genome fasta;"
  >&2 echo "3) no Names are equal to IDs - these could be set using reference strain (see README)."
  exit 1
fi 
