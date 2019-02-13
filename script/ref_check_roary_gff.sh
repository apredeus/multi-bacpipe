#!/bin/bash

WDIR=$1
TAG=$2

cd $WDIR/ref_strains

if [[ -s $TAG.gff ]]
then 
  echo "==> Checking reference strain $TAG:"
  echo "- non-empty file $TAG.gff found."
else 
  >&2 echo "ERROR: file $TAG.gff is not found in refstr directory!"
  exit 1
fi

FASTA=`grep "##FASTA" $TAG.gff`
N1=`grep -c -P "\t" $TAG.gff`
N2=`grep -c -P "\tCDS\t" $TAG.gff`
N3=`perl -ne '@t=split/\t+/; if ($t[8]=~m/ID=(.*?);/) {print "$1\n"}' $TAG.gff | sort | uniq | wc -l`
NAMES=`perl -ne '@t=split/\t+/; if ($t[8]=~m/ID=(.*?);.*Name=(.*?);/) {print "$1\t$2\n"}' $TAG.gff | awk '$1==$2'`

if [[ $FASTA != "" && $N1 == $N2 && $N1 == $N3 && $NAMES == "" ]]
then 
  echo "- $N1 unique locus tags found in $TAG.gff."
else
  >&2 echo "ERROR: it seems like there are some problems with $TAG.gff!"
  >&2 echo "Please make sure that: 1) all CDS IDs are unique;"
  >&2 echo "2) file is in Roary format, i.e. has GTF file followed by ##FASTA delimiter, and genome fasta;"
  >&2 echo "3) no Names are equal to IDs - these could be set using reference strain (see README)."
  exit 1
fi
