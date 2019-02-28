#!/bin/bash 

## this is ran in <wdir>/ref_strains
source activate prokka 
set -euo pipefail

SDIR=$1
WDIR=$2
TAG=$3
REF=$4

FA=$TAG.fa
GFF=$TAG.gff 

if [[ -d "$WDIR/ref_strains/$TAG" ]]
then
  echo "Found $WDIR/ref_strains/$TAG! Will add files to the existing directory."
else 
  echo "Directory $WDIR/ref_strains/$TAG was not found and will be created." 
  mkdir $WDIR/ref_strains/$TAG
fi

if [[ $REF == "" ]] 
then
  echo "==> No extra reference was specified, ncRNAs will be ignored"
else 
  echo "==> External reference fasta $REF will be used to annotate additional ncRNA and CDS"
fi 

######################################################

cp $FA $TAG.genome.fa

if [[ $REF == "" ]]
then
  ## clean up the GFF, get file $TAG.clean.gff (ncRNA/CDS/pseudo == longest CDS) 
  $SDIR/script/reference_gff_cleanup.pl $TAG.gff > $TAG.clean.gff 
 
  grep -P "\tCDS\t" $TAG.clean.gff > $TAG.roary.gff 
  echo "##FASTA" >> $TAG.roary.gff 
  cat $TAG.genome.fa >> $TAG.roary.gff 
  echo "==> Files $TAG.genome.fa, $TAG.cleaned.gff, and $TAG.roary.gff successfully generated"
else 
  ## with REF blasting 

  ## clean up the GFF, get file $TAG.clean.gff (ncRNA/CDS/pseudo == longest CDS) 
  $SDIR/script/reference_gff_cleanup.pl $TAG.gff > $TAG.clean.gff
  
  grep -P "\tCDS\t" $TAG.clean.gff > $TAG.roary.gff 
  echo "##FASTA" >> $TAG.roary.gff 
  cat $TAG.genome.fa >> $TAG.roary.gff 
  
  makeblastdb -dbtype nucl -in $TAG.genome.fa -out ${TAG}.blast &> /dev/null
  blastn -query $REF -db ${TAG}.blast -evalue 1 -task megablast -outfmt 6 > $TAG.ref_blast.out 2> /dev/null

  ## this script matches existing locus tags to Blast results, generating $TAG.match.tsv table 
  ## CDS could only match CDS, and ncRNA could only match ncRNA  
  $SDIR/script/match_reference_gff.pl $TAG.clean.gff $TAG.ref_blast.out $REF $TAG
  echo "==> Files $TAG.genome.fa, $TAG.cleaned.gff, $TAG.ref_blast.out, $TAG.match.tsv, and $TAG.roary.gff successfully generated"
fi
 
## mv all to the ref dir 
mv $TAG.genome.fa $TAG.roary.gff $TAG.clean.gff $TAG.ref_blast.out $TAG.match.tsv $WDIR/ref_strains/$TAG || :

echo "==> All the generated files and indexes have been moved to $WDIR/ref_strains/$TAG."
echo "==> Strain $TAG: all reference files successfully generated!" 
