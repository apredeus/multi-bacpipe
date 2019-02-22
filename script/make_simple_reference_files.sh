#!/bin/bash 

source activate prokka 
set -euo pipefail

SDIR=$1
WDIR=$2
TAG=$3
CPUS=$4

if [[ -d "$WDIR/study_strains/$TAG" ]]
then
  echo "Found $WDIR/study_strains/$TAG! Will add files to the existing directory."
  rm -rf $WDIR/study_strains/$TAG/*.STAR $WDIR/study_strains/$TAG/*.prokka 
else 
  echo "Directory $WDIR/study_strains/$TAG was not found and will be created." 
  mkdir $WDIR/study_strains/$TAG
fi

if [[ $CPUS == "" ]]
  then 
  echo "==> Parallel jobs will be ran on 4 cores (default)."
  CPUS=4
fi


######################################################

## just like with prepare_strain_ref.sh, all this is happening in /study_strains 

cp $TAG.fa $TAG.genome.fa
samtools faidx $TAG.genome.fa
cut -f 1,2 $TAG.genome.fa.fai > $TAG.chrom.sizes 

echo "==> Running Prokka annotation to identify tRNA and rRNA loci"

if [[ -d $TAG.prokka ]]
then
  rm -rf $TAG.prokka
fi
prokka --noanno --cpus $CPUS --outdir $TAG.prokka --prefix $TAG.prokka --locustag ${TAG%%_*} $TAG.genome.fa &> /dev/null 

## this will write $TAG.CDS.gff and $TAG.ncRNA.gff with all the correct attributes, and account for tmRNA, misc_RNA etc
$SDIR/script/simple_gff_cleanup.pl $TAG.gff > $TAG.gene.gff

echo "==> Files $TAG.genome.fa and $TAG.gene.gff successfully generated"

## make STAR reference for small genome size 
if [[ -d $TAG.STAR ]]
then
  rm -rf $TAG.STAR
fi
mkdir $TAG.STAR 
STAR --runThreadN $CPUS --runMode genomeGenerate --genomeDir $TAG.STAR --genomeFastaFiles $TAG.genome.fa --genomeSAindexNbases 10 &> /dev/null
mv Log.out $TAG.star.log 
echo "==> STAR aligner index $TAG.STAR successfully generated"

##make rRNA/tRNA interval file  
$SDIR/script/make_rrna_operon.pl $TAG.prokka/$TAG.prokka.gff $TAG.gene.gff | sort -k1,1 -k2,2n | bedtools merge -i - > $TAG.rRNA.bed
echo "==> rRNA/tRNA operon interval file $TAG.rRNA.bed successfully created"

## mv all to the ref dir 
mv $TAG.genome.fa $TAG.genome.fa.fai $TAG.chrom.sizes $TAG.gene.gff $WDIR/study_strains/$TAG
mv $TAG.prokka ${TAG}.STAR $TAG.star.log $WDIR/study_strains/$TAG
mv $TAG.rRNA.bed $WDIR/study_strains/$TAG

echo "All the generated files and indexes have been moved to $WDIR/study_strains/$TAG."
echo "Strain $TAG: all done generating reference!" 
