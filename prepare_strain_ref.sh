#!/bin/bash 

REFDIR=$1
TAG=$2
FA=$3
CPUS=$4

##   v0.2 - all necessary files are in one place (REFDIR/SPECIES)
##   no Prokka annotation of proteins - just find the CDS, Roary will do the rest. 
##   no Rsem/kallisto - they don't work correctly anyway. 


if [[ $# != 4 ]]
then
  echo "Usage: prepare_strain_ref.sh <ref_directory> <strain_tag> <genome_fa> <CPUs>"
  exit 1
fi

if [[ -d "$REFDIR/$TAG" ]]
then
  echo "Found $REFDIR/$TAG! Going to remove it in order to prevent possible version conflicts." 
  rm -rf $REFDIR/$TAG
fi
mkdir $REFDIR/$TAG
 
cp $FA $TAG.genome.fa
samtools faidx $TAG.genome.fa
cut -f 1,2 $TAG.genome.fa.fai > $TAG.chrom.sizes 

source activate prokka2 
## I will add custom Rfam DB here 
echo "Running Prokka annotation; using --noanno option to only discover CDS; annotate noncoding RNAs with custom Rfam DB." 
prokka --noanno --cpus $CPUS --outdir $TAG.prokka --prefix $TAG.prokka --locustag ${TAG%%_*} --rfam $TAG.genome.fa &> /dev/null 
source deactivate

grep -P "\tCDS\t" $TAG.prokka/$TAG.prokka.gff | sed "s/$/;gene_biotype=protein_coding;/g" > $TAG.CDS.gff
grep -P "\tmisc_RNA\t" $TAG.prokka/$TAG.prokka.gff | sed "s/misc_RNA/ncRNA/g" | sed "s/$/;gene_biotype=noncoding_rna;/g" > $TAG.ncRNA.gff
N_CDS=`grep -c -P "\tCDS\t" $TAG.CDS.gff`
N_NCR=`grep -c -P "\tncRNA\t" $TAG.ncRNA.gff`
echo "Found $N_CDS protein-coding (CDS) and $N_NCR non-coding RNA (misc_RNA/ncRNA) features!"
sed "s/\tCDS\t/\tgene\t/g"   $TAG.CDS.gff   >  $TAG.gene.gff
sed "s/\tncRNA\t/\tgene\t/g" $TAG.ncRNA.gff >> $TAG.gene.gff

echo "Files $TAG.genome.fa, $TAG.CDS.gff, $TAG.ncRNA.gff, and $TAG.gene.gff successfully generated."

## make STAR reference for small genome size 
mkdir ${TAG}.STAR 
STAR --runThreadN $CPUS --runMode genomeGenerate --genomeDir ${TAG}.STAR --genomeFastaFiles $TAG.genome.fa --genomeSAindexNbases 10 &> $TAG.star.log 

## mv all to the ref dir 
mv $TAG.genome.fa $TAG.genome.fa.fai $TAG.chrom.sizes $REFDIR/$TAG
mv $TAG.gene.gff $TAG.CDS.gff $TAG.ncRNA.gff $REFDIR/$TAG
mv $TAG.prokka ${TAG}.STAR $TAG.star.log $REFDIR/$TAG

echo "All the generated files and indexes have been moved to $REFDIR/$TAG."
echo "Strain $TAG: all done generating reference!" 
