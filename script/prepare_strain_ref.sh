#!/bin/bash 

source activate prokka 
set -euo pipefail

SDIR=$1
WDIR=$2
TAG=$3
CPUS=$4
REF=$5

FA=$TAG.fa
PROPHAGE=$TAG.prophage.bed

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

if [[ $REF == "" ]] 
then
  echo "==> No extra reference was specified, ncRNAs will be ignored"
elif [[ $REF == "prokka" ]]
then
  echo "==> Option \"-r prokka\" activated, ncRNAs will be identified using Prokka's built-it Rfam database"
else 
  echo "==> External reference fasta $REF will be used to annotate additional ncRNA and CDS"
fi 

######################################################

cp $FA $TAG.genome.fa
samtools faidx $TAG.genome.fa
cut -f 1,2 $TAG.genome.fa.fai > $TAG.chrom.sizes 

## annotate with Prokka - either with Rfam for ncRNAs, or blasting & parsing extra genes from reference fasta: 

if [[ -d $TAG.prokka ]]
then
  rm -rf $TAG.prokka
fi
 
if [[ $REF == "" ]]
then
  ## Rfam database here is OK, but still quite outdated. 
  ## It finds about 200 ncRNAs for average Salmonella; there are about 280 annotated right now.  
  echo "Running Prokka annotation; using --noanno option to only discover CDS"
  echo "Annotating noncoding RNAs using default Prokka Rfam database!" 
  prokka --noanno --cpus $CPUS --outdir $TAG.prokka --prefix $TAG.prokka --locustag ${TAG%%_*} --rfam $TAG.genome.fa &> /dev/null 
  grep -P "\tCDS\t" $TAG.prokka/$TAG.prokka.gff | sed "s/$/;gene_biotype=protein_coding;/g" > $TAG.CDS.gff
  grep -P "\tmisc_RNA\t" $TAG.prokka/$TAG.prokka.gff | sed "s/misc_RNA/ncRNA/g" | sed "s/$/;gene_biotype=noncoding_rna;/g" > $TAG.ncRNA.gff
  N_CDS=`grep -c -P "\tCDS\t" $TAG.CDS.gff`
  N_NCR=`grep -c -P "\tncRNA\t" $TAG.ncRNA.gff`
  echo
  echo "==> Found $N_CDS protein-coding (CDS) and $N_NCR non-coding RNA (misc_RNA/ncRNA) features."
  echo 
else 
  ## Make sure you have correct ncRNA names in the reference fasta - they will be used as a Name in GFF. 
  echo "Running Prokka annotation; using --noanno option to only discover CDS."
  echo "Annotating noncoding RNAs using blastn and custom reference file $REF!" 
 
  ## find all CDS using Prodigal's default settings - no Rfam search here 
  prokka --noanno --cpus $CPUS --outdir $TAG.prokka --prefix $TAG.prokka --locustag ${TAG%%_*} $TAG.genome.fa &> /dev/null 
  grep -P "\tCDS\t" $TAG.prokka/$TAG.prokka.gff | sed "s/$/;gene_biotype=protein_coding;/g" > $TAG.CDS.gff
  ## find all sORF and ncRNA
  makeblastdb -dbtype nucl -in $TAG.genome.fa -out ${TAG}_blast &> /dev/null 
  blastn -query $NC_REF -db ${TAG}_blast -evalue 1 -task megablast -outfmt 6 > $TAG.ncRNA_blast.out 2> /dev/null 

  ## new version of this script drops all mia- sORFs overlapping a Prokka CDS 
  $SDIR/script/make_gff_from_ref_blast.pl $NC_REF $TAG.CDS.gff ${TAG%%_*} $TAG.ncRNA_blast.out > $TAG.ncRNA.gff
  
  N_CDS=`grep -c -P "\tCDS\t" $TAG.CDS.gff`
  N_NCR=`wc -l $TAG.ncRNA.gff | awk '{print $1}'`
  echo
  echo "==> Found $N_CDS protein-coding (CDS) and $N_NCR non-coding RNA (misc_RNA/ncRNA) features."
  echo
  rm ${TAG}_blast.n* $TAG.ncRNA_blast.out
fi 

sed "s/\tCDS\t/\tgene\t/g"   $TAG.CDS.gff   >  $TAG.gene.gff
sed "s/\tncRNA\t/\tgene\t/g" $TAG.ncRNA.gff >> $TAG.gene.gff

echo "==> Files $TAG.genome.fa, $TAG.CDS.gff, $TAG.ncRNA.gff, and $TAG.gene.gff successfully generated"

## make STAR reference for small genome size

if [[ -d $TAG.STAR ]]
then
  rm -rf $TAG.STAR
fi
mkdir ${TAG}.STAR 
STAR --runThreadN $CPUS --runMode genomeGenerate --genomeDir ${TAG}.STAR --genomeFastaFiles $TAG.genome.fa --genomeSAindexNbases 10 &> /dev/null
mv Log.out $TAG.star.log 
echo "==> STAR aligner index $TAG.STAR successfully generated"

##make rRNA/tRNA interval file  
$SDIR/script/make_rrna_operon.pl $TAG.prokka/$TAG.prokka.gff $TAG.ncRNA.gff | sort -k1,1 -k2,2n | bedtools merge -i - > $TAG.rRNA.bed

## mv all to the ref dir 
mv $TAG.genome.fa $TAG.genome.fa.fai $TAG.chrom.sizes $WDIR/study_strains/$TAG
mv $TAG.gene.gff $TAG.CDS.gff $TAG.ncRNA.gff $WDIR/study_strains/$TAG
mv $TAG.prokka ${TAG}.STAR $TAG.star.log $WDIR/study_strains/$TAG
cp $PROPHAGE $WDIR/study_strains/$TAG/$TAG.prophage.bed
mv $TAG.rRNA.bed $WDIR/study_strains/$TAG

echo "==> All the generated files and indexes have been moved to $WDIR/study_strains/$TAG."
echo "==> Strain $TAG: all done generating reference!" 
