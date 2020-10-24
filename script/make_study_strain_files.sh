#!/bin/bash 

## this is ran in <wdir>/study_strains
set -eo pipefail

SDIR=$1
WDIR=$2
TAG=$3
CPUS=$4
REF=$5

FA=$TAG.fa

if [[ -d "$WDIR/study_strains/$TAG" ]]
then
  echo "==> Found $WDIR/study_strains/$TAG! Will add files to the existing directory."
  rm -rf $WDIR/study_strains/$TAG/*.STAR $WDIR/study_strains/$TAG/*.prokka 
else 
  echo "==> Directory $WDIR/study_strains/$TAG was not found and will be created." 
  mkdir $WDIR/study_strains/$TAG
fi

if [[ $CPUS == "" ]]
then 
  echo "==> Parallel jobs will be ran on 4 cores (default)."
  CPUS=4
fi

if [[ $REF == "" ]] 
then
  echo "==> No extra reference was specified, ncRNAs will be predicted using Prokka's --rfam."
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
  ## Probably does a lot worse for less sequenced/studied organisms
 
  echo "==> Running Prokka annotation; using --noanno option to discover CDS/rRNA/tRNA, and --rfam to discover ncRNA."
  prokka --noanno --cpus $CPUS --outdir $TAG.prokka --prefix $TAG.prokka --locustag ${TAG%%_*} --rfam $TAG.genome.fa &> /dev/null 
 
  grep -P "\tCDS\t" $TAG.prokka/$TAG.prokka.gff | sed "s/$/;gene_biotype=protein_coding;/g" > $TAG.CDS.gff
  grep -P "\tmisc_RNA\t" $TAG.prokka/$TAG.prokka.gff | sed "s/misc_RNA/ncRNA/g" | sed "s/$/;gene_biotype=noncoding_rna;/g" > $TAG.ncRNA.gff
  N_CDS=`grep -c -P "\tCDS\t" $TAG.CDS.gff`
  N_NCR=`grep -c -P "\tncRNA\t" $TAG.ncRNA.gff`
  echo
  echo "==> Found $N_CDS protein-coding (CDS) and $N_NCR non-coding RNA (misc_RNA/ncRNA) features."
  echo 
 
  ## TODO: Change this into reformatting Prokka -> United 
  cat $TAG.ncRNA.gff $TAG.CDS.gff | sort -k1,1 -k4,4n > $TAG.united.gff 
  if [[ -f $TAG.prophage.bed ]] 
  then 
    bedtools intersect -wo -a $TAG.united.gff -b $TAG.prophage.bed > $TAG.prophage_overlap.tsv
  fi 
else 
  ## Make sure you have correct ncRNA names in the reference fasta - they will be used as a Name in GFF. 
  echo "==> Running Prokka annotation; using --noanno option to discover CDS/rRNA/tRNA."

  ## find all CDS using Prodigal's default settings - no Rfam search here 
  prokka --noanno --cpus $CPUS --outdir $TAG.prokka --prefix $TAG.prokka --locustag ${TAG%%_*} $TAG.genome.fa &> /dev/null 

  makeblastdb -dbtype nucl -in $TAG.genome.fa -out ${TAG}.blast &> /dev/null 
  blastn -query $REF -db ${TAG}.blast -evalue 1 -task megablast -outfmt 6 > $TAG.ref_blast.out 2> /dev/null 

  ## starting with v0.6 this gets all new locus tags 
  $SDIR/script/unify_study_gff.pl $TAG.prokka/$TAG.prokka.gff $TAG.ref_blast.out $REF
  grep -P "\tCDS\t"   $TAG.united.gff > $TAG.CDS.gff
  grep -P "\tncRNA\t" $TAG.united.gff > $TAG.ncRNA.gff
  
  N_CDS=`grep -c -P "\tCDS\t" $TAG.united.gff`
  N_NCR=`grep -c -P "\tncRNA\t" $TAG.united.gff`
  echo
  echo "==> Found $N_CDS protein-coding (CDS) and $N_NCR non-coding RNA (misc_RNA/ncRNA) features."
  echo
  rm ${TAG}.blast.n*
  if [[ -f $TAG.prophage.bed ]] 
  then 
    bedtools intersect -wo -a $TAG.united.gff -b $TAG.prophage.bed > $TAG.prophage_overlap.tsv
  fi 
fi
 
## we need this format for featureCounts quant
perl -ne 's/\tCDS\t|\tother\t|\tncRNA\t|\trRNA\t|\ttRNA\t|\tmisc\t/\tgene\t/g; print' $TAG.united.gff > $TAG.gene.gff

echo "==> Files $TAG.genome.fa, $TAG.united.gff, and $TAG.gene.gff successfully generated."

## make STAR reference for small genome size

if [[ -d $TAG.STAR ]]
then
  rm -rf $TAG.STAR
fi
mkdir ${TAG}.STAR 

## STAR needs adjustments for small genomes, 8 is good for ~Ecoli-Salmonella
NSAI=`awk '{sum+=$2} END {printf "%d\n",log(sum)/(2*log(2))-3}' $TAG.genome.fa.fai`

STAR --runThreadN $CPUS --runMode genomeGenerate --genomeDir ${TAG}.STAR --genomeFastaFiles $TAG.genome.fa --genomeSAindexNbases $NSAI &> /dev/null
echo "==> STAR aligner index $TAG.STAR successfully generated"

##make rRNA/tRNA interval file  
$SDIR/script/make_rrna_operon.pl $TAG.prokka/$TAG.prokka.gff $TAG.united.gff | sort -k1,1 -k2,2n | bedtools merge -i - > $TAG.rRNA.bed
N_INT=`wc -l $TAG.rRNA.bed | awk '{print $1}'`
N_OP=`awk '$3-$2>3000' $TAG.rRNA.bed | wc -l`
echo "==> Making combined rRNA operon interval file: Identified $N_OP rRNA operons, $N_INT intervals overall." 

## mv all necessary files to the appropriate ref dirs 
mv $TAG.gene.gff $TAG.united.gff $TAG.CDS.gff $TAG.ncRNA.gff $TAG.rRNA.bed $WDIR/study_strains/$TAG
mv $TAG.prokka ${TAG}.STAR $TAG.genome.fa $TAG.genome.fa.fai $TAG.chrom.sizes $WDIR/study_strains/$TAG
mv $TAG.ref_blast.out $TAG.match.tsv $TAG.prophage_overlap.tsv $WDIR/study_strains/$TAG  || : 

echo "==> All the generated files and indexes have been moved to $WDIR/study_strains/$TAG."
echo "==> Strain $TAG: all reference files successfully generated!"
echo
echo 
