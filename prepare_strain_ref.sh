#!/bin/bash 

##   v0.2 - all necessary files are in one place (REFDIR/SPECIES)
##   no Prokka annotation of proteins - just find the CDS, Roary will do the rest. 
##   no Rsem/kallisto - they don't work correctly anyway. 

## next snippet is adapted from https://medium.com/@Drew_Stokes/bash-argument-parsing-54f3b81a6a8f

GRN='\033[1;32m'
GRN2='\033[0;32m'
RED='\033[1;31m'
BL='\033[0;34m'
NC='\033[0m' # No Color

if [[ $# < 5 ]]
then
  echo 
  printf "Step 1 of reference preparation: prepare reference for a single ${RED}study${NC} strain.\n"
  echo "============================================================================="
  printf "Usage: ${GRN}prepare_strain_ref.sh ${GRN2}<ref_directory> <strain_tag> <genome_fa> <prophage_bed> <rRNA_bed> [-p CPUs]${NC}\n"
  echo "       (to predict ncRNAs using Prokka's Rfam DB)"
  echo "       - or - " 
  printf "       ${GRN}prepare_strain_ref.sh ${GRN2}<ref_directory> <strain_tag> <genome_fa> <prophage_bed> <rRNA_bed> [-p CPUs] [-r ref_ncRNA_fasta]${NC}\n"
  echo "       (to assign ncRNAs by simply blasting the existing reference ncRNA fasta to the genome)"
  echo 
  exit 1
fi

PARAMS=""
NC_REF=""
CPUS=""

while (( "$#" )); do
  case "$1" in
    -r|--ref_ncrna)
      NC_REF=$2
      shift 2
      if [[ $NC_REF == "" ]]
      then
        echo "ERROR: -r flag requires a non-empty argument (reference ncRNA fasta file)!" 
        exit 1 
      fi
      echo "==> Invoking -r option: annotation of ncRNA by blasting the reference ncRNA fasta to genome."
      echo 
      printf "${RED}If you are using a custom-made ncRNA fasta file, please make sure sequence names are correct.${NC}\n"
      printf "${RED}They will be used as ncRNA gene names in the GTF file and the final expression table.${NC}\n"
      echo 
      ;;  
    -p|--cpus)
      CPUS=$2
      shift 2
      if [[ $CPUS == "" ]]
      then
        echo "ERROR: -p flag requires a non-empty argument (number of CPUs)!" 
        exit 1 
      fi
      echo "==> Invoking -p option: parallel jobs will be run on $CPUS cores."
      ;;  
    --) # end argument parsing
      shift
      break
      ;;  
    -*|--*=) # unsupported flags
      echo "ERROR: unsupported flag $1" >&2
      exit 1
      ;;  
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift
      ;;  
  esac
done
eval set -- "$PARAMS"

REFDIR=$1
TAG=$2
FA=$3
PROPHAGE=$4
RRNA=$5

if [[ -d "$REFDIR/$TAG" ]]
then
  echo "Found $REFDIR/$TAG! Will add files to the existing directory."
  rm -rf $REFDIR/$TAG/*.STAR $REFDIR/$TAG/*.prokka 
else 
  echo "Directory $REFDIR/$TAG was not found and will be created." 
  mkdir $REFDIR/$TAG
fi

if [[ $CPUS == "" ]]
  then 
  echo "==> Parallel jobs will be ran on 16 cores (default)."
  CPUS=16
fi


######################################################

source activate prokka 
set -euo pipefail

cp $FA $TAG.genome.fa
samtools faidx $TAG.genome.fa
cut -f 1,2 $TAG.genome.fa.fai > $TAG.chrom.sizes 

## annotate with Prokka - either 
if [[ $NC_REF == "" ]]
then
  ## Rfam database here is OK, but still quite outdated. 
  ## However, it does successfully find quite few leader peptides that are overlooked by Prokka otherwise. 
  echo "Running Prokka annotation; using --noanno option to only discover CDS."
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
  echo "Annotating noncoding RNAs using blastn and custom reference file $NC_REF!" 
  makeblastdb -dbtype nucl -in $TAG.genome.fa -out ${TAG}_blast &> /dev/null 
  blastn -query $NC_REF -db ${TAG}_blast -evalue 1 -task megablast -outfmt 6 > $TAG.ncRNA_blast.out 2> /dev/null 
  make_ncRNA_gff_from_blast.pl $NC_REF ${TAG%%_*} $TAG.ncRNA_blast.out > $TAG.ncRNA.gff
  
  ## note no --rfam option in this case  
  prokka --noanno --cpus $CPUS --outdir $TAG.prokka --prefix $TAG.prokka --locustag ${TAG%%_*} $TAG.genome.fa &> /dev/null 
  grep -P "\tCDS\t" $TAG.prokka/$TAG.prokka.gff | sed "s/$/;gene_biotype=protein_coding;/g" > $TAG.CDS.gff
  N_CDS=`grep -c -P "\tCDS\t" $TAG.CDS.gff`
  N_NCR=`wc -l $TAG.ncRNA.gff | awk '{print $1}'`
  echo
  echo "==> Found $N_CDS protein-coding (CDS) and $N_NCR non-coding RNA (misc_RNA/ncRNA) features."
  echo
  rm ${TAG}_blast.n* $TAG.ncRNA_blast.out
fi 

sed "s/\tCDS\t/\tgene\t/g"   $TAG.CDS.gff   >  $TAG.gene.gff
sed "s/\tncRNA\t/\tgene\t/g" $TAG.ncRNA.gff >> $TAG.gene.gff

echo "Files $TAG.genome.fa, $TAG.CDS.gff, $TAG.ncRNA.gff, and $TAG.gene.gff successfully generated."

## make STAR reference for small genome size 
mkdir ${TAG}.STAR 
STAR --runThreadN $CPUS --runMode genomeGenerate --genomeDir ${TAG}.STAR --genomeFastaFiles $TAG.genome.fa --genomeSAindexNbases 10 &> /dev/null
mv Log.out $TAG.star.log 

## mv all to the ref dir 
mv $TAG.genome.fa $TAG.genome.fa.fai $TAG.chrom.sizes $REFDIR/$TAG
mv $TAG.gene.gff $TAG.CDS.gff $TAG.ncRNA.gff $REFDIR/$TAG
mv $TAG.prokka ${TAG}.STAR $TAG.star.log $REFDIR/$TAG
cp $PROPHAGE $REFDIR/$TAG/$TAG.prophage.bed
cp $RRNA     $REFDIR/$TAG/$TAG.rRNA.bed 

echo "All the generated files and indexes have been moved to $REFDIR/$TAG."
echo "Strain $TAG: all done generating reference!" 
