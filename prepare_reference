#!/bin/bash 

GRN='\033[1;32m'
GRN2='\033[0;32m'
RED='\033[1;31m'
NC='\033[0m' # No Color

if [[ $# < 3 ]]
then
  echo 
  printf "One-command reference preparation:\n" 
  printf "prepare multi-reference for all ${RED}study${NC} and ${GRN}reference${NC} strains listed in the config file, or\n"
  printf "prepare simple reference using reference genome fasta and GFF annotation.\n"
  echo "=========================================================================================="
  printf "Usage: ${GRN}prepare_reference.sh ${GRN2}<working_directory> <config> [-p CPUs] [-r ref_fasta]${NC}\n"
  echo " - or - " 
  printf "Usage: ${GRN}prepare_reference.sh ${GRN2}<working_directory> --simple <genome_fa> <genome_gff> [-p CPUs]${NC}\n"
  echo 
  exit 1
fi

PARAMS=""
CPUS=""
NC_REF="" 
SIMPLE=false
REF_FNA=""
REF_GFF=""  

while (( "$#" )); do
  case "$1" in
    -r|--ref_fasta)
      NC_REF=$2
      shift 2
      if [[ $NC_REF == "" ]]
      then
        echo "ERROR: -r flag requires a non-empty argument (reference ncRNA/CDS fasta file)!" 
        exit 1
      fi
      echo "==> Invoking -r option: annotation of extra ncRNA/CDS by blasting the reference fasta to genome."
      echo 
      printf "${RED}If you are using a custom-made ncRNA/CDS fasta file, please make sure that format is correct.${NC}\n"
      printf "${RED}They will be used as ncRNA/CDS gene names and type in the GTF file and the final expression table.${NC}\n"
      echo 
      ;;
    -p)
      CPUS=$2
      shift 2
      if [[ $CPUS == "" ]]
      then
        echo "ERROR: -p flag requires a non-empty argument (number of CPUs)!" 
        exit 1
      fi
      echo "==> Invoking -p option: parallel jobs will be run on $CPUS cores."
      ;;
    --simple)
      REF_FNA=$2
      REF_GFF=$3 
      SIMPLE=true
      shift 3
      if [[ $REF_FNA == "" || $REF_GFF == "" ]]
      then
        echo "ERROR: --simple flag requires 2 non-empty arguments (reference genome fasta and GFF annotation)!" 
        exit 1
      fi
      echo "==> Invoking --simple option: using established reference genome sequence and GFF annotation."
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

WDIR=$1
CONFIG=$2

## make sure process quits if inconsistencies are discovered 
set -euo pipefail

## config file has to be in $WDIR; throw an error otherwise 
CONFIG=`basename $CONFIG`

if [[ ! -s $WDIR/$CONFIG ]]
then
  >&2 echo "ERROR: Config file must be located in the working directory!" 
  exit 1 
fi

if [[ $CPUS == "" ]]
  then
  echo "==> Parallel jobs will be ran on 16 cores (default)."
  CPUS=16
fi

## TODO: check if $STUDY and $REFSTR come up non-empty
## TODO: check if all works with Rfam workflow (names don't get extracted)  

STUDY=`grep -v "^Reference" $WDIR/$CONFIG | cut -f 2 | sort | uniq`
REFSTR=`grep   "^Reference" $WDIR/$CONFIG | cut -f 2 | sort | uniq`

echo -e "Following study strains will be processed:\n\n$STUDY\n"
echo -e "Following reference strains will be processed:\n\n$REFSTR\n"

source activate roary 

## prepare individual strain references for each study strain
## you have to have a genome file and prophage file in place 

cd $WDIR/study_strains

if [[ $NC_REF == "" ]]
then
  for i in $STUDY
  do
    prepare_strain_ref.sh $WDIR $i.fa $i.prophage.bed -p $CPUS
  done 
else 
  for i in $STUDY 
  do
    prepare_strain_ref.sh $WDIR $i.fa $i.prophage.bed -p $CPUS -r $NC_REF
  done
fi

## now check if that all worked out OK 
for i in $STUDY
do
  ref_check_study_dir.sh $WDIR $i
done  
echo
echo "ALL STUDY STRAIN FILES AND DIRS ARE OK!"  
echo

cd $WDIR/ref_strains

## check if reference strain GTF files in refstr are in Roary-friendly format (e.g. converted from NCBI to). 
## reqs: 1) unique IDs that are locus tags; 2) no names that are equal to ID; 3) only CDS features; 4) ##FASTA and genomic fa are present. 

for i in $REFSTR
do
  ref_check_roary_gff.sh $WDIR $i
done 
echo
echo "ALL REFERENCE STRAIN GFF FILES ARE OK!"  
echo

## make study strain GFF to run roary 
for i in $STUDY
do
  cat $WDIR/study_strains/$i/$i.CDS.gff   >  0_$i.roary.gff
  echo "##FASTA"                          >> 0_$i.roary.gff
  cat $WDIR/study_strains/$i/$i.genome.fa >> 0_$i.roary.gff
done 

## run Roary on all the strains present in the refstr

echo "==> Running Roary on all strains, using $CPUS cores..."
roary -p $CPUS -v -e -f roary *.gff &> roary.log
rm 0*.roary.gff
mv roary roary.log $WDIR

## edit the presence-absence file. You should have dos2unix in your $PATH

echo "==> Roary pan-genome analysis is done! Performing some post-processing..."
cd $WDIR/roary 
dos2unix gene_presence_absence.csv
sed '1 s/.roary//g' gene_presence_absence.csv | sed '1 s/0_//g' > gene_presence_absence.csv.tmp 
mv gene_presence_absence.csv.tmp gene_presence_absence.csv

## make annotated CDS and ncRNA tables
cd $WDIR
echo "==> Generating a table of CDS orthologs."
annotate_CDS.pl $WDIR/roary/gene_presence_absence.csv $WDIR/study_strains $WDIR/$CONFIG > annotated_CDS.tsv
NCDS=`wc -l annotated_CDS.tsv | awk '{print $1}'`
echo "==> Table of $NCDS CDS orthologs successfully generated!"

echo "==> Generating a table of ncRNA orthologs."
annotate_ncRNA.pl $WDIR/study_strains $WDIR/$CONFIG > annotated_ncRNA.tsv
NRNA=`wc -l annotated_ncRNA.tsv | awk '{print $1}'`
echo "==> Table of $NRNA ncRNA orthologs successfully generated!"

## TODO: print stats - how many genes are common for all study strains, how many are unique, how many are seen in 2 and more strains.

echo "==> DONE generating multi-strain reference!" 
