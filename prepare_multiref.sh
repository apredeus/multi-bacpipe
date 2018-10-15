#!/bin/bash 

GRN='\033[1;32m'
GRN2='\033[0;32m'
RED='\033[1;31m'
NC='\033[0m' # No Color

if [[ $# < 3 ]]
then
  echo 
  printf "Step 2 of reference preparation:\n" 
  printf "prepare multi-reference for all ${RED}study${NC} and ${GRN}reference${NC} strains listed in the config file.\n"
  echo "======================================================================================"
  printf "Usage: ${GRN}prepare_multiref.sh ${GRN2}<working_directory> <ref_directory> <config> [-p CPUs]${NC}\n"
  echo 
  exit 1
fi

PARAMS=""
CPUS=""

while (( "$#" )); do
  case "$1" in
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

WDIR=$1
REFDIR=$2
CONFIG=$3

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

STUDY=`grep -v "^Reference" $WDIR/$CONFIG | cut -f 2 | sort | uniq`
REFSTR=`grep   "^Reference" $WDIR/$CONFIG | cut -f 2 | sort | uniq`

source activate roary 
cd $WDIR/refstr
## make sure process quits if inconsistencies are discovered 
set -euo pipefail

## check if reference strain GTF files in refstr are in Roary-friendly format (e.g. converted from NCBI to). 
## reqs: 1) unique IDs that are locus tags; 2) no names that are equal to ID; 3) only CDS features; 4) ##FASTA and genomic fa are present. 

for i in $REFSTR
do
  ref_check_roary_gff.sh $WDIR $i
done 

## check if every study strain has a non-empty folder - take this part from the mbc_check_config.sh 

for i in $STUDY
do
  ref_check_study_dir.sh $REFDIR $i
done  

## make study strain GFF to run roary 
for i in $STUDY
do
  cat $REFDIR/$i/$i.CDS.gff   >  0_$i.roary.gff
  echo "##FASTA"              >> 0_$i.roary.gff
  cat $REFDIR/$i/$i.genome.fa >> 0_$i.roary.gff
done 

## run Roary on all the strains present in the refstr

echo "==> Running Roary on all strains, using $CPUS cores..."
roary -p $CPUS -v -e -f roary *.roary.gff &> roary.log
mv roary $WDIR

## edit the presence-absence file. You should have dos2unix in your $PATH

echo "==> Roary pan-genome analysis is done! Performing some post-processing..."
cd $WDIR/roary 
dos2unix gene_presence_absence.csv
sed '1 s/.roary//g' gene_presence_absence.csv | sed '1 s/0_//g' > gene_presence_absence.csv.tmp 
mv gene_presence_absence.csv.tmp gene_presence_absence.csv

## make annotated CDS and ncRNA tables
cd $WDIR/refstr
echo "==> Generating a table of CDS orthologs."
annotate_CDS.pl $WDIR/roary/gene_presence_absence.csv $REFDIR $WDIR/mult.cfg > annotated_CDS.tsv
NCDS=`wc -l annotated_CDS.tsv | awk '{print $1}'`
echo "==> Table of $NCDS CDS orthologs successfully generated!"

echo "==> Generating a table of ncRNA orthologs."
annotate_ncRNA.pl $REFDIR $WDIR/mult.cfg > annotated_ncRNA.tsv
NRNA=`wc -l annotated_ncRNA.tsv | awk '{print $1}'`
echo "==> Table of $NRNA ncRNA orthologs successfully generated!"

echo "==> DONE generating multi-strain reference!" 
