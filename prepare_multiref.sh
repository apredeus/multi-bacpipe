#!/bin/bash 

WDIR=$1
REFDIR=$2
CONFIG=$3

GRN='\033[1;32m'
GRN2='\033[0;32m'
RED='\033[1;31m'
BL='\033[0;34m'
NC='\033[0m' # No Color

if [[ $# < 3 ]]
then
  echo 
  printf "Step 2 of reference preparation:\n" 
  printf "prepare multi-reference for all ${RED}study${NC} and ${GRN}reference${NC} strains listed in the config file.\n"
  echo "======================================================================================"
  printf "Usage: ${GRN}prepare_multiref.sh ${GRN2}<working_directory> <ref_directory> <config> ${NC}\n"
  echo 
  exit 1
fi

## config file has to be in $WDIR; throw an error otherwise 
CONFIG=`basename $CONFIG`
if [[ ! -s $WDIR/$CONFIG ]]
then
  >&2 echo "ERROR: Config file must be located in the working directory!" 
  exit 1 
fi 

STUDY=`grep -v "^Reference" $WDIR/$CONFIG | cut -f 2 | sort | uniq`
REFSTR=`grep -v "Reference" $WDIR/$CONFIG | cut -f 2 | sort | uniq`

source activate roary 
cd $WDIR/refstr

## check if reference strain GTF files in refstr are in Roary-friendly format (e.g. converted from NCBI to). 
## reqs: 1) unique IDs that are locus tags; 2) no names that are equal to ID; 3) only CDS features; 4) ##FASTA and genomic fa are present. 

for i in $STUDY
do
  ref_check_roary_gff.sh $i
done 

## check if every study strain has a non-empty folder - take this part from the mbc_check_config.sh 

for i in $REFSTR
do
  ref_check_study_dir.sh $REFDIR $i
done  

## make study strain GFF to run roary 
for i in $STUDY
do
  cat $REFDIR/$TAG/$TAG.CDS.gff   >  0_$TAG.roary.gff
  echo "##FASTA"                  >> 0_$TAG.roary.gff
  cat $REFDIR/$TAG/$TAG.genome.fa >> 0_$TAG.roary.gff
done 

## run Roary on all the strains present in the refstr

roary -p $CPUS -v -e -f roary *.roary.gff &> roary.log
mv roary $WDIR

## edit the presence-absence file. You should have dos2unix in your $PATH

cd $WDIR/roary 
dos2unix gene_presence_absence.csv
## !!!!!!!!!!!!!! add extra parsing here 


## make annotated CDS and ncRNA tables 
cd $WDIR/refstr
annotate_CDS.pl $WDIR/roary/gene_presence_absence.csv $REFDIR $WDIR/mult.cfg > annotated_CDS.tsv
annotate_ncRNA.pl $REFDIR $WDIR/mult.cfg > annotated_ncRNA.tsv
