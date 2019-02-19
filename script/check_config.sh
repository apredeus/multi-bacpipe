#!/bin/bash 

## this is where we shall check 
## - if fastq exists and is non-empty; 
## - if all other necessary dirs exist or need to be created; 
## - if all the reference dirs exist and contain all the necessary files; 
## - if roary dir exists and has all the necessary files

set -euo pipefail

SDIR=$1
WDIR=$2
CONFIG=$3
MULTI="" 

cd $WDIR

## this prevents pipefail from stopping script execution on grep coming up with nothing
MULTI=`grep -c "^Reference" $CONFIG || true`

if [[ -d fastqs && "$(ls -A fastqs)" ]]
then
  echo "==> Found non-empty directory named fastqs! Continuing.."
else
  >&2 echo "ERROR: directory fastqs does not exist or is empty!"
  exit 1
fi

if [[ $MULTI != 0 ]]
then
  echo "==> Processing workflow: MULTI-STRAIN. Continuing.." 
  if [[ -d roary && "$(ls -A roary)" ]]
  then ## add extra checks for properly parsed Roary output 
    echo "==> Found non-empty directory named roary! Continuing.."
  else
    >&2 echo "ERROR: directory roary does not exist or is empty! Please make sure you've ran \"prepare_bacterial_reference\" on your config file."
    exit 1
  fi

  STUDY=`grep -v "^Reference" $CONFIG | cut -f 2 | sort | uniq`
  REFSTR=`grep   "^Reference" $CONFIG | cut -f 2 | sort | uniq`

  echo -e "Following study strains are found:\n\n$STUDY\n"
  echo -e "Following reference strains are found:\n\n$REFSTR\n"

  ## let's check again if our study strains check out  
  for i in $STUDY
  do
    $SDIR/script/ref_check_study_dir.sh $WDIR $i
  done
  echo
  echo "ALL STUDY STRAIN FILES AND DIRS ARE OK!"  
  echo

  cd $WDIR/ref_strains

  ## check if reference strain GTF files in refstr are in Roary-friendly format (e.g. all NCBI files were cleaned up). 
  ## reqs: 1) unique IDs that are locus tags; 2) no names that are equal to ID; 3) only CDS features; 4) ##FASTA and genomic fa are present. 

  for i in $REFSTR
  do
    $SDIR/script/ref_check_roary_gff.sh $WDIR $i
  done
  echo
  echo "ALL REFERENCE STRAIN GFF FILES ARE OK!"  
  echo


else 
  echo "==> Processing workflow: SIMPLE SINGLE-STRAIN. Continuing.."
  TAG=`cut -f 2 $CONFIG | uniq`

  echo -e "Following study strain is found:\n\n$TAG\n"

  ## let's check again if our study strains check out  
  $SDIR/script/ref_check_study_dir.sh $WDIR $TAG

  echo
  echo "STUDY STRAIN FILES AND DIRS ARE OK!"  
  echo
fi 

if [[ ! -d bams || ! -d stats || ! -d strand || ! -d tdfs_and_bws || ! -d featureCounts || ! -d FastQC || ! -d exp_tables  ]]  
then
  echo "==> One (or more) of the required directories is missing, I will try to create them.."
  mkdir -p featureCounts FastQC bams stats strand tdfs_and_bws exp_tables
else
  echo "==> All the necessary directories found, continuing." 
fi

echo "ALL DIRECTORY TESTING WAS PASSED SUCCESSFULLY!"
