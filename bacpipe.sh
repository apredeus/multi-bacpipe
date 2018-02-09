#!/bin/bash 

REFDIR=$1
SPECIES=$2 ## e.g. P125109_ncbi
CPUS=$3
WDIR=`pwd`

if [[ $# != "3" ]] 
then 
  echo "What is the airspeed velocity of an unladen swallow?"
  #echo "For more usage information, please see https://github.com/apredeus/bacpipe"
  exit 1
fi

cd $WDIR 

if [[ -d fastqs && "$(ls -A fastqs)" ]]; then
  echo "Found non-empty directory named fastqs! Continuing.."
else
  echo "ERROR: directory fastqs does not exist and is empty!"
  exit 1
fi

if [[ ! -d bams || ! -d stats || ! -d strand || ! -d tdfs_and_bws || \
! -d RSEM || ! -d featureCounts || ! -d FastQC || \
! -d kallisto || ! -d exp_tables || ! -d cleaned_fastqs ]]
then
  echo "One of the required directories is missing, I will try to create them..."
  mkdir bams stats strand tdfs_and_bws RSEM exp_tables
  mkdir featureCounts FastQC kallisto cleaned_fastqs
else
  echo "All the necessary directories found, continuing..." 
fi
echo 
echo "=================================================================================="
echo

if [[ $SPECIES == "" || $REFDIR == "" ]]
then
  echo "ERROR: You have to specify REFDIR and SPECIES!"
  exit 1
fi

if [[ $CPUS == "" ]]
	then 
  echo "Parallel jobs have been set to default - running on 16 cores."
  CPUS=16
else 
  echo "Parallel jobs will be ran on $CPUS cores."
fi

echo "["`date +%H:%M:%S`"] Step 1: Running FastQC.."
prun_fastqc.sh $WDIR $CPUS
echo 
echo "=================================================================================="
echo

echo "["`date +%H:%M:%S`"] Step 2: Running rRNA evaluation and alignment.."
prun_bowtie2.sh $WDIR $REFDIR $SPECIES $CPUS
echo 
echo "=================================================================================="
echo

echo "["`date +%H:%M:%S`"] Step 3: Making TDF and strand-specific bigWig files.." 
prun_coverage.sh $WDIR $REFDIR $SPECIES $CPUS
echo 
echo "=================================================================================="
echo

echo "["`date +%H:%M:%S`"] Step 4: Running featureCounts on all possible strand settings.."
prun_strand.sh $WDIR $REFDIR $SPECIES $CPUS
echo 
echo "=================================================================================="
echo

echo "["`date +%H:%M:%S`"] Step 5: Calculating strandedness and other statistics.."
prun_stats.sh $WDIR
echo 
echo "=================================================================================="
echo

cd stats
cat *.strand | awk 'BEGIN {min=100;max=0} {sum+=$16; if($16>max) \
{max=$16}; if($16<min) {min=$16};} END {print "Average percent of \
reads matching the coding strand: "sum/NR", lowest: "min", highest: "max}'

STRAND=`cat *.strand | awk '{sum+=$16} END {x=sum/NR; if (x<10) \
{print "RF"} else if (x>90) {print "FR"} else if (x>45 && x<55) \
{print "NONE"} else {print "ERROR"}}'` 

if [[ $STRAND == "ERROR" ]]
then
  echo "ERROR: something is very much off with the strand-specificity of your RNA-seq!"
  exit 1
else 
  echo "The strandedness of your experiment was determined to be $STRAND"
fi
cd .. 

echo 
echo "=================================================================================="
echo

echo "["`date +%H:%M:%S`"] Step 6: Running featureCounts on normal and extended annotation.."
prun_fcount.sh $WDIR $REFDIR $SPECIES $CPUS $STRAND
echo 
echo "=================================================================================="
echo

echo "["`date +%H:%M:%S`"] Step 7: Running kallisto on normal and extended annotation.."
prun_kallisto.sh $WDIR $REFDIR $SPECIES $CPUS $STRAND
echo 
echo "=================================================================================="
echo

echo "["`date +%H:%M:%S`"] Step 8: Running RSEM on normal and extended annotation.."
prun_rsem.sh $WDIR $REFDIR $SPECIES $CPUS $STRAND
echo 
echo "=================================================================================="
echo

echo "["`date +%H:%M:%S`"] Step 9: Making final expression tables.."
make_tables.sh $WDIR $REFDIR $SPECIES $CPUS
echo 
echo "=================================================================================="
echo


echo "["`date +%H:%M:%S`"] ALL YOUR BASE ARE BELONG TO US!!!"
