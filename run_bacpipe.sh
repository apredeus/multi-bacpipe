#!/bin/bash 

REFDIR=$1
SPECIES=$2 ## e.g. P125109_ncbi
CPUS=$3

if [[ -d fastqs && "$(ls -A fastqs)" ]]; then
  echo "Found non-empty directory named fastqs! Continuing.."
else
  echo "ERROR: directory fastqs does not exist and is empty!"
  exit 1
fi

if [[ ! -d bams || ! -d stats || ! -d tdfs || ! -d RSEM || ! -d featureCounts ]]
then
  echo "One of the required directories is missing, I will try to create them..."
  mkdir bams tr_bams stats tdfs RSEM kallisto FastQC STAR_logs
else 
  echo "All the necessary directories found, continuing..." 
fi 

cp ~/bacpipe/*sh .

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

echo "Step 1: Running bowtie2.."
./prun_bowtie2.sh $REFDIR $SPECIES $CPUS
echo "Step 2: Making TDF files.." 
./pmake_tdf.sh $SPECIES
echo "Step 3: Running Picard..:"
./prun_picard.sh $REFDIR 
echo "Step 4: Running featureCounts.."
./prun_fcount.sh $REFDIR $SPECIES 
echo "Step 5: Running bacstat.."
./prun_bacstat.sh 

cd stats 
FLAG=`strand | tr "%" " " | awk '{sum+=$2} END {x=sum/NR; if (x<3) {print "RF"} else if (x>97) {print "FR"} else if (x>47 && x<53) {print "NONE"} else {print "ERROR"}}'`
if [[ $FLAG == "ERROR" ]]
then
  echo "ERROR: something is very much off with the strand-specificity of your RNA-seq!"
  exit 1
else 
  echo "The strandedness of your experiment was determined to be $FLAG"
fi
cd .. 
echo "Step 6: Running RSEM.." 

echo "Step 6: Running RSEM.."
./prun_rsem.sh $REFDIR $SPECIES $FLAG $CPUS

echo "All processing is now complete!"
