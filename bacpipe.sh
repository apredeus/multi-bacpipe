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

if [[ ! -d bams || ! -d stats || ! -d strand || ! -d tdfs_and_bws || ! -d RSEM || ! -d featureCounts || ! -d FastQC || ! -d rtRNA  ]]
then
  echo "One of the required directories is missing, I will try to create them..."
  mkdir bams stats strand tdfs_and_bws RSEM featureCounts FastQC rtRNA 
else 
  echo "All the necessary directories found, continuing..." 
fi 

#cp ~/bacpipe/*sh .

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

echo "Step 1: Running FastQC.."
#./prun_fastqc.sh
echo 
echo "=================================================================================="
echo

echo "Step 2: Running rRNA evaluation and alignment.."
#./prun_bowtie2.sh $REFDIR $SPECIES $CPUS
echo 
echo "=================================================================================="
echo

echo "Step 3: Making TDF and strand-specific bigWig files.." 
#./prun_coverage.sh $REFDIR $SPECIES $CPUS
echo 
echo "=================================================================================="
echo

echo "Step 4: Running featureCounts on all possible strand settings.."
#./prun_strand.sh $REFDIR $SPECIES $CPUS
echo 
echo "=================================================================================="
echo

echo "Step 5: Calculating strandedness and other statistics.."
#./prun_stat.sh $CPUS
echo 
echo "=================================================================================="
echo

cd stats
cat *rnastat | awk 'BEGIN {min=100;max=0} {sum+=$16; if($16>max) \
{max=$16}; if($16<min) {min=$16};} END {print "Average percent of \
reads matching the coding strand: "sum/NR", lowest: "min", highest: "max}'

STRAND=`cat *rnastat | awk '{sum+=$16} END {x=sum/NR; if (x<10) \
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

echo "Step 6: Running featureCounts on normal and extended annotation.."
./prun_fcount.sh $REFDIR $SPECIES $CPUS $STRAND
echo 
echo "=================================================================================="
echo

echo "Step 7: Running RSEM on normal and extended annotation.."
#./prun_rsem.sh $REFDIR $SPECIES $CPUS $STRAND
echo 
echo "=================================================================================="
echo

echo "All processing is now complete!"
