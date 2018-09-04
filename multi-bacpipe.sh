#!/bin/bash 

REFDIR=$1
CONFIG=$2 ## tab-separated config file: tag to species ID 
CPUS=$3
WDIR=`pwd`

if [[ $# != "3" ]] 
then 
  echo "Usage: multi-bacpipe.sh <ref_dir> <config> <CPUs>" 
  echo "See github.com/apredeus/multi-bacpipe for more information."
  exit 1
fi

echo "=================================================================================="
echo "=================================================================================="
echo "===                                                                            ==="
echo "===                             Welcome to BACPIPE!                            ==="
echo "===                              PAY OR I PLAY (R)                             ==="
echo "===  For more information, please visit https://github.com/apredeus/bacpipe    ==="
echo "===                          Publication in preparation.                       ==="
echo "===                                                                            ==="
echo "=================================================================================="
echo "=================================================================================="
echo
echo

cd $WDIR 

if [[ $CONFIG == "" || $REFDIR == "" ]]
then
  echo "ERROR: You have to specify REFDIR and CONFIG!"
  exit 1
fi

if [[ $CPUS == "" ]]
	then 
  echo "Parallel jobs have been set to default - running on 16 cores."
  CPUS=16
else 
  echo "Parallel jobs will be ran on $CPUS cores."
fi

echo 
echo "=================================================================================="
echo

echo "["`date +%H:%M:%S`"] Step 0: Checking configuration file and reference availability.."
cd $WDIR
## add actually QUITTING on errors!
mbc_check_config.sh $WDIR $REFDIR $CONFIG
echo 
echo "=================================================================================="
echo

echo "["`date +%H:%M:%S`"] Step 1: Running FastQC.."
cd $WDIR/fastqs 
mbc_prun_fastqc.sh $WDIR $CPUS
echo 
echo "=================================================================================="
echo

echo "["`date +%H:%M:%S`"] Step 2: Running STAR alignment.."
cd $WDIR/fastqs 
mbc_prun_star.sh $WDIR $REFDIR $CONFIG $CPUS
echo 
echo "=================================================================================="
echo

echo "["`date +%H:%M:%S`"] Step 3: Making TDF and strand-specific bigWig files.." 
cd $WDIR/bams 
mbc_prun_coverage.sh $WDIR $REFDIR $CONFIG $CPUS
echo 
echo "=================================================================================="
echo

echo "["`date +%H:%M:%S`"] Step 4: Running featureCounts on all possible strand settings.."
cd $WDIR/bams 
mbc_prun_strand.sh $WDIR $REFDIR $CONFIG $CPUS
echo 
echo "=================================================================================="
echo

echo "["`date +%H:%M:%S`"] Step 5: Calculating strandedness and other statistics.."
cd $WDIR/fastqs 
mbc_prun_stats.sh $WDIR $CONFIG $CPUS
echo 
echo "=================================================================================="
echo

cd $WDIR/stats
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
cd $WDIR

echo 
echo "=================================================================================="
echo

echo "["`date +%H:%M:%S`"] Step 6: Running featureCounts on normal and extended annotation.."
cd $WDIR/bams 
mbc_prun_fcount.sh $WDIR $REFDIR $CONFIG $CPUS $STRAND
echo 
echo "=================================================================================="
echo

echo "["`date +%H:%M:%S`"] Step 7: Making final expression tables.."
cd $WDIR/featureCounts 
mbc_make_tables.sh $WDIR $REFDIR $CONFIG $CPUS
echo 
echo "=================================================================================="
echo


echo "["`date +%H:%M:%S`"] ALL YOUR BASE ARE BELONG TO US!!!"
