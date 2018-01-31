#!/bin/bash 

## can be used for both single and paired-end
## archived FASTQ is assumed

TAG=$1
REFDIR=$2
SPECIES=$3
CPUS=$4

READS=""
REF=$REFDIR/bowtie2/${SPECIES}
RRNA=$REFDIR/bowtie2/${SPECIES}.rRNA
TRNA=$REFDIR/bowtie2/${SPECIES}.tRNA
WDIR=`pwd`

if [[ $TAG == "" || $SPECIES == "" || $CPUS == "" ]]
then 
  echo "ERROR: Please provide TAG, SPECIES handle, and # of CPUS! for individual bowtie2 runs!"
  exit 1
fi 

if [[ -e $TAG.fastq.gz ]]
then 
  echo "bowtie2: processing sample $TAG as single-ended, using bowtie2 index $REF."
  READS="-U $WDIR/$TAG.fastq.gz"
elif [[ -e $TAG.R1.fastq.gz && -e $TAG.R2.fastq.gz ]]
then
  echo "bowtie2: processing sample $TAG as paired-ended, using bowtie2 index $REF."
  READS="-1 $WDIR/$TAG.R1.fastq.gz -2 $WDIR/$TAG.R2.fastq.gz"
else
  echo "ERROR: The reqiured fastq.gz files were not found!" 
  exit 1
fi

if [[ ! -d $REFDIR/bowtie2 ]]
then 
  echo "ERROR: bowtie2 index $REF does not exist!"
  exit 1
fi

bowtie2 --very-sensitive-local -t -p $CPUS -x $RRNA $READS --un $TAG.norrna.fastq > /dev/null 2> $TAG.bowtie2_rrna.log 
bowtie2 --very-sensitive-local -t -p $CPUS -x $TRNA -U $TAG.norrna.fastq --un $TAG.nortrna.fastq > /dev/null 2> $TAG.bowtie2_trna.log 
bowtie2 --very-sensitive-local -t -p $CPUS -S $TAG.sam -x $REF -U $TAG.nortrna.fastq &> $TAG.bowtie2.log
samtools view -bS -q 10 $TAG.sam > $TAG.bam
samtools sort -@ $CPUS -T $TAG -o $TAG.sorted.bam $TAG.bam &> /dev/null 
mv $TAG.sorted.bam $TAG.bam
rm $TAG.sam $TAG.norrna.fastq $TAG.nortrna.fastq
samtools index $TAG.bam 

