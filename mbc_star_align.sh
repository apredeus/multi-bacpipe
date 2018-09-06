#!/bin/bash 

## can be used for both single and paired-end
## archived FASTQ is assumed

TAG=$1
WDIR=$2
REFDIR=$3
SPECIES=$4
CPUS=$5

READS=""
PAIRED=""
REF=$REFDIR/${SPECIES}/${SPECIES}.STAR
RRNA=$REFDIR/${SPECIES}/${SPECIES}.rRNA.bed ## use whole operons here
FQDIR=$WDIR/fastqs 

if [[ $TAG == "" || $SPECIES == "" || $CPUS == "" ]]
then
  echo "Usage: bc_star_align.sh <tag> <wdir> <refdir> <species> <cpus>"  
  echo "ERROR: Please provide sample ID (TAG), SPECIES handle, and # of CPUS for individual STAR runs."
  exit 1
fi 

cd $FQDIR
mkdir ${TAG}_star
cd ${TAG}_star

if [[ -e $FQDIR/$TAG.fastq.gz ]]
then 
  echo "STAR: processing sample $TAG as single-ended."
  READS="$FQDIR/$TAG.fastq.gz"
  PAIRED="false"
elif [[ -e $FQDIR/$TAG.R1.fastq.gz && -e $FQDIR/$TAG.R2.fastq.gz ]]
then
  echo "STAR: processing sample $TAG as paired-ended."
  READS="$FQDIR/$TAG.R1.fastq.gz $FQDIR/$TAG.R2.fastq.gz"
  PAIRED="true"
else
  echo "ERROR: The reqiured fastq.gz files were not found!" 
  exit 1
fi

if [[ ! -d $REF || ! -e $RRNA ]]
then 
  echo "ERROR: either STAR genome index $REF or rRNA interval file $RRNA does not exist!"
  exit 1
fi

echo -e "command: STAR --genomeDir $REF --readFilesIn $READS --alignIntronMin 20 --alignIntronMax 19 --readFilesCommand zcat --runThreadN $CPUS --outFilterMultimapNmax 20 --outReadsUnmapped Fastx"
STAR --genomeDir $REF --readFilesIn $READS --alignIntronMin 20 --alignIntronMax 19 --readFilesCommand zcat --runThreadN $CPUS --outFilterMultimapNmax 20 --outSAMtype BAM Unsorted &> /dev/null 

## only take the reads that don't overlap rRNA operons 
bedtools intersect -nonamecheck -v -a Aligned.out.bam -b $RRNA | samtools sort -@ $CPUS - > $TAG.bam 
samtools index $TAG.bam
rm Aligned.out.bam
mv Log.final.out $TAG.star.log 
