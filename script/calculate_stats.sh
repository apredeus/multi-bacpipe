#!/bin/bash 

TAG=$1
WDIR=$2
SPECIES=$3
PAIRED="" 

cd $WDIR/fastqs 

if [[ -e $TAG.R1.fastq.gz && -e $TAG.R2.fastq.gz ]]
then
  PAIRED="true" 
elif [[ -e $TAG.fastq.gz ]] 
then 
  PAIRED="false" 
else 
  echo "ERROR: the experiment $TAG does not appear to be neither paired nor single!" 
  exit 1
fi 

echo "Calculating read alignment statistics for sample $TAG; library is PE: $PAIRED"

cd $WDIR/stats 

if [[ $PAIRED == "true" ]]
then
  N_ALL=`grep "Number of input reads" $WDIR/bams/logs/$TAG.star.log | awk '{print $6}'`
  UNI1=`grep "Uniquely mapped reads number" $WDIR/bams/logs/$TAG.star.log | awk '{print $6}'`
  MULT1=`grep "Number of reads mapped to multiple loci" $WDIR/bams/logs/$TAG.star.log | awk '{print $9}'`
  N_UNMP=$((N_ALL-UNI1-MULT1))

  ## note that BAM file already has rRNA removed

  UNI2=`samtools view $WDIR/bams/$TAG.bam  | grep    -P "\tNH:i:1\t" | cut -f 1 | sort | uniq | wc -l`
  MULT2=`samtools view $WDIR/bams/$TAG.bam | grep -v -P "\tNH:i:1\t" | cut -f 1 | sort | uniq | wc -l`
  N_RRNA=$((N_ALL-N_UNMP-UNI2-MULT2))
  N_UNIQ=$UNI2
  N_MULT=$MULT2

  ## from the non-strand-specific assignment - should be ok for now
  N_ASSG=`grep "Successfully assigned fragments" $WDIR/strand/logs/$TAG.fc.s0.log | awk '{print $6}'`
  
  echo -e "$TAG\t$SPECIES\t$N_ALL\t$N_RRNA\t$N_UNMP\t$N_MULT\t$N_UNIQ\t$N_ASSG" > $TAG.counts 
  
  P_RRNA=`echo $N_ALL | awk -v v=$N_RRNA '{printf "%.3f\n",v*100/$0}'` 
  P_UNMP=`echo $N_ALL | awk -v v=$N_UNMP '{printf "%.3f\n",v*100/$0}'` 
  P_UNIQ=`echo $N_ALL | awk -v v=$N_UNIQ '{printf "%.3f\n",v*100/$0}'` 
  P_MULT=`echo $N_ALL | awk -v v=$N_MULT '{printf "%.3f\n",v*100/$0}'` 
  P_ASSG=`echo $N_ALL | awk -v v=$N_ASSG '{printf "%.3f\n",v*100/$0}'` 
  
  echo -e "$TAG\t$SPECIES\t$P_RRNA\t$P_UNMP\t$P_MULT\t$P_UNIQ\t$P_ASSG" > $TAG.percent
   
  R0=`grep "Successfully assigned fragments" $WDIR/strand/logs/$TAG.fc.s0.log | awk '{print $6}'`
  R1=`grep "Successfully assigned fragments" $WDIR/strand/logs/$TAG.fc.s1.log | awk '{print $6}'`
  R2=`grep "Successfully assigned fragments" $WDIR/strand/logs/$TAG.fc.s2.log | awk '{print $6}'`
  PCT=`echo "" | awk '{printf "%.3f\n",v1*100/v2}' v1=$R1 v2=$R0`
  
  echo "Strandedness evaluation: R0 is $R0, R1 is $R1, R2 is $R2; percent first strand is $PCT" > $TAG.strand
else 
  N_ALL=`grep "Number of input reads" $WDIR/bams/logs/$TAG.star.log | awk '{print $6}'`
  UNI1=`grep "Uniquely mapped reads number" $WDIR/bams/logs/$TAG.star.log | awk '{print $6}'`
  MULT1=`grep "Number of reads mapped to multiple loci" $WDIR/bams/logs/$TAG.star.log | awk '{print $9}'`
  N_UNMP=$((N_ALL-UNI1-MULT1))
  
  ## note that BAM file already has rRNA removed

  UNI2=`samtools view $WDIR/bams/$TAG.bam | grep -c -P "\tNH:i:1\t"`
  MULT2=`samtools view $WDIR/bams/$TAG.bam | grep -v -P "\tNH:i:1\t" | cut -f 1 | sort | uniq | wc -l`
  N_RRNA=$((N_ALL-N_UNMP-UNI2-MULT2))
  N_UNIQ=$UNI2
  N_MULT=$MULT2
  N_ASSG=`grep "Successfully assigned reads" $WDIR/strand/logs/$TAG.fc.s0.log | awk '{print $6}'`
  
  echo -e "$TAG\t$SPECIES\t$N_ALL\t$N_RRNA\t$N_UNMP\t$N_MULT\t$N_UNIQ\t$N_ASSG" > $TAG.counts 
  
  P_RRNA=`echo $N_ALL | awk -v v=$N_RRNA '{printf "%.3f\n",v*100/$0}'` 
  P_UNMP=`echo $N_ALL | awk -v v=$N_UNMP '{printf "%.3f\n",v*100/$0}'` 
  P_UNIQ=`echo $N_ALL | awk -v v=$N_UNIQ '{printf "%.3f\n",v*100/$0}'` 
  P_MULT=`echo $N_ALL | awk -v v=$N_MULT '{printf "%.3f\n",v*100/$0}'` 
  P_ASSG=`echo $N_ALL | awk -v v=$N_ASSG '{printf "%.3f\n",v*100/$0}'` 
  
  echo -e "$TAG\t$SPECIES\t$P_RRNA\t$P_UNMP\t$P_MULT\t$P_UNIQ\t$P_ASSG" > $TAG.percent
   
  R0=`grep "Successfully assigned reads" $WDIR/strand/logs/$TAG.fc.s0.log | awk '{print $6}'`
  R1=`grep "Successfully assigned reads" $WDIR/strand/logs/$TAG.fc.s1.log | awk '{print $6}'`
  R2=`grep "Successfully assigned reads" $WDIR/strand/logs/$TAG.fc.s2.log | awk '{print $6}'`
  PCT=`echo "" | awk '{printf "%.3f\n",v1*100/v2}' v1=$R1 v2=$R0`
  
  echo "Strandedness evaluation: R0 is $R0, R1 is $R1, R2 is $R2; percent first strand is $PCT" > $TAG.strand
fi 
