#!/bin/bash 

TAG=$1
WDIR=$2
cd $WDIR/stats 

N_ALL=`grep "reads; of these:" $WDIR/bams/logs/$TAG.bowtie2_rrna.log | awk '{print $1}'`
N2=`grep "aligned 0 times"     $WDIR/bams/logs/$TAG.bowtie2_rrna.log | awk '{print $1}'`
N_RRNA=$((N_ALL-N2))
N3=`grep "aligned 0 times"     $WDIR/bams/logs/$TAG.bowtie2_trna.log | awk '{print $1}'`
N_TRNA=$((N2-N3))
N_UNMP=`grep "aligned 0 times" $WDIR/bams/logs/$TAG.bowtie2.log | awk '{print $1}'`
N_UNIQ=`grep "aligned exactly 1 time" $WDIR/bams/logs/$TAG.bowtie2.log | awk '{print $1}'`
N_MULT=`grep "aligned >1 times" $WDIR/bams/logs/$TAG.bowtie2.log | awk '{print $1}'`
N_ASSG=`grep "Successfully assigned reads" $WDIR/strand/logs/$TAG.fc.s0.log | awk '{print $6}'`

echo -e "$TAG\t$N_ALL\t$N_RRNA\t$N_TRNA\t$N_UNMP\t$N_MULT\t$N_UNIQ\t$N_ASSG" > $TAG.counts 

P_RRNA=`echo $N_ALL | awk -v v=$N_RRNA '{printf "%.3f\n",v*100/$0}'` 
P_TRNA=`echo $N_ALL | awk -v v=$N_TRNA '{printf "%.3f\n",v*100/$0}'` 
P_UNMP=`echo $N_ALL | awk -v v=$N_UNMP '{printf "%.3f\n",v*100/$0}'` 
P_UNIQ=`echo $N_ALL | awk -v v=$N_UNIQ '{printf "%.3f\n",v*100/$0}'` 
P_MULT=`echo $N_ALL | awk -v v=$N_MULT '{printf "%.3f\n",v*100/$0}'` 
P_ASSG=`echo $N_ALL | awk -v v=$N_ASSG '{printf "%.3f\n",v*100/$0}'` 

echo -e "$TAG\t$P_RRNA\t$P_TRNA\t$P_UNMP\t$P_MULT\t$P_UNIQ\t$P_ASSG" > $TAG.percent
 
R0=`grep "Successfully assigned reads" $WDIR/strand/logs/$TAG.fc.s0.log | awk '{print $6}'`
R1=`grep "Successfully assigned reads" $WDIR/strand/logs/$TAG.fc.s1.log | awk '{print $6}'`
R2=`grep "Successfully assigned reads" $WDIR/strand/logs/$TAG.fc.s2.log | awk '{print $6}'`
PCT=`echo "" | awk '{printf "%.3f\n",v1*100/v2}' v1=$R1 v2=$R0`

echo "Strandedness evaluation: R0 is $R0, R1 is $R1, R2 is $R2; percent first strand is $PCT" > $TAG.strand
