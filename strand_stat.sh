#!/bin/bash 

TAG=$1
 
R0=`grep "Successfully assigned reads" $TAG.fc.s0.log | awk '{print $6}'`
R1=`grep "Successfully assigned reads" $TAG.fc.s1.log | awk '{print $6}'`
R2=`grep "Successfully assigned reads" $TAG.fc.s2.log | awk '{print $6}'`
PCT=`echo "" | awk '{printf "%.3f\n",v1*100/v2}' v1=$R1 v2=$R0`

echo "Strandedness evaluation: R0 is $R0, R1 is $R1, R2 is $R2; percent first strand is $PCT" > $TAG.rnastat
