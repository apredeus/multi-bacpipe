#!/bin/bash 

## PIPELINE VERSION

cd strand/logs


for i in *.fc.s0.log 
do
  TAG=${i%%.fc.s0.log}
  echo "Calculating basic stats for tag $TAG.." 
  ../../strand_stat.sh $TAG 
done

mv *.rnastat ../../stats 

echo "ALL RNA-STAT CALCULATIONS ARE DONE!"
