#!/bin/bash 

TAG=$1
REFDIR=$2
SPECIES=$3
CHROM=$REFDIR/Assemblies/${SPECIES}.chrom.sizes

SCALE=`samtools view $TAG.bam | wc -l | awk '{printf "%.3f\n",1000000000/$1}'`
## this is the scaling equivalent to TPM at every given nucleotide
## thus, the average height per gene should be equal to that gene's TPM

echo "Making bigWig and TDF files for sample $TAG.." 
echo "$TAG: scaling bigWig files by 10^9/n.aligned.reads, coefficient is $SCALE"

bedtools genomecov -scale $SCALE -ibam $TAG.bam -bg -strand + > $TAG.plus.bedGraph
bedtools genomecov -scale $SCALE -ibam $TAG.bam -bg -strand - > $TAG.minus.bedGraph
bedGraphToBigWig $TAG.plus.bedGraph  $CHROM $TAG.plus.bw
bedGraphToBigWig $TAG.minus.bedGraph $CHROM $TAG.minus.bw
rm $TAG.plus.bedGraph $TAG.minus.bedGraph
igvtools count -z 5 -w 1 -e 0 $TAG.bam $TAG.tdf $CHROM &> $TAG.tdf.log  
