#!/bin/bash 

TAG=$1
REFDIR=$2
SPECIES=$3
CHROM=$REFDIR/Assemblies/${SPECIES}.chrom.sizes

bedtools genomecov -ibam $TAG.bam -bg -strand + > $TAG.plus.bedGraph
bedtools genomecov -ibam $TAG.bam -bg -strand - > $TAG.minus.bedGraph
bedGraphToBigWig $TAG.plus.bedGraph  $CHROM $TAG.plus.bw
bedGraphToBigWig $TAG.minus.bedGraph $CHROM $TAG.minus.bw
rm $TAG.plus.bedGraph $TAG.minus.bedGraph
igvtools count -z 5 -w 1 -e 0 $TAG.bam $TAG.tdf $CHROM
