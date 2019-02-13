#!/bin/bash 

TAG=$1
WDIR=$2
REFDIR=$3
SPECIES=$4
CHROM=$REFDIR/$SPECIES/${SPECIES}.chrom.sizes
cd $WDIR/bams 

## scaling is done to a coefficient of (L*10^6)/(n*l)
## where n - aligned read number, l - read length, L - average feature (gene) length (est 1kb). 
## we have lot's of softclips, so instead of n*l let's calculate actual number of bases available in the bedGraph (both strands)

BPS=`samtools view $TAG.bam | cut -f 6 | perl -ne 'm/(\d+)M/; print "$1\n";' | awk '{sum+=$1} END {print sum}'`
SCF=`echo $BPS | awk '{printf "%.3f\n",1000000000/$1}'`

## this is the scaling equivalent to TPM at every given nucleotide
## thus, the average height per gene should be equal to that gene's TPM

echo "Making bigWig and TDF files for sample $TAG.." 
echo "$TAG: scaling bigWig files by 10^9/n.aligned.bases, coefficient is $SCF"

bedtools genomecov -scale  $SCF -ibam $TAG.bam -bg -strand + > $TAG.plus.bedGraph 
bedtools genomecov -scale -$SCF -ibam $TAG.bam -bg -strand - > $TAG.minus.bedGraph 

bedGraphToBigWig $TAG.plus.bedGraph  $CHROM $TAG.plus.bw
bedGraphToBigWig $TAG.minus.bedGraph $CHROM $TAG.minus.bw
rm $TAG.plus.bedGraph $TAG.minus.bedGraph

echo -e "command: igvtools count -z 5 -w 1 -e 0 $TAG.bam $TAG.tdf $CHROM &> $TAG.tdf.log\n" > $TAG.tdf.log
igvtools count -z 5 -w 1 -e 0 $TAG.bam $TAG.tdf $CHROM &>> $TAG.tdf.log  
