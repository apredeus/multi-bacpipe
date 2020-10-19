#!/bin/bash 

TAG=$1
WDIR=$2
REFDIR=$3
SPECIES=$4
CPUS=$5
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

## deeptools does it faster and result is the same - and no need for an extra utility
bamCoverage --scaleFactor  $SCF -p $CPUS --filterRNAstrand reverse -b $TAG.bam -bs 1 -o $TAG.plus.bw  &>  $TAG.cov.log
bamCoverage --scaleFactor -$SCF -p $CPUS --filterRNAstrand forward -b $TAG.bam -bs 1 -o $TAG.minus.bw &>> $TAG.cov.log
igvtools count -z 5 -w 1 -e 0 $TAG.bam $TAG.tdf $CHROM &>> $TAG.cov.log  

echo
echo -e "command: bamCoverage --scaleFactor  $SCF -p $CPUS --filterRNAstrand reverse -b $TAG.bam -bs 1 -o $TAG.plus.bw  &> $TAG.cov.log"
echo -e "command: bamCoverage --scaleFactor -$SCF -p $CPUS --filterRNAstrand forward -b $TAG.bam -bs 1 -o $TAG.minus.bw &>> $TAG.cov.log"
echo -e "command: igvtools count -z 5 -w 1 -e 0 $TAG.bam $TAG.tdf $CHROM &>> $TAG.cov.log\n"
