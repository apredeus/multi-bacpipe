#!/bin/bash 

##   this script would generate everything in the working directory, which generally should be $WDIR/study_strains
##   no Prokka annotation of proteins - just find the CDS, Roary will do the rest. 

## next snippet is adapted from https://medium.com/@Drew_Stokes/bash-argument-parsing-54f3b81a6a8f

PARAMS=""
CPUS=""

while (( "$#" )); do
  case "$1" in
    -p|--cpus)
      CPUS=$2
      shift 2
      if [[ $CPUS == "" ]]
      then
        echo "ERROR: -p flag requires a non-empty argument (number of CPUs)!" 
        exit 1 
      fi
      echo "==> Invoking -p option: parallel jobs will be run on $CPUS cores."
      ;;  
    --) # end argument parsing
      shift
      break
      ;;  
    -*|--*=) # unsupported flags
      echo "ERROR: unsupported flag $1" >&2
      exit 1
      ;;  
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift
      ;;  
  esac
done
eval set -- "$PARAMS"

WDIR=$1
FA=$2
GFF=$3
PROPHAGE=$4
TAG=${FA%%.fa}

if [[ -d "$WDIR/study_strains/$TAG" ]]
then
  echo "Found $WDIR/study_strains/$TAG! Will add files to the existing directory."
  rm -rf $WDIR/study_strains/$TAG/*.STAR $WDIR/study_strains/$TAG/*.prokka 
else 
  echo "Directory $WDIR/study_strains/$TAG was not found and will be created." 
  mkdir $WDIR/study_strains/$TAG
fi

if [[ $CPUS == "" ]]
  then 
  echo "==> Parallel jobs will be ran on 4 cores (default)."
  CPUS=4
fi


######################################################

source activate prokka 
set -euo pipefail

## just like with prepare_strain_ref.sh, all this is happening in /study_strains 

cp $FA $TAG.genome.fa
samtools faidx $TAG.genome.fa
cut -f 1,2 $TAG.genome.fa.fai > $TAG.chrom.sizes 

## However, it does successfully find quite few leader peptides that are overlooked by Prokka otherwise. 
echo "Running Prokka annotation to unbiasedly identify tRNA and rRNA loci."

prokka --noanno --cpus $CPUS --outdir $TAG.prokka --prefix $TAG.prokka --locustag ${TAG%%_*} $TAG.genome.fa &> /dev/null 

## this will write $TAG.CDS.gff and $TAG.ncRNA.gff with all the correct attributes 
simple_gff_reannotation.pl $FA $GFF 
N_CDS=`grep -c -P "\tCDS\t" $TAG.CDS.gff`
N_NCR=`grep -c -P "\tncRNA\t" $TAG.ncRNA.gff`

echo
echo "==> Found $N_CDS protein-coding (CDS) and $N_NCR non-coding RNA (misc_RNA/ncRNA) features."
echo 

sed "s/\tCDS\t/\tgene\t/g"   $TAG.CDS.gff   >  $TAG.gene.gff
sed "s/\tncRNA\t/\tgene\t/g" $TAG.ncRNA.gff >> $TAG.gene.gff

echo "Files $TAG.genome.fa, $TAG.CDS.gff, $TAG.ncRNA.gff, and $TAG.gene.gff successfully generated."

## make STAR reference for small genome size 
mkdir ${TAG}.STAR 
STAR --runThreadN $CPUS --runMode genomeGenerate --genomeDir ${TAG}.STAR --genomeFastaFiles $TAG.genome.fa --genomeSAindexNbases 10 &> /dev/null
mv Log.out $TAG.star.log 

##make rRNA/tRNA interval file  
make_rrna_operon.pl $TAG.prokka/$TAG.prokka.gff $TAG.ncRNA.gff | sort -k1,1 -k2,2n | bedtools merge -i - > $TAG.rRNA.bed

## mv all to the ref dir 
mv $TAG.genome.fa $TAG.genome.fa.fai $TAG.chrom.sizes $WDIR/study_strains/$TAG
mv $TAG.gene.gff $TAG.CDS.gff $TAG.ncRNA.gff $WDIR/study_strains/$TAG
mv $TAG.prokka ${TAG}.STAR $TAG.star.log $WDIR/study_strains/$TAG
cp $PROPHAGE $WDIR/study_strains/$TAG/$TAG.prophage.bed
mv $TAG.rRNA.bed $WDIR/study_strains/$TAG

echo "All the generated files and indexes have been moved to $WDIR/study_strains/$TAG."
echo "Strain $TAG: all done generating reference!" 
