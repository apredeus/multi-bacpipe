#!/bin/bash 

SDIR=$1
WDIR=$2
REF=$3
NEWREF=$4

cd $WDIR

## concatenate all study genomes, and blast ref against it 
echo
echo "Validating results by blasting external ncRNA/CDS/pseudogene reference to all study strain genomes.."
echo -e "\n--------------------------------------------------------------\n" 
cat $WDIR/study_strains/*.fa > All_study_genomes.fa 
samtools faidx All_study_genomes.fa
makeblastdb -dbtype nucl -in All_study_genomes.fa -out All_study_genomes.blast &> /dev/null 
blastn -query $REF -db All_study_genomes.blast -evalue 1 -outfmt 6 > All_study_genomes.ref_blast.out

## filter blast hits, convert them to GFF/BED, extract nucleotide sequences and convert to proteins
$SDIR/script/blast_to_gff.pl All_study_genomes.ref_blast.out $REF > All_study_genomes.blast.gff
grep -P "\tCDS\t" All_study_genomes.blast.gff | tr '=' '\t' | tr ';' '\t' | awk '{print $1"\t"$4-1"\t"$5"\t"$10"\t.\t"$7}' > All_study_genomes.CDS.bed
## this is here because bedtools 2.26 and 2.27 produce two different outputs with -name - this makes it version-independent (was I going insane at this point?)
bedtools getfasta -s -name -fi All_study_genomes.fa -fo - -bed All_study_genomes.CDS.bed | perl -ne 'if (m/^>([\w-]+)/) {print ">$1\n"} else {print}' > All_study_genomes.CDS.fa
transeq -clean -table 11 -sequence All_study_genomes.CDS.fa -outseq All_study_genomes.translated.fa &> /dev/null 

## find problematic CDS and change their type to pseudogene, keep new ref in $WDIR 
$SDIR/script/modify_reference_fasta.pl All_study_genomes.translated.fa $REF > $NEWREF

rm All_study_genomes.blast.n* All_study_genomes.fa All_study_genomes.ref_blast.out All_study_genomes.CDS.bed All_study_genomes.CDS.fa All_study_genomes.translated.fa All_study_genomes.blast.gff All_study_genomes.fa.fai 
echo -e "\n--------------------------------------------------------------\n" 
