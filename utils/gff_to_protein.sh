#!/bin/bash 

## get protein sequences using bacterial GFF3 and nucleotide fasta 
## requires EMBOSS tool transeq (conda install emboss) 

GFF=$1
FA=$2
TAG=${GFF%%.CDS.gff}

grep -v -P "#" $GFF | tr '=' '\t' | tr ';' '\t' | awk '{print $1"\t"$4-1"\t"$5"\t"$10"\t.\t"$7}' > $TAG.genes.bed 

bedtools getfasta -s -name -fi $FA -fo $TAG.genes.fa -bed $TAG.genes.bed

transeq -trim -clean -table 11 -sequence $TAG.genes.fa -outseq $TAG.protein.fa

perl -ne 's/\(\+\)_1//g; s/\(-\)_1//g; print' $TAG.protein.fa > $TAG.protein.fa.tmp 

mv $TAG.protein.fa.tmp $TAG.protein.fa
