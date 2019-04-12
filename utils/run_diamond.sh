#!/bin/bash 

## have to first make a reference, eg 
## diamond makedb -in /pub37/alexp/reference/IPG/IPG.06_jun_2018.protein.faa -db ipg 

FA=$1
TAG=${FA%%.faa}
DB=/pub37/alexp/reference/IPG/ipg

diamond blastp -d $DB -q $FA -f 6 qseqid sseqid stitle pident mismatch gapopen length qlen slen qstart qend sstart send evalue bitscore -o ${TAG}_vs_ipg.tsv 
