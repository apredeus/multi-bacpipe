#!/bin/bash 

## inspect expernal reference, making sure you don't have stops in the CDS sequences 
REF=$1
TAG=${REF%%.fa}

transeq -clean -table 11 -sequence $REF -outseq $TAG.translated.fa

seqtk seq -l 0 $TAG.translated.fa | grep -A1 "\.CDS" | grep -vP "\-\-" | perl -ne 'if (m/^>/) {chomp; s/^>//g; s/_1$//g; print "$_\t"} else {printf "%d\t%s",length($_)-1,$_}' > $TAG.CDS.tsv

N1=`grep ">" $REF | grep -c ncRNA`
N2=`grep ">" $REF | grep -c CDS`
N3=`grep ">" $REF | grep -c pseudogene`
NR=`grep -c ">" $REF`

echo "Reference contains $NR sequences: $N1 ncRNA, $N2 CDS, and $N3 pseudogenes"

if (( NR != N1 + N2 + N3 )) 
then
  echo "ERROR: Nucleotide reference file must contain only ncRNA, CDS, and pseudogenes. Your numbers don't add up." 
  exit 1
fi 

NP=`cat $TAG.CDS.tsv | wc -l`

echo "Number of features identified as CDS: $NP"
echo
echo "-----------------------------------------------------------------" 
echo 

KK=`cut -f 1 $TAG.CDS.tsv`
PP=`cut -f 3 $TAG.CDS.tsv`

a=( $KK ) 
b=( $PP ) 

for i in `seq 0 $((NP-1))`
do
  FIRST=`echo ${b[$i]} | awk '{print substr($0,1,1)}'`
  LAST=`echo  ${b[$i]} | awk '{print substr($0,length($0),1)}'`
  MIDX=`echo  ${b[$i]} | awk '{print substr($0,1,length($0)-1)}' | grep X`
  SIZE=`echo  ${b[$i]} | awk '{print length($0)}'`

  if [[ $LAST != "X" ]] 
  then 
    echo "WARNING: ${a[$i]} (sequence ${b[$i]}) - no stop codon."
  fi

  if [[ $FIRST != "M" && $FIRST != "V" && $FIRST != "L" ]]
  then
    echo "WARNING: ${a[$i]} (sequence ${b[$i]}) - unusual start codon."
  fi 

  if [[ $MIDX != "" ]]
  then 
    echo "WARNING: ${a[$i]} (sequence ${b[$i]}) - stop codon in the middle of the sequence."
  fi

  if (( SIZE < 12 )) 
  then
    echo "WARNING: ${a[$i]} (sequence ${b[$i]}) - CDS size is less than 36nt (12aa)."
  fi
done 
