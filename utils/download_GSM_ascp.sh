#!/bin/bash 

## esearch does not return SRR anymore - SRA stopped including GSM IDs 
## so we'll just use the most up-to-date SRR/GSM dump 

TAB=/pub37/alexp/reference/SRA/only_SRR_GSM.21-oct-2020.tsv

GSM_LIST=$1
KK=`cat $GSM_LIST`
COUNT=1

for i in $KK
do

  echo "==> ["`date +%H:%M:%S`"] List $GSM_LIST, sample # $COUNT, sample $i"
  SRRS=`grep -w $i $TAB | cut -f 1`

  for j in $SRRS
  do
    echo "==> ["`date +%H:%M:%S`"] $i - downloading SRR $j"
    LEN=${#j}
    if (( $LEN < 10 )) 
    then
      MASK=`echo $j | perl -ne 'printf "%s\n",substr $_,0,6'`
      ascp -QT -l 300m -P33001 -i /pub37/alexp/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/$MASK/$j/ .
      if [[ $? == "1" ]]
      then
        echo "==> ["`date +%H:%M:%S`"] Failed to download $j via ascp; will try fastq-dump instead."
        fastq-dump --split-3 $j & 
      fi 
      mv $j/* .
      rm -rf $j 
    else 
      MASK=`echo $j | perl -ne 'printf "%s\n",substr $_,0,6'`
      SUB=`echo $j | awk '{print "00"substr($1,10,1)}'`
      ascp -QT -l 300m -P33001 -i /pub37/alexp/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/$MASK/$SUB/$j/ .
      if [[ $? == "1" ]]
      then
        echo "==> ["`date +%H:%M:%S`"] Failed to download $j via ascp; will try fastq-dump instead."
        fastq-dump --split-3 $j & 
      fi 
      mv $j/* .
      rm -rf $j
    fi 
  done
  COUNT=$((COUNT+1))
done

wait

echo "ALL DONE!"  
