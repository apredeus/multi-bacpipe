#!/bin/bash 

## in 2020 no GSMs are available from SRA 
## thus we will rely on just another table 

GSM_NAME=$1
GSM_SRR=$2

if (( $# < 2 )) 
then
  echo "Usage: ./names_from_SRR.sh <gsm_to_name> <gsm_to_srr>" 
  exit 1
fi

N1=`cat $GSM_NAME | wc -l`
N2=`cat $GSM_SRR | wc -l`

if (( $N1 != $N2 )) 
then 
  echo "ERROR: Must have the same number of rows in $GSM_NAME and $GSM_SRR!"
  exit 1
fi 

KK=`cut -f1 $GSM_NAME`
PP=`cut -f2 $GSM_NAME`

a=( $KK )
b=( $PP ) 

for i in `seq 0 $((N2-1))`
do
  echo "Processing sample: records: ${a[$i]} ${b[$i]}"
  SRRS=`grep ${a[$i]} $GSM_SRR | cut -f 2`
  c=( $SRRS ) 
  if [[ ${c[1]} != "" ]] 
  then 
    echo "Multiple SRRS per GSM: ${a[$i]} :: $SRRS; will be concatenating the files!" 
    if [[ -s ${c[0]}.fastq ]] 
    then 
      echo "UNgzipped single-end files detected! Concatenating and archiving..."  
      cat `echo $SRRS | perl -ne '@t=split/\s+/; foreach $t (@t) {print "$t.fastq "}'` | gzip > ${b[$i]}.fastq.gz &
    elif [[ -s ${c[0]}.fastq.gz ]]
    then 
      echo "Gzipped single-end files detected! Concatenating and archiving..."
      zcat `echo $SRRS | perl -ne '@t=split/\s+/; foreach $t (@t) {print "$t.fastq.gz "}'` | gzip > ${b[$i]}.fastq.gz & 
    elif [[ -s ${c[0]}_1.fastq && -s ${c[0]}_2.fastq ]]
    then
      echo "UNgzipped paired-end files detected! Concatenating and archiving..." 
      cat `echo $SRRS | perl -ne '@t=split/\s+/; foreach $t (@t) {print "${t}_1.fastq "}'` | gzip > ${b[$i]}.R1.fastq.gz & 
      cat `echo $SRRS | perl -ne '@t=split/\s+/; foreach $t (@t) {print "${t}_2.fastq "}'` | gzip > ${b[$i]}.R2.fastq.gz & 
    elif [[ -s ${c[0]}_1.fastq.gz && -s ${c[0]}_2.fastq.gz ]]
    then
      echo "Gzipped paired-end files detected! Concatenating and archiving..." 
      zcat `echo $SRRS | perl -ne '@t=split/\s+/; foreach $t (@t) {print "${t}_1.fastq.gz "}'` | gzip > ${b[$i]}.R1.fastq.gz &
      zcat `echo $SRRS | perl -ne '@t=split/\s+/; foreach $t (@t) {print "${t}_2.fastq.gz "}'` | gzip > ${b[$i]}.R2.fastq.gz &
    fi
  else 
    echo "Single SRR per GSM: ${a[$i]} :: $SRRS; will be moving the files!" 
    if [[ -s ${c[0]}.fastq ]] 
    then 
      echo "UNgzipped single-end file detected! Archiving..."  
      gzip -c ${c[0]}.fastq > ${b[$i]}.fastq.gz & 
    elif [[ -s ${c[0]}.fastq.gz ]]
    then 
      echo "Gzipped single-end file detected! Moving..."
      mv ${c[0]}.fastq.gz ${b[$i]}.fastq.gz
    elif [[ -s ${c[0]}_1.fastq && -s ${c[0]}_2.fastq ]]
    then
      echo "UNgzipped paired-end files detected! Archiving..." 
      gzip -c ${c[0]}_1.fastq > ${b[$i]}.R1.fastq.gz & 
      gzip -c ${c[0]}_2.fastq > ${b[$i]}.R2.fastq.gz &
    elif [[ -s ${c[0]}_1.fastq.gz && -s ${c[0]}_2.fastq.gz ]]
    then
      echo "Gzipped paired-end files detected! Moving..." 
      mv ${c[0]}_1.fastq.gz ${b[$i]}.R1.fastq.gz
      mv ${c[0]}_2.fastq.gz ${b[$i]}.R2.fastq.gz
    fi
  fi
done   

wait
