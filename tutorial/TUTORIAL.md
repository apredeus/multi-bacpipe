# bacpipe tutorial

Below is the tutorial illustrating the way of processing of bacterial RNA-seq experiments in `--simple` and `--multi` modes.

## Single-strain ("simple") workflow 

In order to illustrate the processing of RNA-seq using `--simple` workflow, we will use the dataset [GSE46391](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46391), describing 1 replicate of 6 *Salmonella enterica* serovar Enteritidis strains (3 invasive and 3 non-invasive). 

First of all, let's choose a directory with plenty of space and enter it. 

```bash 
mkdir ent_tutorial
cd ent_tutorial
```

After this, download the fastq files and rename them. You'll need working SRA tools (could be installed with `conda install sra-tools`). If you get an error about unknown option, replace `--split-3` with `--split-e`. 

```bash 
fastq-dump --split-3 SRR1033744 & 
fastq-dump --split-3 SRR1033745 & 
fastq-dump --split-3 SRR1033746 & 
fastq-dump --split-3 SRR1033747 & 
fastq-dump --split-3 SRR1033748 & 
fastq-dump --split-3 SRR1033749 & 
fastq-dump --split-3 SRR1033750 & 
fastq-dump --split-3 SRR1033751 & 
fastq-dump --split-3 SRR1033752 & 
fastq-dump --split-3 SRR1033753 & 
fastq-dump --split-3 SRR1033754 & 
fastq-dump --split-3 SRR1033755 &
wait 
```
After all of the downloading is done, concatenate, rename, and gzip the 
```bash
cat SRR1033744.fastq SRR1033745.fastq | gzip > Sample_264G1LB.fastq.gz  &  
cat SRR1033746.fastq SRR1033747.fastq | gzip > Sample_264UK1LB.fastq.gz &  
cat SRR1033748.fastq SRR1033749.fastq | gzip > Sample_264BC8LB.fastq.gz &  
cat SRR1033750.fastq SRR1033751.fastq | gzip > Sample_264C19LB.fastq.gz &  
cat SRR1033752.fastq SRR1033753.fastq | gzip > Sample_264C45LB.fastq.gz &  
cat SRR1033754.fastq SRR1033755.fastq | gzip > Sample_264G45LB.fastq.gz &  
wait
rm SRR*fastq
```

While read download and processing is happening, let's download and rename the reference genome and its annotation. We will use strain [P125109](https://www.ncbi.nlm.nih.gov/assembly/GCA_000009505.1) GenBank assembly and annotation. 

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/009/505/GCA_000009505.1_ASM950v1/GCA_000009505.1_ASM950v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/009/505/GCA_000009505.1_ASM950v1/GCA_000009505.1_ASM950v1_genomic.gff.gz
gzip -d GCA*gz
mv GCA_000009505.1_ASM950v1_genomic.fna P125109.fa
mv GCA_000009505.1_ASM950v1_genomic.gff P125109.gff
```

After all this is done, set up the directory structure and move files to appropriate sub-directories:
```bash 
mkdir fastqs study_strains
mv Sample*fastq.gz fastqs 
mv P125109* study_strains
```

We should be all set to go. Let's run the reference preparation.   
