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
mv  Sample*  fastqs 
mv  P125109* study_strains
```

We should be all set to go. Let's run the reference preparation in `--simple` mode. Make sure you are in your main working directory - **ent_tutorial**. 

`prepare_bacterial_reference --simple . P125109 -p 8`

Your logs should look something like this: 

> ==> Initiating bacpipe reference preparation for SIMPLE SINGLE-STRAIN workflow!  
> ==> Following variables were set:  
>        WDIR: /pub37/alexp/data/rnaseq/other_Salmonella/6strains_SE  
>        TAG: P125109  
>        CPUS: 8  
> ==> Parallel jobs will be ran on 8 cores.  
> Directory /pub37/alexp/data/rnaseq/other_Salmonella/6strains_SE/study_strains/P125109 was not found and will be created.  
> ==> Running Prokka annotation to identify tRNA and rRNA loci  
> P125109.gff annotation file processed; found 22 gene entries without a locus tag, for which new locus tags were generated.  
> Parsed GFF stats: 4200 protein coding, 125 pseudogenes, 0 noncoding RNAs, 0 tRNAs, 22 rRNAs, 0 others.  
> ==> Files P125109.genome.fa and P125109.gene.gff successfully generated  
> ==> STAR aligner index P125109.STAR successfully generated  
> ==> rRNA/tRNA operon interval file P125109.rRNA.bed successfully created  
> All the generated files and indexes have been moved to /pub37/alexp/data/rnaseq/other_Salmonella/6strains_SE/study_strains/P125109.  
> Strain P125109: all done generating reference!  
> ==> Checking study strain P125109:  
> \- all necessary annotation files and directories for P125109 exist and are not empty.  
> ==> Writing simple config file, simple.cfg  
> ==> DONE generating simple single-strain reference!  

You will find that inside */study_strains* subdirectory, there is now a directory named P125109 that contains all the neccessary reference files (STAR index, Prokka annotation, rRNA intervals, etc). You will also find *simple.cfg* config file written in your main working directory. 

We're all ready to go now! Go to **ent_tutorial** and run bacpipe like this (adjust the number of cores used according to what's available on your system): 

`bacpipe . simple.cfg -p 64 | tee bacpipe.log`

This should follow the steps outlined in the main README, while printing logs to the screen. If you don't want to monitor the progress real-time, you can simply replace the command with 

`bacpipe . simple.cfg -p 64 &> bacpipe.log`

After this, you should obtain several folders with the processed RNA-seq data. 

* **FastQC** contains the output of fastqc intial read quality control;
* **bams** contains rRNA/tRNA-filtered indexed BAM files used for visualization, as well as STAR alignment logs;
* **tdfs_and_bws** contain TDF files (to be visualized in IGV), and straind-specific, normalized bigWig files (can be visualized in almost any genome browser); 
* **strand** and **featureCounts** contain the results of `featureCounts` quantification. The former contains tables for all possible strand-specificities, and the latter contains only the quantificaion for the correct strand; 
* **exp_tables** contains annotated and non-annotated expression tables in the form of raw counts and TPMs; 
* **stats** contains statistics calculated per individual sample. A summary of these statistics is given in the final log. 

## Multi-strain ("multi") workflow

Processing differs significantly for multiple strains. In order to demonstrate how it works, we will use a small part of [GSE119724](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119724) dataset, in which expression in *Salmonella* Typhimurium strains 4/74 and D23580 was compared under various *in vitro* conditions and inside macrophages. For simplicity, we will only use 6 samples of bacteria grown inside macrophages. 

First of all, let's choose a directory with plenty (~ 100 Gb) of space and enter it.

```bash 
mkdir stm_tutorial
cd stm_tutorial
```

After this, download the fastq files using `sra-tools` from NCBI. These tools could be installed with `conda install sra-tools`. If you get an error about unknown option, replace `--split-3` with `--split-e`.


Let's download the data: 

```bash 
fastq-dump --split-3 SRR7814119 & 
fastq-dump --split-3 SRR7814120 & 
fastq-dump --split-3 SRR7814121 & 
fastq-dump --split-3 SRR7814152 & 
fastq-dump --split-3 SRR7814153 & 
fastq-dump --split-3 SRR7814154 &
wait
```

After all of the downloading is done, let's give the files more informative names and archive them:
```bash 
mv SRR7814119.fastq D23_mac_rep1.fastq
mv SRR7814120.fastq D23_mac_rep2.fastq
mv SRR7814121.fastq D23_mac_rep3.fastq
mv SRR7814152.fastq 474_mac_rep1.fastq
mv SRR7814153.fastq 474_mac_rep2.fastq 
mv SRR7814154.fastq 474_mac_rep3.fastq
gzip *fastq
```

While the files are compressed, let's download the reference files we shall need for the multi-strain processing. First, create **study_strains** and **ref_strains** directories inside our main working directory, **stm_tutorial**. After this, let's download genome assemblies and GFF annotations (GenBank versions) of D23580 and 4/74. 

```bash
mkdir study_strains ref_strains
```
