# bacpipe tutorial

Below is the tutorial illustrating the way of processing of bacterial RNA-seq experiments in `--simple` and `--multi` modes.

## Which mode to use? 

If you have RNA-seq experiments done in 1 strain, or strains that differ by a mutation/gene deletion/plasmid, use `--simple` with the reference that includes all elements of interest (e.g. if strains differ by a plasmid, make sure the plasmid is included in your reference sequence). If your strains differ substantially (e.g. multiple prophages and plasmids), use `--multi` mode. 

## Single-strain ("simple") workflow 
### Downloading the necessary files

In order to illustrate the processing of RNA-seq using `--simple` workflow, we will use the dataset [GSE50184](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50184), studying the influence of *misR* regulon in *Neisseria gonorrhoeae*. The dataset consists of 6 samples - 3 replicates of the wildtype and *misR* knockout. 

First of all, let's choose a directory with plenty of space and enter it. 

```bash 
mkdir Ngo_tutorial
cd Ngo_tutorial
```

After this, download the fastq files and rename them. You'll need working SRA tools (could be installed with `conda install sra-tools`). If you get an error about unknown option, replace `--split-3` with `--split-e`. 

```basi
for i in `seq 50 55`
do 
  fastq-dump --split-3 SRR9570$i & 
done
wait 
```
After all of the downloading is done, rename and gzip the fastq files:
```bash
gzip -c SRR957050.fastq > WT_rep1.fastq.gz & 
gzip -c SRR957051.fastq > WT_rep2.fastq.gz & 
gzip -c SRR957052.fastq > WT_rep3.fastq.gz & 
gzip -c SRR957053.fastq > misR_KO_rep1.fastq.gz & 
gzip -c SRR957054.fastq > misR_KO_rep2.fastq.gz & 
gzip -c SRR957055.fastq > misR_KO_rep3.fastq.gz & 
wait
rm SRR*fastq
```

While read download and processing is happening, let's download and rename the reference genome and its annotation. We will use strain [FA19](https://www.ncbi.nlm.nih.gov/assembly/GCA_001047225.1) GenBank assembly and annotation. 

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/047/225/GCA_001047225.1_ASM104722v1/GCA_001047225.1_ASM104722v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/047/225/GCA_001047225.1_ASM104722v1/GCA_001047225.1_ASM104722v1_genomic.gff.gz
gzip -d GCA*gz
mv GCA_001047225.1_ASM104722v1_genomic.fna FA19.fa
mv GCA_001047225.1_ASM104722v1_genomic.gff FA19.gff
```

Additionally, let's make a prophage interval file. As a crude estimate, we can use [PHASTER](http://phaster.ca/) results, that predict the presence of 4 prophage-like elements in FA19: 

  | Chromosome | Begin   |  End    | Name |
  |------------|---------|---------|------|
  | CP012026.1 | 191612  | 239575  | Ph1  |
  | CP012026.1 | 451221  | 474055  | Ph2  | 
  | CP012026.1 | 851813  | 880939  | Ph3  |
  | CP012026.1 | 1367842 | 1402902 | Ph4  |

Simply use any text editor and make the file called FA19.prophage.bed. Take the data from the table (without the header) and replace all spaces with a tab. Alternatively, download the file from [here](https://www.dropbox.com/s/luamz7n9m8wpe67/FA19.prophage.bed).

### Reference preparation

After all the files are downloaded and renamed, set up the directory structure and move files to appropriate sub-directories:
```bash 
mkdir fastqs      study_strains
mv    *.fastq.gz  fastqs 
mv    FA19.*      study_strains
```

We should be all set to go. Let's run the reference preparation in `--simple` mode. Make sure you are in your main working directory - **Ngo_tutorial**. 

`prepare_bacterial_reference --simple . FA19 -p 8`

Your logs should look something like this: 

> ==> Initiating bacpipe reference preparation for SIMPLE SINGLE-STRAIN workflow!
> ==> Following variables were set:
> 
>               WDIR: .
>                TAG: FA19
>               CPUS: 8
> 
> ==> Parallel jobs will be ran on 8 cores.
> Directory /pub37/alexp/data/rnaseq/for_bacpipe_paper/ngon_misR_GSE50184/study_strains/FA19 was not found and will be created.
> ==> Running Prokka annotation to identify tRNA and rRNA loci
> FA19.gff annotation file processed; found 0 gene entries without a locus tag, for which new locus tags were generated.
> Parsed GFF stats: 2220 protein coding, 0 pseudogenes, 1 noncoding RNAs, 56 tRNAs, 12 rRNAs, 0 others.
> ==> Files FA19.genome.fa and FA19.gene.gff successfully generated
> ==> STAR aligner index FA19.STAR successfully generated
> ==> rRNA/tRNA operon interval file FA19.rRNA.bed successfully created
> All the generated files and indexes have been moved to /pub37/alexp/data/rnaseq/for_bacpipe_paper/ngon_misR_GSE50184/study_strains/FA19.
> Strain FA19: all done generating reference!
> ==> Checking study strain FA19:
> - all necessary annotation files and directories for FA19 exist and are not empty.
> ==> Writing simple config file, simple.cfg
> ==> DONE generating simple single-strain reference!

Note the statistics of the GFF file (2220 protein coding, 0 pseudogenes, 1 noncoding RNAs, 56 tRNAs, 12 rRNAs, 0 others). If this is different from what you expect for your strain, examine the GFF file for the reason of such discrepancies.

You will find that inside **study_strains** subdirectory, there is now a directory named FA19 that contains all the neccessary reference files (STAR index, Prokka annotation, rRNA intervals, etc). You will also find *simple.cfg* config file written in your main working directory. 

### Running the main pipeline 

We're all ready to go now! Go to **Ngo_tutorial** and run bacpipe like this (adjust the number of cores used according to what's available on your system): 

`bacpipe . simple.cfg -p 64 | tee bacpipe.log`

This will write the *bacpipe.log* file while printing logs to the screen at the same time. If you don't want to monitor the progress real-time, you can simply replace the command with 

`bacpipe . simple.cfg -p 64 &> bacpipe.log`

After this, you should obtain several folders with the processed RNA-seq data. 

* **FastQC** contains the output of fastqc intial read quality control;
* **bams** contains rRNA/tRNA-filtered indexed BAM files used for visualization, as well as STAR alignment logs;
* **tdfs_and_bws** contain TDF files (to be visualized in IGV), and straind-specific, normalized bigWig files (can be visualized in almost any genome browser); 
* **strand** and **featureCounts** contain the results of `featureCounts` quantification. The former contains tables for all possible strand-specificities, and the latter contains only the quantificaion for the correct strand; 
* **exp_tables** contains annotated and non-annotated expression tables in the form of raw counts and TPMs; 
* **stats** contains statistics calculated per individual sample. A summary of these statistics is given in the final log. 

Strand specificity is reported in the log: 

> Average percent of reads matching the coding strand: 50.5977, lowest: 49.963, highest: 51.220

This particular RNA-seq experiment was not done in a strand-specific manner. 

Read number statistics are reported above: 

 | Sample             | Strain  |  All     |  rRNA    |  unmapped | multimap | 1_loc   | assigned |
 |--------------------|---------|----------|----------|-----------|----------|---------|----------|
 | misR_KO_rep1       | FA19    | 18692365 | 15551153 | 301948    | 117585   | 2721679 | 2740861  |
 | misR_KO_rep2       | FA19    | 7644453  | 5448011  | 285936    | 76984    | 1833522 | 1839844  |
 | misR_KO_rep3       | FA19    | 9183811  | 6997228  | 196893    | 78893    | 1910797 | 1930790  |
 | WT_rep1            | FA19    | 8293165  | 4336837  | 605150    | 141683   | 3209495 | 3316358  |
 | WT_rep2            | FA19    | 7730590  | 4930697  | 119528    | 119860   | 2560505 | 2584079  |
 | WT_rep3            | FA19    | 11372201 | 5711041  | 1346176   | 206595   | 4108389 | 4211370  |

This shows us the raw read numbers (read pairs for paired-end experiments). In this case we see that the experiments were not fully rRNA-depleted, with 8-18M initial reads and 2-4M reads assigned to annotated genes. Percentages of the initial reads are given below: 

 | Sample             |  Strain |   rRNA   | unmapped | multimap  | 1_loc  | assigned |
 |--------------------|---------|----------|----------|-----------|--------|----------|
 | misR_KO_rep1       | FA19    | 83.195   | 1.615    | 0.629     | 14.560 | 14.663   |
 | misR_KO_rep2       | FA19    | 71.268   | 3.740    | 1.007     | 23.985 | 24.068   |
 | misR_KO_rep3       | FA19    | 76.191   | 2.144    | 0.859     | 20.806 | 21.024   |
 | WT_rep1            | FA19    | 52.294   | 7.297    | 1.708     | 38.700 | 39.989   |
 | WT_rep2            | FA19    | 63.782   | 1.546    | 1.550     | 33.122 | 33.427   |
 | WT_rep3            | FA19    | 50.219   | 11.837   | 1.817     | 36.127 | 37.032   |

### Downstream processing

After you've successfully completed the processing, you can use the obtained expression tables to get PCA, clustering, and differential expression analysis. A quick and interactive processing can be done in [Phantasus](http://genome.ifmo.ru/phantasus-dev/). Taking annotated table of counts (*FA19.annotated.counts.tsv*), we can upload the table to the web tool, normalize it (Adjust, log2-transform, quantile normalize), and run PCA analysis (Tools -> Plots -> PCA plot). This should give us the following plot, indicating clear separation of mutants from WT samples: 

<img align="center" src="https://github.com/apredeus/multi-bacpipe/blob/master/img/PCA.png">

Quick differential expression using limma (Tools -> Differential expression -> limma) generates the following picture, identifying genes up- and down-regulated in misR mutants: 

<img align="center" src="https://github.com/apredeus/multi-bacpipe/blob/master/img/heatmap.png">

Overall, the pipeline streamlines the processing of simple, one-strain RNA-seq experiments.

## Multi-strain ("multi") workflow

Processing differs significantly for multiple strains. In order to demonstrate how it works, we will use a small part of [GSE119724](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119724) dataset, in which expression in *Salmonella* Typhimurium strains 4/74 and D23580 was compared under various *in vitro* conditions and inside macrophages. For simplicity, we will only use 6 samples of bacteria grown under *in vitro* conditions that induce SPI-2 gene expression (low pH, low phosphorus) - InSPI2. 

First of all, let's choose a directory with plenty (~ 50 Gb) of space and enter it.

```bash 
mkdir Stm_tutorial
cd Stm_tutorial
```

After this, download the fastq files using `sra-tools` from NCBI. These tools could be installed with `conda install sra-tools`. If you get an error about unknown option, replace `--split-3` with `--split-e`.


Let's download the data: 

```bash 
fastq-dump --split-3 SRR7814112 & 
fastq-dump --split-3 SRR7814113 & 
fastq-dump --split-3 SRR7814114 & 
fastq-dump --split-3 SRR7814145 &
fastq-dump --split-3 SRR7814146 &
fastq-dump --split-3 SRR7814147 &
wait
```

After all of the downloading is done, let's give the files more informative names and archive them:

```bash 
gzip -c SRR7814112.fastq > D23_InSPI2_rep2.fastq.gz &
gzip -c SRR7814113.fastq > D23_InSPI2_rep3.fastq.gz &
gzip -c SRR7814114.fastq > D23_InSPI2_rep4.fastq.gz &
gzip -c SRR7814145.fastq > 474_InSPI2_rep2.fastq.gz &
gzip -c SRR7814146.fastq > 474_InSPI2_rep3.fastq.gz &
gzip -c SRR7814147.fastq > 474_InSPI2_rep4.fastq.gz &
wait
gzip *fastq
```

While the files are compressed, let's download the reference files we shall need for the multi-strain processing. First, create **study_strains** and **ref_strains** directories inside our main working directory, **Stm_tutorial**. After this, let's download genome assemblies and GFF annotations (GenBank versions) of D23580 and 4/74. 

```bash
mkdir study_strains ref_strains fastqs

```
