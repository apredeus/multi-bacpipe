# bacpipe tutorial

Below is the tutorial illustrating the way of processing of bacterial RNA-seq experiments in `--simple` and `--multi` modes.

## Which mode to use? 

If you have RNA-seq experiments done in 1 strain, or strains that differ by a mutation/gene deletion/plasmid, use `--simple` with the reference that includes all elements of interest (e.g. if strains differ by a plasmid, make sure the plasmid is included in your reference sequence). If your strains differ substantially (e.g. multiple prophages and plasmids), use `--multi` mode. Both workflows require at least one well-annotated reference strain, and work best with complete genome assemblies (i.e. full assembled chromosomes, and not contigs). 

## Installation

Before starting this tutorial, install `bacpipe` and all dependencies as described [here](https://github.com/apredeus/multi-bacpipe#installation-and-dependencies). 

## Single-strain ("simple") workflow 
### Downloading the necessary files

In order to illustrate the processing of RNA-seq using `--simple` workflow, we will use the dataset [GSE50184](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50184), studying the influence of *misR* regulon in *Neisseria gonorrhoeae*. The dataset consists of 6 samples - 3 replicates each of the wildtype and *misR* knockout. 

First of all, let's choose a directory with plenty of space and enter it. 

```bash 
mkdir Ngo_tutorial
cd Ngo_tutorial
```
After this, download the fastq files using `fastq-dump` tool from NCBI `sra-tools` from NCBI (If you get an error about unknown option, replace `--split-3` with `--split-e` in the commands below):

```bash
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

```bash
CP012026.1	191612	239575	Ph1		
CP012026.1	451221	474055	Ph2			
CP012026.1	851813	880939	Ph3		
CP012026.1	1367842	1402902	Ph4		
```

Simply use any text editor, copy the table above, and save it as a text file named *FA19.prophage.bed*. Alternatively, download the file from [here](https://www.dropbox.com/s/luamz7n9m8wpe67/FA19.prophage.bed).

### Reference preparation

After all the files are downloaded and renamed, set up the directory structure and move files to appropriate sub-directories:
```bash 
mkdir fastqs study_strains
mv *.fastq.gz fastqs 
mv FA19.* study_strains
```

We should be all set to go. Let's run the reference preparation in `--simple` mode. Make sure you are in your main working directory - **Ngo_tutorial**:

`prepare_bacterial_reference --simple . FA19 -p 8`

Your logs should look something like this: 

>&nbsp;==> Initiating bacpipe reference preparation for SIMPLE SINGLE-STRAIN workflow!  
>&nbsp;==> Following variables were set:  
>&nbsp;  
>&nbsp;              WDIR: .  
>&nbsp;               TAG: FA19  
>&nbsp;              CPUS: 8  
>&nbsp;  
>&nbsp;==> Parallel jobs will be ran on 8 cores.  
>&nbsp;Directory /pub37/alexp/data/rnaseq/for_bacpipe_paper/ngon_misR_GSE50184/study_strains/FA19 was not found and will be created.  
>&nbsp;==> Running Prokka annotation to identify tRNA and rRNA loci  
>&nbsp;FA19.gff annotation file processed; found 0 gene entries without a locus tag, for which new locus tags were generated.  
>&nbsp;Parsed GFF stats: 2220 protein coding, 0 pseudogenes, 1 noncoding RNAs, 56 tRNAs, 12 rRNAs, 0 others.  
>&nbsp;==> Files FA19.genome.fa and FA19.gene.gff successfully generated  
>&nbsp;==> STAR aligner index FA19.STAR successfully generated  
>&nbsp;==> rRNA/tRNA operon interval file FA19.rRNA.bed successfully created  
>&nbsp;All the generated files and indexes have been moved to /pub37/alexp/data/rnaseq/for_bacpipe_paper/ngon_misR_GSE50184/study_strains/FA19.  
>&nbsp;Strain FA19: all done generating reference!  
>&nbsp;==> Checking study strain FA19:  
>&nbsp;- all necessary annotation files and directories for FA19 exist and are not empty.  
>&nbsp;==> Writing simple config file, simple.cfg  
>&nbsp;==> DONE generating simple single-strain reference!  

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

This shows us the raw read numbers (read pairs for paired-end experiments). In this case we see that the experiments were not fully rRNA-depleted, with 8-18 M initial reads and 2-4 M reads assigned to annotated genes. Percentages of the initial reads are given below: 

 | Sample             |  Strain |   rRNA   | unmapped | multimap  | 1_loc  | assigned |
 |--------------------|---------|----------|----------|-----------|--------|----------|
 | misR_KO_rep1       | FA19    | 83.195   | 1.615    | 0.629     | 14.560 | 14.663   |
 | misR_KO_rep2       | FA19    | 71.268   | 3.740    | 1.007     | 23.985 | 24.068   |
 | misR_KO_rep3       | FA19    | 76.191   | 2.144    | 0.859     | 20.806 | 21.024   |
 | WT_rep1            | FA19    | 52.294   | 7.297    | 1.708     | 38.700 | 39.989   |
 | WT_rep2            | FA19    | 63.782   | 1.546    | 1.550     | 33.122 | 33.427   |
 | WT_rep3            | FA19    | 50.219   | 11.837   | 1.817     | 36.127 | 37.032   |

The percentage of reads that don't map to rRNA/tRNA (multimap + 1_loc) is pretty close to "assigned" - meaning that the annotation does not miss any strongly expressed features. A larger difference would indicated that the annotation is missing some strongly expressed features (most commonly, ncRNA). 

### Downstream processing

After you've successfully completed the processing, you can use the obtained expression tables to get PCA, clustering, and differential expression analysis. A quick and interactive processing can be done in [Phantasus](http://genome.ifmo.ru/phantasus-dev/). Taking annotated table of counts (*FA19.annotated.counts.tsv*), we can upload the table to the web tool, normalize it (Adjust -> log2-transform -> quantile normalize), and run PCA analysis (Tools -> Plots -> PCA plot). This should give us the following plot, indicating clear separation of mutants from WT samples: 

<p align="center">
  <img src="https://github.com/apredeus/multi-bacpipe/blob/master/img/PCA.png?raw=true" alt="PCA plot for FA19 WT/misR KO RNA-seq"/>
</p>

Quick differential expression using limma (Tools -> Differential expression -> limma) generates the following picture, identifying genes up- and down-regulated in misR mutants: 

<p align="center">
  <img src="https://github.com/apredeus/multi-bacpipe/blob/master/img/heatmap.png?raw=true" alt="Limma differential expression analysis for FA19 WT/misR KO RNA-seq"/>
</p>

Overall, `bacpipe` streamlines the processing of simple, single-strain RNA-seq experiments.

## Multi-strain ("multi") workflow
### Downloading the necessary files

Processing differs significantly for multiple strains. In order to demonstrate how it works, we will use a small part of [GSE119724](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119724) dataset, in which expression in *Salmonella* Typhimurium strains 4/74 and D23580 (often shortened to D23) was compared under various *in vitro* conditions and inside macrophages. For simplicity, we will only use 6 samples of bacteria grown under *in vitro* conditions that induce SPI-2 gene expression (low pH, low phosphorus) - InSPI2. 

First of all, let's choose a directory with plenty (~ 50 Gb) of space and enter it.

```bash 
mkdir Stm_tutorial
cd Stm_tutorial
```

After this, download the fastq files using `fastq-dump` tool from NCBI `sra-tools` from NCBI (If you get an error about unknown option, replace `--split-3` with `--split-e` in the commands below):

```bash 
for i in `seq 12 14` `seq 45 47`
do
  fastq-dump --split-3 SRR78141$i & 
done
wait
```

After all of the downloading is done, let's give the files more informative names and compress them, deleting the original SRR files:

```bash 
gzip -c SRR7814112.fastq > D23_InSPI2_rep1.fastq.gz &
gzip -c SRR7814113.fastq > D23_InSPI2_rep2.fastq.gz &
gzip -c SRR7814114.fastq > D23_InSPI2_rep3.fastq.gz &
gzip -c SRR7814145.fastq > 474_InSPI2_rep1.fastq.gz &
gzip -c SRR7814146.fastq > 474_InSPI2_rep2.fastq.gz &
gzip -c SRR7814147.fastq > 474_InSPI2_rep3.fastq.gz &
wait
rm SRR*fastq
```

After this, let's download genome assemblies and GFF annotations (GenBank versions) for reference strains. We will use [old D23580](https://www.ncbi.nlm.nih.gov/assembly/GCA_000027025.1), [LT2](https://www.ncbi.nlm.nih.gov/assembly/GCA_000006945.2), and [SL1344](https://www.ncbi.nlm.nih.gov/assembly/GCA_000210855.2) (in order to be used as a reference strain, it needs to have a genomic GFF3 file available). For study strains, we will use the [new assembly of D23580](https://www.ncbi.nlm.nih.gov/assembly/GCA_900538085.1), and [4/74](https://www.ncbi.nlm.nih.gov/assembly/GCF_000188735.1/). It's also OK to use the same strain as a study and reference strain, but it's best to give them different strain tags. Here we will add "\_ref" to reference strain tags: 

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/027/025/GCA_000027025.1_ASM2702v1/GCA_000027025.1_ASM2702v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/027/025/GCA_000027025.1_ASM2702v1/GCA_000027025.1_ASM2702v1_genomic.gff.gz
gzip -d GCA*gz 
mv GCA_000027025.1_ASM2702v1_genomic.fna D23_ref.fa
mv GCA_000027025.1_ASM2702v1_genomic.gff D23_ref.gff
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/210/855/GCA_000210855.2_ASM21085v2/GCA_000210855.2_ASM21085v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/210/855/GCA_000210855.2_ASM21085v2/GCA_000210855.2_ASM21085v2_genomic.gff.gz
gzip -d GCA*gz
mv GCA_000210855.2_ASM21085v2_genomic.fna SL1344_ref.fa
mv GCA_000210855.2_ASM21085v2_genomic.gff SL1344_ref.gff
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/945/GCA_000006945.2_ASM694v2/GCA_000006945.2_ASM694v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/006/945/GCA_000006945.2_ASM694v2/GCA_000006945.2_ASM694v2_genomic.gff.gz
gzip -d GCA*gz
mv GCA_000006945.2_ASM694v2_genomic.fna LT2_ref.fa
mv GCA_000006945.2_ASM694v2_genomic.gff LT2_ref.gff 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/188/735/GCA_000188735.1_ASM18873v1/GCA_000188735.1_ASM18873v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/538/085/GCA_900538085.1_D23580_liv/GCA_900538085.1_D23580_liv_genomic.fna.gz
gzip -d GCA*gz
mv GCA_000188735.1_ASM18873v1_genomic.fna 474.fa 
mv GCA_900538085.1_D23580_liv_genomic.fna D23.fa 
```

Also download pre-made prophage BED files for [D23](https://www.dropbox.com/s/u75us6n1pgy3tka/D23.prophage.bed) and [474](https://www.dropbox.com/s/j8vhsi2qc74zntv/474.prophage.bed). These are the intervals generated by curated prophage annotation; a cruder, much faster estimate can be generated using [Phaster](http://phaster.ca/) online server. 

Finally, let's download a file of Salmonella-specific ncRNAs, small ORFs, and pseudogenes/selenocysteine genes - [St_ncNRA_sORF.fa](https://www.dropbox.com/s/q8tqsp7rz56wilt/St_ncNRA_sORF.fa).

### Reference preparation 

After all the necessary files have been downloaded, renamed, and compressed, let's place them in the appropriate folders: 

```bash
mkdir study_strains ref_strains fastqs
mv *fastq.gz fastqs
mv *_ref.fa *ref.gff ref_strains
mv 474.fa D23.fa *.bed study_strains
```
The extra reference file St_ncNRA_sORF.fa we've downloaded previously should remain in the main working directory. 

Now, let's make a config file. Using your favourite text editor, make a text file named *mult.cfg* and put the following (tab-separated) lines in it: 

```bash
D23_InSPI2_rep1	D23
D23_InSPI2_rep2	D23
D23_InSPI2_rep3	D23
474_InSPI2_rep1	474
474_InSPI2_rep2	474
474_InSPI2_rep3	474
Reference	LT2_ref
Reference	D23_ref
Reference	SL1344_ref
```

We should now be ready to make all the necessary reference files. Let's run the reference preparation in `--multi` mode. Make sure you are in your main working directory - **Stm_tutorial**: 

`prepare_bacterial_reference --multi . mult.cfg -p 16 -r St_ncNRA_sORF.fa`

The logs will indicate many performed operations and quality controls. Most important things to note are the following messages: 

>&nbsp;====> Merging blast-based annotation with Prokka-predicted features for strain 474:  
>&nbsp;...  
>&nbsp;==> Found 4867 protein-coding (CDS) and 280 non-coding RNA (misc_RNA/ncRNA) features.  
>&nbsp;====> Merging blast-based annotation with Prokka-predicted features for strain D23:  
>&nbsp;...  
>&nbsp;==> Found 4897 protein-coding (CDS) and 287 non-coding RNA (misc_RNA/ncRNA) features.  
>&nbsp;...  
>&nbsp;==> Directory /pub37/alexp/data/rnaseq/for_bacpipe_paper/jay_tut/ref_strains/D23_ref was not found and will be created.  
>&nbsp;...  
>&nbsp;Output GFF: 4446 CDS, 75 pseudogenes (as CDS), 0 tRNA, 0 rRNA, 0 ncRNA, and 0 others.  
>&nbsp;Total number of features in the output: 4521.  
>&nbsp;...  
>&nbsp;==> Directory /pub37/alexp/data/rnaseq/for_bacpipe_paper/jay_tut/ref_strains/LT2_ref was not found and will be created.  
>&nbsp;...  
>&nbsp;Output GFF: 4554 CDS, 39 pseudogenes (as CDS), 85 tRNA, 22 rRNA, 11 ncRNA, and 6 others.  
>&nbsp;Total number of features in the output: 4717.  
>&nbsp;...  
>&nbsp;==> Directory /pub37/alexp/data/rnaseq/for_bacpipe_paper/jay_tut/ref_strains/SL1344_ref was not found and will be created.  
>&nbsp;...  
>&nbsp;Output GFF: 4672 CDS, 72 pseudogenes (as CDS), 0 tRNA, 0 rRNA, 0 ncRNA, and 0 others.  
>&nbsp;Total number of features in the output: 4744.  
>&nbsp;...  
>&nbsp;==>  Performing Roary output reformatting; results will be written to roary/presence_absence_unix.csv  
>&nbsp;==> Generating a table of CDS orthologs.  
>&nbsp;...  
>&nbsp;==> DONE generating multi-strain reference!  

After this is complete, you should find the file *orthologs.tsv* in your main working directory. The file contains 5667 lines (+header), of which 
  *  5367 are CDS 
  *  290 are ncRNA
  *  7 are misc (pseudogenes etc)
  *  3 are repeat_region (CRISPR repeats identified by Prokka)

Of these, 
  *  4936 are located on (non-prophage part of) the chromosome
  *  398 are located in a prophage region; 
  *  333 are located on a plasmid. 

### Running the main pipeline 

We're all ready to go now! Go to **Stm_tutorial** and run bacpipe like this (adjust the number of cores used according to what's available on your system): 

`bacpipe . mult.cfg -p 64 -k LT2_ref | tee bacpipe.log`

This will write the *bacpipe.log* file while printing logs to the screen at the same time. If you don't want to monitor the progress real-time, you can simply replace the command with 

`bacpipe . mult.cfg -p 64 -k LT2_ref &> bacpipe.log`

After this, you should obtain several folders with the processed RNA-seq data. 

* **FastQC** contains the output of fastqc intial read quality control;
* **bams** contains rRNA/tRNA-filtered indexed BAM files used for visualization, as well as STAR alignment logs;
* **tdfs_and_bws** contain TDF files (to be visualized in IGV), and straind-specific, normalized bigWig files (can be visualized in almost any genome browser); 
* **strand** and **featureCounts** contain the results of `featureCounts` quantification. The former contains tables for all possible strand-specificities, and the latter contains only the quantificaion for the correct strand; 
* **exp_tables** contains annotated and non-annotated expression tables in the form of raw counts and TPMs; 
* **stats** contains statistics calculated per individual sample. A summary of these statistics is given in the final log. 

Strand specificity is reported in the log: 

> Average percent of reads matching the coding strand: 95.419, lowest: 95.134, highest: 95.740  
> The strandedness of your experiment was determined to be FR  

This RNA-seq experiment was done using a strand-specific protocol - particularly, the variety in which reads align to the same strand as the gene which gave rise to the read. Read number statistics are reported above:


 | Sample           | Strain  |         All  |     rRNA   | unmapped |  multimap |     1_loc  | assigned | 
 |------------------|---------|--------------|------------|----------|-----------|------------|----------|
 | 474_InSPI2_rep1  |    474  |     14333158 | 10478217   |  687949  |    17933  |   3149059  |  2948423 |
 | 474_InSPI2_rep2  |    474  |     15502120 | 11326243   |  498271  |    19149  |   3658457  |  3405101 |
 | 474_InSPI2_rep3  |    474  |     13405916 |  9592274   |  679457  |    16818  |   3117367  |  2895934 |
 | D23_InSPI2_rep1  |    D23  |      7785113 |  5747933   |  139763  |    11681  |   1885736  |  1773365 |
 | D23_InSPI2_rep2  |    D23  |      7057773 |  5296264   |  184865  |     9703  |   1566941  |  1473220 |
 | D23_InSPI2_rep3  |    D23  |      6534047 |  4859759   |   83644  |     9826  |   1580818  |  1478296 |

This shows us the raw read numbers (read pairs for paired-end experiments). In this case we see that the experiments were not fully rRNA-depleted, with 7-16 M initial reads and 1.5-3.4 M reads assigned to annotated genes. Percentages of the initial reads are given below: 

 | Sample            |  Strain  |  rRNA    | unmapped | multimap | 1_loc  | assigned |
 |-------------------|----------|----------|----------|----------|--------|----------|
 | 474_InSPI2_rep1   |    474   |   73.105 | 4.800    | 0.125    | 21.970 | 20.571   |
 | 474_InSPI2_rep2   |    474   |   73.063 | 3.214    | 0.124    | 23.600 | 21.965   |
 | 474_InSPI2_rep3   |    474   |   71.553 | 5.068    | 0.125    | 23.254 | 21.602   |
 | D23_InSPI2_rep1   |    D23   |   73.832 | 1.795    | 0.150    | 24.222 | 22.779   |
 | D23_InSPI2_rep2   |    D23   |   75.042 | 2.619    | 0.137    | 22.202 | 20.874   | 
 | D23_InSPI2_rep3   |    D23   |   74.376 | 1.280    | 0.150    | 24.194 | 22.625   |

The percentage of reads that don't map to rRNA/tRNA (multimap + 1_loc) is pretty close to "assigned" - meaning that the annotation does not miss any strongly expressed features. A larger difference would indicated that the annotation is missing some strongly expressed features (most commonly, ncRNA). 

### Downstream processing

After you've successfully completed the processing, you can use the obtained expression tables to get PCA, clustering, and differential expression analysis. A quick and interactive processing can be done in [Phantasus](http://genome.ifmo.ru/phantasus-dev/). Taking annotated table of counts (*Master_table.counts.tsv*), we can upload the table to the web tool, normalize it (Adjust -> log2-transform -> quantile normalize), and run PCA analysis (Tools -> Plots -> PCA plot). This should give us the following plot, indicating clear separation of mutants from WT samples: 

<p align="center">
  <img src="https://github.com/apredeus/multi-bacpipe/blob/master/img/PCA2.png?raw=true" alt="PCA plot, D23 vs 4/74 InSPI2 RNA-seq"/>
</p>

Quick differential expression using limma (Tools -> Differential expression -> limma) generates the following picture, identifying genes up- and down-regulated in misR mutants: 

<p align="center">
  <img src="https://github.com/apredeus/multi-bacpipe/blob/master/img/heatmap2.png?raw=true" alt="Limma differential expression analysis of D23 vs 4/74 InSPI2 RNA-seq"/>
</p>

Overall, `bacpipe` streamlines the processing of simple, single-strain RNA-seq experiments.


