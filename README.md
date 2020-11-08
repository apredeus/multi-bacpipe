# multi-bacpipe
A pipeline for uniform multi-strain RNA-seq processing.
<img align="right" width="281" height="330" src="https://github.com/apredeus/multi-bacpipe/blob/master/img/pay_or_i_play.jpg">

## Author
[Alexander Predeus](https://www.researchgate.net/profile/Alexander_Predeus), [Jay Hinton Laboratory](http://www.hintonlab.com/), [University of Liverpool](https://www.liverpool.ac.uk/)

(c) 2018-2020, GPL v3 license

## Why use bacpipe?
### 1. Transcriptomic comparison of multiple bacterial strains
Processing of multi-strain RNA-sequencing experiments is hard, because discrepancies between custom strain annotations introduce many false positive results. No previous approach allowed to analyze both orthologous and non-orthologous genes simultaneously. **Bacpipe** does just that: it processes RNA-seq data of several related strains in a uniform manner that minimizes the number of false positives. The results are generated in a form that's streamlined for downstream heatmap visualization/clustering analysis/PCA/etc in [Phantasus](https://artyomovlab.wustl.edu/phantasus/) or a similar tool, as well as differential gene expression analysis in DESeq2, limma-voom, or edgeR. 

What's even better is that you don't need to have individual annotation for each strain! Annotation files in GFF format could be generated with minimal effort after running the pipeline. All you need is several closely related, well-annotated strains available from NCBI. 
### 2. Informative RNA-seq statistics independent on annotation quirks or used library prep
Bacterial RNA-seq library prep could be done in several ways - without rRNA depletion, with partial depletion, or with depletion using Ribo-Zero or other kits. NCBI genome annotations often do not contain rRNA/tRNA genes, which makes alignment statistics very uninformative. 

This pipeline takes advantage of Prokka's uniform annotation of rRNA and tRNA genes, and calculates the alignment statitics of non-rRNA/tRNA reads, in addition to rRNA and tRNA. This gives a much fuller picture of the alignment and annotation. 
### 3. Easy processing and visualization for single-strain experiments
Even if you do not have a multi-strain RNA-seq experiment, the pipeline should allow for a quick and intuitive processing of RNA-seq experiment. Just choose `simple` workflow during the reference prep, and run the pipeline as usual. The pipeline automatically handles gzipped/non-gzipped/paired-end/single-end/mixed `fastq` files, and automatically estimates the strand-specificity of the RNA-seq library preparation method.  

## Input and output 

Pipeline requires the following input: 
* Genome sequences of each of the studied strains; 
* Genome sequences and GFF annotations of one or several reference strains;  
* Configuration file describing study and reference strains (for multi-strain mode only). 

When successfully applied, this should generate (for each strain):
* GFF annotation file for each study strain; 
* genomic BAM files for read-resolution visualization and analysis;
* TDF files for quick and scalable visualization in IGV;
* scaled, strand-specific bigWig files for visualization in JBrowse (see [doi:10.1128/mBio.01442-14](http://mbio.asm.org/content/5/4/e01442-14.full) for description of scaling); 
* raw read and TPM expression tables using [featureCounts](http://subread.sourceforge.net/) and accounting for multi-mapping reads; 

Additionally, the following things will be generated:
* master count/TPM tables with full list of coding and non-coding gene orthologs; 
* informative alignment statistics, including rRNA/tRNA content, number of multimappers, fraction of successfully assigned reads, etc.; 
* a single [MultiQC](http://multiqc.info/) quality control report.


## Installation and dependencies
Clone the pipeline scripts into your home directory and add them to $PATH variable in bash: 

```bash
cd ~
git clone https://github.com/apredeus/multi-bacpipe
echo "export PATH=~/multi-bacpipe:\$PATH" >> ~/.bashrc
```

To install the requirements, use [Bioconda](https://bioconda.github.io/). Below are the programs that need to be installed. Management of virtual environments in Bioconda is described [here](https://conda.io/docs/user-guide/tasks/manage-environments.html). First, create a separate virtual environment for all the tools:

```bash
conda update conda
conda create -n bacpipe
conda activate bacpipe
```  
After this, install all the necessary tools in the **bacpipe** environment:

```bash
conda install -c conda-forge -c bioconda -c defaults prokka=1.14.6
conda install -c bioconda roary=3.13.0
conda install deeptools=3.5.0
conda install fastqc=0.11.9
conda install star=2.7.6a
conda install samtools=1.11
conda install bedtools=2.29.2
conda install igvtools=2.3.93
conda install subread=2.0.1
conda install emboss=6.6.0
conda install diamond=2.0.4
```
## Tutorial

For the easiest introduction, please follow the bacpipe [tutorial](https://github.com/apredeus/multi-bacpipe/blob/master/tutorial/TUTORIAL.md). 

## Used terms
We use the following terms during the multi-strain RNA-seq analysis: 

* **study strains** - strains used for the experiment that is being processed;
* **reference strains** - related, well-annotated strains used as a reference for annotation; 
* **key strain** - one of reference strains which locus tags could be used when there's no gene name, e.g. STM IDs from LT2 for *Salmonella enterica* serovar Typhimurium;
* **working directory** - a directory where analysis is done; any writeable directory with enough space; 
* **strain tag** - an unique identifier describing a study or a reference strain; 

## Setting up the working directory
After you have chosen your working directory (**with enough space and where you have write permissions**), you need to create a sub-directory named **fastqs**. 
Place your archived RNA-seq *fastq* files here. Stick to the following naming rules: 

* Single-end sequencing: **Sample_ID.fastq.gz**; 
* Paired-end sequencing: **Sample_ID.R1.fastq.gz**, **Sample_ID.R2.fastq.gz**. 

If you have more than one file per sample (e.g. if it was sequenced twice), you must merge these files first.

## Reference preparation 
Please refer to [Tutorial](https://github.com/apredeus/multi-bacpipe/blob/master/tutorial/TUTORIAL.md) for an easy example of reference preparation in both `--simple` and `--multi` modes.
### 1. Simple mode 
In simple mode (one strain used for all experiments), you'll need to create a directory named **study_strains** and place two appropriately named files in it: 1) full genome assembly, *strain_tag.fa*; 2) GFF3 annotation, *strain_tag.gff*. Additionally, you can add a simple 4-column BED file that defines prophages and other regions of interest - e.g. pathogenicity islands or AR cassettes. For example, if you're using *Salmonella enterica* strain LT2, the files need to be called *LT2.fa* and *LT2.gff* (and, optionally, *LT2.prophage.bed*). After this, simply run the following command from your working directory (in which **fastqs** and **study_strains** are located):

`prepare_bacterial_reference --simple . <strain_tag>`

This will generate **strain_tag** subdirectory inside **study_strains**, as well as config file *simple.cfg*, which consists of simple, tab-separated table of samples and the reference strain: 

``` 
Sample_ID1	<strain_tag>
Sample_ID2	<strain_tag>
Sample_ID3	<strain_tag>
Sample_ID4	<strain_tag>
```

### 2. Multi-strain mode 
In multi-strain mode, more preparations are needed. First, create a directory named **ref_strains** and place two appropriately named files (full genome assembly, *strain_tag.fa*, and GFF3 annotation, *strain_tag.gff*) for each strain you want to use as a reference. Second, create a directory named **study_strains** and place 

#### Config file
In contrast to `--simple` mode, here you'll need to make the config file yourself. Take your favourite text editor and make a simple tab-separated text file. Each sample should be matched by appropriate study strain tag. Below all samples, list all available reference strains as "Reference" entries: 

```
Sample_ID1	<study_tag1>
Sample_ID2	<study_tag1>
Sample_ID3	<study_tag2>
Sample_ID4	<study_tag2>
Reference	<ref_tag1>
Reference	<ref_tag2>
Reference	<ref_tag3>
```
Optionally, you can also add prophage BED file to each of the study strains, and place the files (named *study_tag.prophage.bed*) into **study_strains** directory. 

#### Extra ncRNA/small CDS files 
In multi-strain mode, you can have an option to add a fasta file of non-coding RNAs and small ORFs identified in a closely related strain. The file should be a simple FASTA, formatted as follows: 

```
>STnc710.ncRNA
GTTTGGGTTCAAGTTAGCCCCCGTCAGGTTGCCAGGTTTATACCAGTCAACGTGCGGGGGTTTTCTCT
>SdsR.ncRNA
GCAAGGCGATTTAGCCTGCATTAATGCCAACTTTTAGCGCACGGCTCTCTCCCAAGAGCCATTTCCCTGGACCGAATACAGGAATCGTATTCGGTCTCTTTTT
>STnc3090.ncRNA
GCACGACGGTGGCCTTGACGACAACCTTCCTGGTGAATCCGGCTAAGGAGTAGAGTGGATTTCCCTTGGCCACCTCTGGCTTTGGCCTCTACTTTTCTCCAGGTCGTTTGCTGCCAAGACACCGCCGTGCGTT
>IsrB2.ncRNA
AAAACGCCCACCGAAGCGGGCGTGCCCTGTCCGGTCCAACCGACCAAAGCGAACCGGACCTAACAACCAGATATATCGGGGTGCTGTTAAGGCA
>mia-3.CDS
ATGACACTTTATTCTCTGAACGCACTTTGCAGACCTTTCCAGGATTAA
>mia-127.misc
ATGCTGGAGAATGTCATCATCTGA
```
Each feature can belong to one of the three categories: 1) ncRNA; 2) CDS; 3) misc. Categories differ in their downstream processing - e.g. ncRNAs overlapping CDS on the same strand are allowed, while a shorter CDS will be removed. "Misc" category is useful for pseudogenes, selenocysteine-coding genes, etc. 

## GFF file processing 
There is an impressive variety of GFF files that can be downloaded from GenBank or RefSeq. In order to make the analysis more robust, GFF files of each reference strain (or the study strain in case of `--simple` workflow) is pre-processed. However, pre-processing is done differently for `--multi` and `--simple` workflows.

### Simple workflow
In simple workflow, all features with defined `locus_tag` are preserved. Additionally, all non-rRNA/tRNA features are given a surrogate locus tag as well. A new GFF file is generated in which each feature is listed as a `gene` for `featureCounts` quantification. 
 
### Multi workflow
In multi workflow, only protein-coding and pseudogene features that have a defined `locus_tag` are used. Pseudogenes are converted to a smaller CDS using the longest annotated ORF. Resulting features are used to generate protein ortholog table using `Roary`. 

## One-command RNA-seq processing
After all the references are successfully created and placed appropriately, run 

`bacpipe <ref_dir> <config> [-p CPUs] [-k key_strain] &> bacpipe.log`

Bacpipe should to be ran in a writeable directory with **fastqs**, **study_strains**, and (in case of `--multi`), **ref_strains** directories, as well as config file. Key strain option can only be used in `--multi` mode. It defines a reference strain which locus tag could be used to replace a generic *group_number* ID in the final expression table. 

## Main pipeline steps

The following steps are performed during the pipeline execution: 
* FastQC is ran on all of the fastq files; 
* STAR is used to align the fastq files to the reference genome; 
* Unsorted bam alignments are filtered using rRNA BED file, sorted, and indexed; 
* tdf files are prepared for visualization in IGV; 
* Normalized bigWig (bw) files are prepared for vizualization in majority of other genomic browsers; 
* featureCounts is ran on genomic bam to evaluate the strandedness of the experiment; 
* strandedness and advanced alignment statistics are calculated; 
* featureCounts output is chosen based on the established strandedness; 
* appropriately formatted logs are generated; 
* multiqc is ran to summarize everything as a nicely formatted report. 
    
In the end you are expected to obtain a number of new directories: FastQC, bams, tdfs_and_bws, stats, strand, featureCounts, exp_tables. Each directory would contain the files generated by its namesake, as well as all appropriate logs. The executed commands with all of the versions and exact options are recorded in the master log. 
