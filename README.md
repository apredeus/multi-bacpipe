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
echo "export PATH=~/multi-bacpipe:~/multi-bacpipe/utils:\$PATH" >> ~/.bashrc
```

To install the requirements, use [Bioconda](https://bioconda.github.io/). Below are the programs that need to be installed. Management of virtual environments in Bioconda is described [here](https://conda.io/docs/user-guide/tasks/manage-environments.html). First, create a separate virtual environment for all the tools:

```bash
conda update conda
conda create -n bacpipe
conda activate bacpipe
```  
After this, install all the necessary tools in the **bacpipe** environment:

```bash
conda install -c conda-forge -c bioconda -c defaults prokka
conda install -c bioconda roary
conda install deeptools
conda install fastqc
conda install star
conda install samtools
conda install bedtools
conda install igvtools
conda install subread
conda install emboss
conda install diamond
```
**/utils** contains Linux-compiled utilities necessary for file processing: [dos2unix](https://linux.die.net/man/1/dos2unix) converts Roary CSV output to Unix format, and [bedGraphToBigWig](https://github.com/ENCODE-DCC/kentUtils) is called upon when generating bigWig files. 

## Tutorial

For the easiest introduction, please follow the bacpipe [tutorial](https://github.com/apredeus/multi-bacpipe/blob/master/tutorial/TUTORIAL.md). 

## Setting up the working directory
After you have chosen your working directory (**with enough space and where you have write permissions**), you need to create a sub-directory named **fastqs**. 
Place your archived RNA-seq *fastq* files here. Stick to the following naming rules: 

* Single-end sequencing: **Sample_ID.fastq.gz**; 
* Paired-end sequencing: **Sample_ID.R1.fastq.gz**, **Sample_ID.R2.fastq.gz**. 

If you have more than one file per sample (e.g. if it was sequenced twice), you must merge these files first.

## Used terms
We use the following terms during the multi-strain RNA-seq analysis: 

* **study strains** - strains used for the experiment that is being processed;
* **reference strains** - related, well-annotated strains used as a reference for annotation; 
* **key strain** - one of reference strains which locus tags could be used when there's no gene name, e.g. STM IDs from LT2 for *Salmonella enterica* serovar Typhimurium;
* **working directory** - a directory where analysis is done; any writeable directory with enough space; 
* **strain tag** - an unique identifier describing a study or a reference strain; 
* **roary-friendly GFF** - GFF file formatted accordingly to criteria listed below, in [Reference strain requirements](https://github.com/apredeus/multi-bacpipe#reference-strain-requirements).

## Configuration file
Multi-strain processing is relying onto single-strain processing principles, described in [Bacpipe README](https://github.com/apredeus/bacpipe). In order to process multiple strains, you would need to make a simple tab-separated configuration file. The file should contain one tab-separated sample ID - strain ID pair per line. It should also include all of the strains that you want to use as a reference: 

```
Sample_ID1	Ecoli_O157_H7
Sample_ID2	Ecoli_O157_H7
Sample_ID3	Ecoli_O104_H4
Sample_ID4	Ecoli_O104_H4
Reference	MG1655
Reference	DH5a
Reference	DE3
```

See details on reference preparation below. 

## Reference strain requirements 
For each reference strain, you need to have a Roary-style GFF file placed into **ref_strains** subdirectory in your working directory. Expected name format is <tag>.roary.gff. Each file needs to satisfy the following criteria: 

* File consists of three parts, similar to Prokka output: 1) tab-separated GTF annotation; 2) ##FASTA delimiter; 3) strain's genomic fasta.
* Only coding sequence features (*CDS*) are present - no *gene*, *sequence_feature*, etc.; 
* Each coding seqeunce is annotated with a **unique** ID - locus tag for this strain; 
* If there's a common name for a gene associated with this CDS, it's given in Name= field; 
* If there is no common name reported, no Name is given. 

I provide a script called **ncbi_to_roary.pl** that converts NCBI-style genomic GFF and FNA (nucleotide FASTA) into required Roary GFF *most of the time*. Few NCBI annotations, like CT18, are truly horrific and require manual parsing, if you insist on using them. 

## Reference preparation 
### 1. Simple mode 
In simple mode, things are simple. 
### 2. Multi-strain mode 
In multi mode, things are not so simple.

## Prophage and rRNA/tRNA regions 
In order to provide additional information, all study strains need a BED file with putative prophage intervals. These can be generated by parsing the output of [PHASTER](http://phaster.ca/).  
Another important file is rRNA operon and tRNA BED file. It's used to accurately calculate the stats, especially the percentage of multi-mapping and un-assigned reads **after rRNA removal**. This file would be automatically generated from Prokka annotation when reference is prepared by `prepare_bacterial_reference`.

## One-command RNA-seq processing
After all the references are successfully created and placed appropriately, run 

`bacpipe <ref_dir> <config> [-p CPUs] [-k key_strain]`

Bacpipe needs to be ran in a writeable directory with fastqs, study_strains, ref_strains, and (in case of multi-strain mode) a roary folders in it. 

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
