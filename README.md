# multi-bacpipe
A pipeline for uniform multi-strain RNA-seq processing.
<img align="right" width="281" height="330" src="https://github.com/apredeus/multi-bacpipe/blob/master/img/pay_or_i_play.jpg">

## Author
[Alexander Predeus](https://www.researchgate.net/profile/Alexander_Predeus), [Jay Hinton Laboratory](http://www.hintonlab.com/), [University of Liverpool](https://www.liverpool.ac.uk/)

(c) 2018-2019, GPL v3 license

## Motivation
Processing of multi-strain RNA-sequencing experiments is hard, because discrepancies between custom strain annotations introduce many false positive results. No previous approach allowed to analyze both orthologous and non-orthologous genes simultaneously. **Bacpipe** does just that: it processes RNA-seq data of several related strains in a uniform manner that minimizes the number of false positives. The results are generated in a form that's streamlined for downstream heatmap visualization/clustering analysis/PCA/etc in [Phantasus](https://artyomovlab.wustl.edu/phantasus/) or a similar tool, as well as differential gene expression analysis in DESeq2, limma-voom, or edgeR. 

What's even better is that you don't need to have individual annotation for each strain! Annotation files in GFF format could be generated with minimal effort after running the pipeline. All you need is several closely related, well-annotated strains available from NCBI. 

Pipeline requires the following input: 
* Genome sequences of each of the studied strains; 
* Genome sequences and GFF annotations of several reference strains;  
* Configuration file describing study and reference strains. 

When successfully applied, this should generate (for each strain):
* GFF annotation file for each study strain; 
* genomic BAM files for read-resolution visualization and analysis;
* TDF files for quick and scalable visualization in IGV;
* scaled bigWig files for visualization in JBrowse (see [doi:10.1128/mBio.01442-14](http://mbio.asm.org/content/5/4/e01442-14.full) for description of scaling); 
* raw read and TPM expression tables using [featureCounts](http://subread.sourceforge.net/) and accounting for multi-mapping reads; 

Additionally, the following things will be generated:
* master count/TPM tables with full coding and non-coding gene orthology; 
* informative alignment statistics, including rRNA/tRNA content, number of multimappers, fraction of successfully assigned reads, etc.; 
* a single [MultiQC](http://multiqc.info/) quality control report.

## Installation and requirements 
Clone the pipeline scripts into your home directory and add them to $PATH variable in bash: 

```bash
cd ~
git clone https://github.com/apredeus/multi-bacpipe
echo "export PATH=~/multi-bacpipe:~/multi-bacpipe/utils:\$PATH" >> ~/.bashrc
```

To install the requirements, use [Bioconda](https://bioconda.github.io/). Below are the programs that need to be installed. Management of virtual environments in Bioconda is described [here](https://conda.io/docs/user-guide/tasks/manage-environments.html). First, create a separate virtual environment for all the tools:

```bash
conda create -n bacpipe
conda activate bacpipe
```  
After this, install all the necessary tools in the **bacpipe** environment:

```bash
conda install -c conda-forge -c bioconda -c defaults prokka
conda install -c bioconda roary
conda install fastqc
conda install star
conda install samtools
conda install bedtools
conda install igvtools
conda install subread
```
**/utils** contains Linux-compiled utilities necessary for file processing: [dos2unix](https://linux.die.net/man/1/dos2unix) converts Roary CSV output to Unix format, and [bedGraphToBigWig](https://github.com/ENCODE-DCC/kentUtils) is called upon when generating bigWig files. 

Make sure you have added **/utils** to your $PATH, or made these utilities available in any other way. Provided binaries should work on any 64-bit Unix machine. 

You also need to have Perl installed. Sorry. 

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
* **key strain** - one of reference strains which locus tags could be used instead of gene names, e.g. LT2 for *Salmonella enterica* serovar Typhimurium, or MG1655 for *Escherichia coli*;
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
For each reference strain, you need to have a Roary-style GFF file placed into **refstr** subdirectory in your working directory. Expected name format is <tag>.roary.gff. Each file needs to satisfy the following criteria: 

* File consists of three parts, similar to Prokka output: 1) tab-separated GTF annotation; 2) ##FASTA delimiter; 3) strain's genomic fasta.
* Only coding sequence features (*CDS*) are present - no *gene*, *sequence_feature*, etc.; 
* Each coding seqeunce is annotated with a **unique** ID - locus tag for this strain; 
* If there's a common name for a gene associated with this CDS, it's given in Name= field; 
* If there is no common name reported, no Name is given. 

I provide a script called **ncbi_to_roary.pl** that converts NCBI-style genomic GFF and FNA (nucleotide FASTA) into required Roary GFF *most of the time*. Few NCBI annotations, like CT18, are truly horrific and require manual parsing, if you insist on using them. 

## Two-step reference preparation
Reference preparation includes two steps: 
* Prepare reference for all strains you are using in your RNA-seq using **prepare_strain_ref.sh** (in the config example above, these are Ecoli_O157_H7 and Ecoli_O104_H4). **All you need to have for your study strains is genome assembly, prophage BED file, and rRNA operon BED file.** Annotations are not required; if a good annotation of a strain exists, it could be used as a reference strain, using a different name. 

Using the example config above, the command would have to be applied twice, once for each of the study strains (Ecoli_O157_H7, Ecoli_O104_H4):

```bash
prepare_strain_ref.sh <ref_directory> Ecoli_O157_H7 Ecoli_O157_H7.fa <prophage_bed> <rRNA_bed>
prepare_strain_ref.sh <ref_directory> Ecoli_O104_H4 Ecoli_O104_H4.fa <prophage_bed> <rRNA_bed>
```

* Prepare combined reference using select reference strains (in the config example above, tags are MG1655, DH5a, and DE3) by running **prepare_multiref.sh**. 

```bash
prepare_multiref.sh <working_directory> <ref_directory> <config>
```

## Prophages and rRNA operons
In order to provide additional information, all study strains would need two additional BED files:
* Putative prophage intervals - these can be generated by parsing the output of [PHASTER](http://phaster.ca/);  
* rRNA operon and tRNA BED file. This is very important in order to accurately calculate the stats, especially the percentage of multi-mapping and un-assigned reads **after rRNA removal**. This file would be automatically generated when overall reference is prepared.

## One-command RNA-seq processing
After all the references are successfully created and placed appropriately, run 

`bacpipe <ref_dir> <config> [-p CPUs] [-k key_strain]`

Bacpipe needs to be ran in a writeable directory with fastqs, study_strains, ref_strains, and roary folders in it (see above for reference preparation). 

Bacpipe:
* handles archived (.gz) and non-archived fastq files; 
* handles single-end and paired-end reads; 
* automatically detects strand-specificity of the experiment; 
* performs quantification according to the calculated parameters. 

The following steps are performed during the pipeline execution: 
* FastQC is ran on all of the fastq files; 
* STAR is used to align the fastq files to the rRNA and tRNA reference to accurately estimate rRNA/tRNA content; 
* Unsorted bam alignments are filtered using rRNA BED file, sorted, and indexed; 
* tdf files are prepared for visualization in IGV; 
* bigWig (bw) files are prepared for vizualization in majority of other genomic browsers; 
* featureCounts is ran on genomic bam to evaluate the strandedness of the experiment; 
* strandedness and basic alignment statistics are calculated; 
* featureCounts output is chosen based correct settings of strandedness; 
* appropriately formatted logs are generated; 
* multiqc is ran to summarize everything as a nicely formatted report. 
    
In the end you are expected to obtain a number of new directories: FastQC, bams, tdfs_and_bws, stats, strand, featureCounts, exp_tables. Each directory would contain the files generated by its namesake, as well as all appropriate logs. The executed commands with all of the versions and exact options are recorded in the master log. 
