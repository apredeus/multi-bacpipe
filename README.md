# multi-bacpipe
A pipeline for uniform multi-strain RNA-seq processing.
<img align="right" width="610" height="343" src="http://static.bnr.bg/sites/en/music/publishingimages/630/12-06-21-81065_2.jpg">

## Author
[Alexander Predeus](https://www.researchgate.net/profile/Alexander_Predeus), [Jay Hinton Laboratory](http://www.hintonlab.com/), [University of Liverpool](https://www.liverpool.ac.uk/)

(c) 2018, GPL v3 license

## Motivation
This is an upgrade to [bacpipe](https://github.com/apredeus/bacpipe) that is capable of performing uniform RNA-seq appropriate annotation of muliple strains and/or genomes, and then create master expression tables streamlined for a detailed strain comparison. 

When successfully applied, this should generate (for each strain):
* genomic bam files for read-resolution visualization and analysis;
* TDF files for visualization in IGV;
* scaled bigWig files for visualization in JBrowse (see [doi:10.1128/mBio.01442-14](http://mbio.asm.org/content/5/4/e01442-14.full) for description of scaling); 
* raw read and TPM expression tables - from [featureCounts](http://subread.sourceforge.net/), [rsem](https://deweylab.github.io/RSEM/), and [kallisto](https://pachterlab.github.io/kallisto/); 
* a single [MultiQC](http://multiqc.info/) quality control report.

This should also generate several master tables for all sorts of downstream analysis (clustering, PCA, differential expression, etc).

## Installation and requirements 
Clone the pipeline scripts into your home directory and add them to $PATH variable in bash: 

```bash
cd ~
git clone https://github.com/apredeus/multi-bacpipe
echo "export ~/multi-bacpipe:$PATH" >> .bashrc
```

To install the requirements, use [Bioconda](https://bioconda.github.io/). Below are the programs that need to be installed. In our experience, it heps to have roary, prokka, and rsem installed in their own virtual environments. Management of virtual environments in Bioconda is described [here](https://conda.io/docs/user-guide/tasks/manage-environments.html).

```bash
conda create -n roary 
conda create -n prokka 
conda create -n rsem 
``` 
After you've created virtual environments, install all the pre-requisites as follows:

```bash 
conda install -n roary roary 
conda install -n prokka prokka 
conda install -n rsem rsem

conda install fastqc
conda install bowtie2
conda install samtools
conda install bedtools
conda install igvtools
conda install rsem
conda install kallisto
conda install subread
```

You would need to clone and compile Jim Kent's utilites: 

```bash 
cd ~
git clone https://github.com/ENCODE-DCC/kentUtils
cd kentUtils
make
echo "export ~/kentUtils:$PATH" >> .bashrc
```

You also need to have Perl installed. Sorry. 

## Reference preparation
In order to start using the pipeline, you would need two things: a genomic *fasta* file, and genome annotation in *gff3* format. You can also use Prokka-style *gff3* file that has both annotation and sequence. The reference preparation script will try and guess which one are you using. 

It is very much recommended to develop a system of "tags" that you would use to identify references; for example, if you are processing data for P125109 strain of Salmonella enterica, and intend to use the assembly and annotation available from NCBI, rename the downloaded files to **P125109_ncbi.fa** and **P125109_ncbi.gff3**. After you set the reference directory, and run the reference-maker script, all of the reference indexes etc would be appropriately named and placed. For example, rsem reference would be in $REFDIR/rsem/P125109_ncbi_rsem, bowtie2 reference in $REFDIR/bowtie2/P125109_ncbi.\*.bt2, and so on. 

After you have procured the *fasta* and the *gff3* and selected a (writeable) reference directory, simply run 

`prepare_bacpipe_reference.sh <reference_dir> <tag> <name>.gff <name>.fa`

or 

`prepare_bacpipe_reference.sh <reference_dir> <tag> <name>.gff` 

for Prokka-style GFF file.

## One-command RNA-seq processing
After all the references are successfully created, simply run 

`bacpipe.sh <reference_dir> <tag> <CPUs>`

Bacpipe needs to be ran in a writeable directory with fastqs folder in it. 

Bacpipe:
* handles archived (.gz) and non-archived fastq files; 
* handles single-end and paired-end reads; 
* automatically detects strand-specificity of the experiment; 
* performs quantification according to the calculated parameters. 

The following steps are performed during the pipeline execution: 
* FastQC is ran on all of the fastq files; 
* bowtie2 is used to align the fastq files to the rRNA and tRNA reference to accurately estimate rRNA/tRNA content; 
* bowtie2 is used to align the fastq files to the genomic reference using --very-sensitive-local mode;
* sam alignments are filtered by quality (q10), sorted, converted to bam, and indexed; 
* tdf files are prepared for visualization in IGV; 
* bigWig (bw) files are prepared for vizualization in majority of other genomic browsers; 
* featureCounts is ran on genomic bam to evaluate the strandedness of the experiment; 
* strandedness and basic alignment statistics are calculated; 
* featureCounts output is chosen based correct settings of strandedness; 
* rsem is ran for EM-based quantification; 
* kallisto is ran to validate the RSEM results; 
* appropriately formatted logs are generated; 
* multiqc is ran to summarize everything as a nicely formatted report. 
    
In the end you are expected to obtain a number of new directories: FastQC, bams, tdfs_and_bws, RSEM, kallisto, stats, strand, featureCounts, exp_tables. Each directory would contain the files generated by its namesake, as well as all appropriate logs. The executed commands with all of the versions and exact options are recorded in the master log. 
    
    
