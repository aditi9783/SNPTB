# SNPTB Tutorial

_A Bioinformatics Pipeline for Analyses of Mycobacterium tuberculosis Whole Genome Sequencing Data_


**Aditi Gupta**

Center for Emerging Pathogens, Rutgers University, Newark NJ 07103

## 0. Disclaimer and Prerequisites

All analyses are done with respect to the H37Rv reference genome of M. tuberculosis (NCBI accession ID: AL123456.3). The versions of the open-source software used in this pipeline are indicated and no guarantees are given that the older or newer versions of these software will perform as outlined in this tutorial. The pipeline has been tested on the Illumina HiSeq and MiSeq paired end data. The pipeline in its current form is thus incompatible with single-end sequencing data and data from non-Illumina sequencing technologies.

This tutorial is designed to teach biologists how to run a bioinformatics pipeline for analyzing whole genome sequencing datasets of Mycobacterium tuberculosis. Familiarity with the UNIX environment (in particular, traversing in and out of directories, copying & moving & deleting files, creating and deleting directories) and the ability to edit files in the UNIX environment is required. Even though the user will be executing python scripts, familiarity with python is NOT required.

If you use the SNPTB pipeline or any of its associated materials, please cite SNPTB as outlined in the title page. Your support encourages open-source science and is appreciated. 

## 1. Introduction

SNPTB is a collection of python scripts that are to be run in a specific order to get the following final outputs from your Illumina whole-genome sequencing data (in “fastq” format) of M. tuberculosis isolates:

•	The Single Nucleotide Polymorphisms (SNPs) in the sequencing data with respect to the reference genome H37Rv.
•	Annotation of the SNPs: name of the gene that has the mutation (or if the mutation is in the intergenic region), the functional effect of the mutation (synonymous or non-synonymous). 
•	Assessment of quality of read mapping: Total number of reads, % of reads mapped, average depth (number of reads that map to each position in the reference genome). 

This tutorial is designed with two things in mind:
1. The user has little computational experience. If not, I strongly recommend using the Googles to explore the wonderful world of UNIX. It is very straightforward and very little is required to use this tutorial.
2. The user has access to a high-performance computing cluster (usually available via your university, more on this later) that has the following open-source software pre-installed:

a) python/2.7.11   
b) bowtie2/2.2.6   
c) samtools/1.2   
d) bcftools/1.2  
e) java/1.8.0_66

### 1.1 Why do I need a high-performance computing cluster (HPCC)? Wait, what is a HPCC?

Each computer (including your laptop) has 2 things: memory, so that you can save a lot of data, and computing power, so that you can access data quickly and run programs. A computing cluster is, as the name suggests, a group of computers without the monitors (also called as servers). Using multiple servers together increases the memory as well as computing power of your “computer”. And when these servers are equipped for computationally intensive tasks such as storing and analyzing gigabases of sequencing data, the cluster is affectionately called a high-performance computing cluster (HPCC). 

Most universities provide support for computationally intensive research by establishing and maintaining a computational facility (under names such as “office for information technology” or “office for advanced computing research”) that purchases and maintains servers, installs commonly used scientific software on the servers, establish rules for fair access to the computing resources, troubleshoot when things go awry, and answer questions from researchers.

At Rutgers, anyone can request access to HPCC by filling out a form on the website of Office of Advanced Research Computing. Please familiarize yourself with the high performance computing center at your university by spending some time on their website.

At Rutgers, every new user is allotted 1 terabyte (=1000 GB) of space, although more space can be requested with justification. The Office for Advanced Research Computing (OARC) at Rutgers requires that programs that require longer to run (half an hour or more, for example) should be run via a job scheduler (more on this later). This not only ensures fair-usage of computing resources but also keeps your program running even if you log-out of your session with the HPCC.

Access to HPCC greatly simplifies genomics analyses because:
1. You don’t have to install bioinformatics software needed for analyses of next generation sequencing data.
2. If you have a lot of samples for analyses, the bioinformatics pipeline will be running for days. The job scheduler allows you to run your program “in the background” so that you can log-off after executing your program and check back later.
3. Sequencing data files are big: from hundreds of megabases to gigabases. Bioinformatics analyses of this data further generates several intermediate and final output files that themselves are tens of hundred of MB. Not only does the HPCC give you the space to store and do the analyses, all data is periodically backed up.

### 1.2 Connecting to your university’s HPCC.

In Mac, open Applications>Terminal to access the window where you would see a prompt (the $ sign). 

At the prompt, use your NetID and password to log in to the HPCC. For example, user with NetID “abcdefg” can access into the Perceval cluster of HPCC by typing:

```
$ ssh abcdefg@perceval.rutgers.edu
```

The user will be prompted for password, and if logging in for the very first time, will see a message asking for accepting “encryption key”. Type yes, and your first login is complete!

In Windows, you would need to install PUTTY to login to HPCC. 


### 1.3 Downloading SNPTB.

SNPTB is available on GitHub. It can be downloaded from https://github.com/aditi9783/SNPTB as a zip file or directly on your UNIX system by typing the following at the prompt:

```
$ git clone https://github.com/aditi9783/SNPTB.git
```

### 1.4 Is directory structure important?

Yes. Very. Because each sequencing project has multiple samples and the genomics analysis generates multiple files per sample, proper directory structure is essential to keep your data and results organized.

When you download SNPTB, you will see one directory called “SNPTB” that has the following sub-directories: “code”, “data”, "H37Rv", "docs", "output", "scheduler", "software", and "test". The directory “SNPTB” has a README file that describes the contents of each of the sub-directories in SNPTB. This tutorial itself can be found in the “docs”.

I would highly recommend that for each new sequencing project, you create a new directory under data (or any other location where you have space) and download your raw sequence data in that directory. It also pays dividends to follow a naming convention when creating these directories. For example, if I got sequencing data for my persistence project in October of 2017, then I would create a directory with the following name: 102017_AG_Persistence. If you expect to receive data from multiple projects, following a naming convention and having a helpful directory structure will keep you organized.

**Note:** Most scientific journals now require that you submit your raw sequenced data as well as some meta-analyses files to a public database. Thus, do not delete your raw sequencing data till you publish your results.

### 1.5. Downloading your sequencing data to your HPCC account from a hyperlink or from your local computer.

Now that you have downloaded SNPTB and created a new directory for your project, you can download your data in this new directory.

#### 1.5.1 Downloading data from hyperlinks
Most genomics cores send you hyperlinks to your sequencing data that you can directly download to your data directory on HPCC. Lets say your genomics core sent you an html page with links to your sequencing data that looks like this.

https://genome.med.nyu.edu/results/external/rutgers/2017-06-16/fastq/

Copy the link for the data file you want to use, and download the data directly to your HPCC account using the “wget” command (make sure you are in your project data directory before you download the data):
```
$ wget https://genome.med.nyu.edu/results/external/rutgers/2017-06-16/fastq/1_S1_L003_R1_001.fastq.gz . 
```
The dot above is a symbol for “current directory”. Thus, the data will be downloaded into the directory where you currently are (if not sure where you are, type 'pwd', short for _print working directory_, and you will see the full path of your current directory).

#### 1.5.2 Downloading data from your local computer

If you have sequencing data on your local computer, you can copy it to your project data directory on the HPCC using the 'scp' command. Scp syntax is same as that of the UNIX cp command.
```
$ scp <local directory> abcdefg@perceval.rutgers.edu:<your project data directory on Perceval>
```

Files can be copied from your directory on HPCC to your local computer by the following command:
```
$ scp abcdefg@perceval.rutgers.edu:<your project data directory on Perceval> <local directory>
```

### 1.6 Organizing your data.

Now that you have created a new project directory and downloaded the data, you should organize the data such that data for each sample is in a folder that is named after the sample. This is important because multiple files will be generated for each sample during the analyses and it is best to keep them organized.

First, lets move to the directory that has all the python scripts. If you are in "SNPTB", just type the following at the prompt.
```
$ cd code
```

Now, run the python script organize_data.py to create folders for each of your samples and copy raw sequence data into them. Provide full path to your data directory.

**NOTE:** To run a python script, you either have to be in the directory that contains the script (recommended) or specify directory path to the script.
```
$ python organize_data.py <full path to data directory>
```
If the code runs successfully, you will see the following message on the prompt:
```
$ python organize_data.py ../test
Data organization complete.
0:00:00.479501
```
The output shows that the code successfully ran without errors in 0.48 seconds.

## 2. Running the pipeline
### 2.1 From raw sequence data to SNP-calling: an overview

Once you have the paired-end raw sequenced reads in the fastq format (as is given by Illumina and some other sequencing technologies), you first need to remove low-quality reads or low-quality regions of good reads. This ‘quality control’ step improves the quality of the downstream analyses. The remaining high-quality reads are then mapped to the reference genome, and SNP calls are made from the mapped reads.

### 2.2 Quality control

During the library preparation of DNA samples for Illumina sequencing, small sequences called ‘adapters’ are ligated to the DNA fragments 

### 2.3 Mapping reads to the reference genome
```
$ bowtie2-build H37Rv.fa H37Rv
```
```
$ ls 
H37Rv.1.bt2  H37Rv.3.bt2  H37Rv.fa                 H37Rv_genes.txt   H37Rv_proteins.txt  H37Rv.rev.2.bt2
H37Rv.2.bt2  H37Rv.4.bt2  H37Rv_gene_features.txt  H37Rv_genpept.gp  H37Rv.rev.1.bt2
```
### 2.4 SNP calling

## 3. Data Analysis: SNP Annotation
### 3.1 Assessing data quality 
```
$ python get_avgdepth_genomecov.py /home/ag1349/SNPTB/test/
```

### 3.2 Understanding the SNP annotation output

### 3.3 Compare mutations between two samples 
```
python getpairedmut.py ../output/2017112916_test/SNPs_qscore_200.0.txt ../test/paired_ids.txt
```

## 4. Doing it all via the job scheduler
```
$ sbatch run_SNPTB_pipeline.sh
```

```
$ squeue -u <netid>
```

## 5. Getting ready to publish
### 5.1 What to include in the “Material and Methods” and how to cite the pipeline and its dependencies
### 5.2 Depositing raw and meta-data to databases ……………………………………



