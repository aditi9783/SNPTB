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
$ wget https://genome.med.nyu.edu/results/external/rutgers/2017-06-16/fastq/1_S1_L003_R1_001.fastq.gz 
```
The data will be downloaded into the directory where you currently are (if not sure where you are, type 'pwd', short for _print working directory_, and you will see the full path of your current directory).

#### 1.5.2 Downloading data from your local computer

If you have sequencing data on your local computer, you can copy it to your project data directory on the HPCC using the 'scp' command. Scp syntax is same as that of the UNIX cp command.
```
$ scp <local directory> abcdefg@perceval.rutgers.edu:<your project data directory on Perceval>
```

Files can be copied from your directory on HPCC to your local computer by the following command:
```
$ scp abcdefg@perceval.rutgers.edu:<your project data directory on Perceval> <local directory>
```

## 2. Running the pipeline
### 2.1 Organizing your data.

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
Lets say you have some data in the fastq.gz format in the "test" folder, and it looks liks this:
```
$ ls ../test/
A31_S8_R1_001.fastq.gz  A31_S8_R2_001.fastq.gz  D23_S15_R1_001.fastq.gz  D23_S15_R2_001.fastq.gz
```

To organize the data, type the following at the prompt (assuming you are in the "code" directory):
```
$ python organize_data.py ../test
Data organization complete.
0:00:00.479501
```
The output shows that the code successfully ran without errors in 0.48 seconds.

The data is now organized into subfolders that are named after the sample names:
```
$ ls ../test/
A31_S8  A31_S8_R1_001.fastq.gz  A31_S8_R2_001.fastq.gz  D23_S15  D23_S15_R1_001.fastq.gz  D23_S15_R2_001.fastq.gz

$ ls ../test/D23_S15
D23_S15_R1_001.fastq.gz  D23_S15_R2_001.fastq.gz

$ ls ../test/A31_S8
A31_S8_R1_001.fastq.gz  A31_S8_R2_001.fastq.gz
```
Two new subfolders are created: A31_S8 and D23_S15. The corresponding fastq.gz files are copied into sample folders and are not deleted from the /test/ folder because you want to save a copy of your data till you submit it to the public databases at the time of manuscript preparation.


### 2.2 From raw sequence data to SNP-calling: an overview

Once you have the paired-end raw sequenced reads in the fastq format (as is given by Illumina and some other sequencing technologies), you first need to remove low-quality reads or low-quality regions of good reads. This ‘quality control’ step improves the quality of the downstream analyses. The remaining high-quality reads are then mapped to the reference genome, and SNP calls are made from the mapped reads.

### 2.3 Quality control, mapping reads to the reference genome, and calling SNPs

#### 2.3.1 Quality Control
During the library preparation of DNA samples for Illumina sequencing, small sequences called ‘adapters’ are ligated to the DNA fragments. To remove these adapters, and to remove low quality reads, SNPTB uses the open-source software Trimmomatic. The quality of the filtered reads are assessed using FastQC.

#### 2.3.2 Read Mapping
Next step is to map the reads to the reference genome (_Mycobacterium tuberculosis_ H37Rv, NCBI Accession: AL123456.3) using Bowtie2. Before we can map the reads, Bowtie2 requires that you index the reference genome for faster lookup of sequence regions (indexes you see at the end of the books is a good analogy I think). This step has already been done, and the reference genome (.fa file) and the index files (.bt2 files) are stored in the "H37Rv" folder. 

This folder also contains additional files for the H37Rv (all downloaded from the NCBI) that are used in the SNP annotation step: the GenBank file, gene and protein sequence files, and so forth. The file 'gene_loci.txt' contains the gene start and end coordiantes. 
You can look at the contents of the "H37Rv" folder from the "code" folder location by typing the following at the prompt:

```
$ ls ../H37Rv/
genes_loci.txt  H37Rv.2.bt2  H37Rv.4.bt2  H37Rv.fa.fai             H37Rv_genes.txt   H37Rv_proteins.txt  H37Rv.rev.2.bt2
H37Rv.1.bt2     H37Rv.3.bt2  H37Rv.fa     H37Rv_gene_features.txt  H37Rv_genpept.gp  H37Rv.rev.1.bt2
```
#### 2.3.3 SNP Calling
After we map the reads to the reference genome, high confidence SNPs (probability that a SNP call is incorrect <1e-20) are identified using SAMTools and BCFtools. The output is stored in the Variant Calling Format (VCF file). 

A single script (qcmap.py) does all of these three steps: quality control, read mapping, and SNP calling. It _creates_ folders "qc" and "mapped" that contains the output files. Let's go over this step-by-step.

First, lets get access to the open-source software that SNPTB uses. Trimmomatic and FastQC are pre-packaged with SNPTB and thus you don't have to worry about it. But you need to get access to Bowtie2, SAMtools, and BCFtools. You have to do this before you run the qcmap.py script because qcmap.py _depends_ on these software. The benefit of working on the HPCC is that these are pre-installed (if not, you can request the HPCC staff and they will either install it for the entire cluster or will help you install it locally in your directory.

Because a lot of software are installed on the HPCC, you explicitly _load_ the ones that you need.

To see what softwares are installed on your HPCC cluster, type the following:
```
$ module avail

----------------------------------------------------------- /opt/sw/modulefiles/Core -----------------------------------------------------------
   ARACNE/20110228         bowtie2/2.2.6             gcc/4.9.4               intel_mkl/17.0.1        pgi/16.10           (D)
   HISAT2/2.0.4            bowtie2/2.2.9      (D)    gcc/5.3                 intel_mkl/17.0.2        python/2.7.11
   HISAT2/2.1.0     (D)    bwa/0.7.12                gcc/5.4          (D)    java/1.7.0_79           python/2.7.12
   MATLAB/R2017a           bwa/0.7.13         (D)    hdf5/1.8.16             java/1.8.0_66           python/3.5.0
   Mathematica/11.1        cuda/7.5                  intel/16.0.1            java/1.8.0_73           python/3.5.2        (D)
   OpenCV/2.3.1            cuda/8.0                  intel/16.0.3     (D)    java/1.8.0_121          samtools/0.1.19
   STAR/2.5.2a             cuda/9.0           (D)    intel/17.0.0            java/1.8.0_141   (D)    samtools/1.2
   Trinotate/2.0.2         cudnn/7.0.3               intel/17.0.1            modeller/9.16           samtools/1.3.1      (D)
   bamtools/2.4.0          cufflinks/2.2.1           intel/17.0.2            moe/2016.0802           trinityrnaseq/2.1.1
   bcftools/1.2            delly/0.7.6               intel/17.0.4            mvapich2/2.1
   bedtools2/2.25.0        gaussian/g03revE01        intel_mkl/16.0.1        mvapich2/2.2     (D)
   blast/2.6.0             gaussian/09revD01  (D)    intel_mkl/16.0.3 (D)    openmpi/2.1.1
   blat/35                 gcc/4.9.3                 intel_mkl/17.0.0        pgi/16.9

--------------------------------------------------- /opt/sw/admin/lmod/lmod/modulefiles/Core ---------------------------------------------------
   lmod/6.0.1    settarg/6.0.1

  Where:
   (D):  Default Module

Use "module spider" to find all possible modules.
Use "module keyword key1 key2 ..." to search for all possible modules matching any of the "keys".

```
To see which software are loaded in _your working environment_:
```
$ module list
No modules loaded
```
Now lets load the ones we would need:
```
$ module load python/2.7.11
$ module load bowtie2/2.2.6
$ module load samtools/1.2
$ module load bcftools/1.2
$ module load java/1.8.0_66
```

We are now ready to run the qcmap.py script. 
```
$ python qcmap.py <full path to data directory>
```
You will see output that looks like this:
```
Number of samples to be processed: 2

/home/ag1349/SNPTB/test/D23_S15
Trimming adapters and doing quality control ...
Running FastQC ...
Completed QC analysis. Files are in  /home/ag1349/SNPTB/test/D23_S15/qc
Mapping filtered reads to the reference ...
SNP calling ...
Completed read mapping and SNP calling. Files are in  /home/ag1349/SNPTB/test/D23_S15/mapped

/home/ag1349/SNPTB/test/A31_S8
Trimming adapters and doing quality control ...
Running FastQC ...
Completed QC analysis. Files are in  /home/ag1349/SNPTB/test/A31_S8/qc
Mapping filtered reads to the reference ...
SNP calling ...
Completed read mapping and SNP calling. Files are in  /home/ag1349/SNPTB/test/A31_S8/mapped
Quality control and read mapping complete.
0:42:03.958316
```
Thus, the script took 42 minutes to do quality control analysis, read mapping, and SNP calling. It also tells you where the quality control output is (the "qc" sub-folder within each sample's directory), and where the read mapping output is (the "mapped" sub-folder within each sample's directory). The SNP calls (the VCF files) are also saved in the "mapped" folder.

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
### 5.2 Depositing raw and meta-data to databases


