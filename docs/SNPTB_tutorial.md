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

<your_dir_path>/SNPTB/test/D23_S15
Trimming adapters and doing quality control ...
Running FastQC ...
Completed QC analysis. Files are in  <your_dir_path>/SNPTB/test/D23_S15/qc
Mapping filtered reads to the reference ...
[mpileup] 1 samples in 1 input files
<mpileup> Set max per-file depth to 8000
SNP calling ...
[mpileup] 1 samples in 1 input files
<mpileup> Set max per-file depth to 8000
SNP calling ...
Completed read mapping and SNP calling. Files are in  <your_dir_path>/SNPTB/test/D23_S15/mapped

<your_dir_path>/SNPTB/test/A31_S8
Trimming adapters and doing quality control ...
Running FastQC ...
Completed QC analysis. Files are in  <your_dir_path>/test/A31_S8/qc
Mapping filtered reads to the reference ...
[mpileup] 1 samples in 1 input files
<mpileup> Set max per-file depth to 8000
SNP calling ...
[mpileup] 1 samples in 1 input files
<mpileup> Set max per-file depth to 8000
SNP calling ...
Completed read mapping and SNP calling. Files are in  <your_dir_path>/test/A31_S8/mapped
Quality control and read mapping complete.
0:42:03.958316
```
Thus, the script took 42 minutes to do quality control analysis, read mapping, and SNP calling. It also tells you where the quality control output is (the "qc" sub-folder within each sample's directory), and where the read mapping output is (the "mapped" sub-folder within each sample's directory). The SNP calls (the VCF files) are also saved in the "mapped" folder.

Lets check out the directory structre for sample D23_S15:

```
$ ls ../test/D23_S15
mapped  qc
$ ls ../test/D23_S15/qc/
fastqc.run  s1_pe_fastqc.html  s1_se_fastqc.html  s2_pe_fastqc.html  s2_se_fastqc.html
qc.run      s1_pe_fastqc.zip   s1_se_fastqc.zip   s2_pe_fastqc.zip   s2_se_fastqc.zip
$ ls ../test/D23_S15/mapped/
aln.sorted.bam.mpileup  D23_S15.aln.sorted.bam  D23_S15.aln.sorted.bam.bai  D23_S15.vcf  map.run
```

Note that the after the analysis is done, the script deletes the raw data fastq.gz files that were stored in the sample directory are deleted (but they are still there in the data directory!). The script also deletes a few other temporarily generated files (such as the filtered reads after quality control). This is important because usually you will have a lot of sequencing data and the script _will_ fail if you run out of space in your HPCC account.

Thus, only two folders remain in each sample directory.

The "qc" folder contains the output of Trimmomatic and FastQC. File 'qc.run' tells you how many reads were there in your data files, and how many low-quality reads were dropped. It also contains some .html files that you can save to you local computer and see the quality of your filtered data, an example can be found here: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.

The "mapped" folder contains the output of Bowtie2 (the mapped reads) and that of SAMtools and BCFtools (VCF file - SNP calls). The file 'map.run' describes how many reads were successfully mapped and so forth:

```
$ more ../test/D23_S15/mapped/map.run 
718668 reads; of these:
  688540 (95.81%) were paired; of these:
    59007 (8.57%) aligned concordantly 0 times
    619168 (89.92%) aligned concordantly exactly 1 time
    10365 (1.51%) aligned concordantly >1 times
    ----
    59007 pairs aligned concordantly 0 times; of these:
      45017 (76.29%) aligned discordantly 1 time
    ----
    13990 pairs aligned 0 times concordantly or discordantly; of these:
      27980 mates make up the pairs; of these:
        21023 (75.14%) aligned 0 times
        4714 (16.85%) aligned exactly 1 time
        2243 (8.02%) aligned >1 times
  30128 (4.19%) were unpaired; of these:
    1242 (4.12%) aligned 0 times
    28237 (93.72%) aligned exactly 1 time
    649 (2.15%) aligned >1 times
98.42% overall alignment rate
```

The variant calling file looks like this:
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##samtoolsVersion=1.2+htslib-1.2.1
##samtoolsCommand=samtools mpileup -go D23_S15.bcf -f /home/ag1349/SNPTB/H37Rv/H37Rv.fa aln.sorted.bam
##reference=file:///home/ag1349/SNPTB/H37Rv/H37Rv.fa
##contig=<ID=AL123456.3,length=4411532>
##ALT=<ID=X,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">
##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
##INFO=<ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">
##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=ICB,Number=1,Type=Float,Description="Inbreeding Coefficient Binomial test (bigger is better)">
##INFO=<ID=HOB,Number=1,Type=Float,Description="Bias in the number of HOMs number (smaller is better)">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
##bcftools_callVersion=1.2+htslib-1.2.1
##bcftools_callCommand=call -vm -o D23_S15.vcf D23_S15.bcf
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	aln.sorted.bam
AL123456.3	1	.	T	G	30.3387	.	DP=2;VDB=0.02;SGB=-0.453602;MQ0F=0;AC=2;AN=2;DP4=0,0,2,0;MQ=42	GT:PL	1/1:58,6,0
AL123456.3	3	.	G	C	16.986	.	DP=4;VDB=0.02;SGB=-0.453602;RPB=1;MQB=1;BQB=1;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=1,0,2,0;MQ=36GT:PL	0/1:50,0,15
AL123456.3	4	.	A	G	26.6778	.	DP=4;VDB=0.0249187;SGB=-0.511536;MQ0F=0;AC=2;AN=2;DP4=0,0,3,0;MQ=36	GT:PL	1/1:59,11,5
AL123456.3	4	.	A	ATTGA	49.7709	.	INDEL;IDV=2;IMF=0.5;DP=4;VDB=0.02;SGB=-0.453602;MQ0F=0;AC=2;AN=2;DP4=0,0,2,0;MQ=42	GT:PL	1/1:77,6,0
```
Thus, at genome position 1, reference base T is subsituted by base G, and this SNP call is made with a quality score (QUAL) of 30.3387.

Congratulations! You have successfully identified SNPs in your sequencing data! :+1:

## 3. Data Analysis: SNP Annotation

### 3.1 Assessing data quality 

Now that you have completed data analyses, you want to see how good the read mapping was. The script get_avgdepth_genomecov.py gets the mapping rate from map.run files for each sample, calculates average depth in each sample, as well as % genome that has at least 5, at least 10, and at least 20 reads at a position.
```
$ python get_avgdepth_genomecov.py <full path to data directory>
```

The output is saved in the data directory (here, SNPTB/test) and looks like this:
```
$ more ../test/avg_depth_maprate.txt 
StrainID	Nreads	Map_rate	Npos_withreads	AvgDepth	%GenomeCov_5x	%GenomeCov_10x	%GenomeCov_20x
D23_S15	718668	98.42%	4371991	33.7377407685	99.1399113127	97.9120725546	93.0430323393
A31_S8   2382455  98.51%   4380656  113.939117566  99.6821024066  99.3841105076  98.7632902469
```

### 3.2 SNP annotation

While you can look at the SNPs in the VCF file, the format is cumbersome to look at and there is a lot of information in there that you don't need.

Thus, to extract only high-confidence SNPs from the VCF file (QUAL score > 200, _i.e._ probability that a SNP call is incorrect <1e-20), and to annotate these high-confidence SNPs with information such as the gene name and to know if the mutation is synonymous or non-synonymous, we will run the script snp_annotation.py:

```
$ python snp_annotation.py <full path to data directory>
```

This will _create_ a new folder in the SNPTB/output with the following naming convention:
<Year><Month><Date><Hour>_<data directory name>
   
Thus, for our test case that has the data in SNPTB/test, the folder "2017113013_test" is created in the "output" folder of SNPTB.

This folder contains a file called 'SNPs_qscore_200.0.txt' that looks like this (please note this file is a dummy file, and not true output of D23_S15 and A31_S8):
```
SNP,D23_S15,A31_S8,gname,mut_effect
A1977G,0,1,intergenic,NA
T4013C,1,1,recF_Rv0003,nonsyn_I245T
G7362C,1,0,gyrA_Rv0006,nonsyn_E21Q
G7585C,1,1,gyrA_Rv0006,nonsyn_S95T
G9304A,0,1,gyrA_Rv0006,nonsyn_G668D
G9516C,1,1,gyrA_Rv0006,nonsyn_E739Q
C11370T,1,1,intergenic,NA
A11879G,1,1,_Rv0008c,nonsyn_S145P
T14785C,1,1,_Rv0012,nonsyn_C233R
C18064T,1,0,pknA_Rv0015c,syn_
A24716G,1,1,fhaA_Rv0020c,syn_
```
The first line is the header that describes what each column is. Since this file is in a comma separated format, you can also open it in Excel. 

The first column shows the nucleotide mutation and the genomic position where mutation happened. The subsequent columns represent each sample in your data directory, with "1" denoting that the mutation is present in that sample and "0" indicating that the mutation is absent. Second to last column lists the gene name(s) where the mutation is (or if it is in an intergenic region), and the last column shows the functional effect of the mutation: synonymous or non-synonymous. If the mutation is non-synonymous, the WT residue, the codon number, and the mutated residue are indicated.

**Note:** The funtional effect of a mutation is calculated by translating the gene sequence from start to finish and then determining if the mutation changed any residue. _If_ the gene sequence has some untranslated region, it is not taken into consideration.

### 3.3 Compare mutations between two samples 

If you want to identify SNPs in one sample _relative_ to another sample (for example, if you sequenced your lab wild-type strain and want to find out mutations relative to this wild-type), then you need to create a file like the one in test/paired_ids.

In brief, each line should identify a pair of samples separated by a space. Script getpairedmut.py will then identify mutations in the 'SNPs_qscore_200.0.txt' file that are different in that pair of samples.

```
python getpairedmut.py <full path to your SNPs_qscore_200.0.txt file> <full path to your paired_ids.txt file>
```

This script will save another file called 'SNPs_qscore_200.0_pair_mutations.txt' in the same folder that has the 'SNPs_qscore_200.0.txt' file. This file looks like this (again, a dummy output):
```

Pair A31_S8:D23_S15	Num. SNPs different:4
A1977G	intergenic	NA
G7362C	gyrA_Rv0006	nonsyn_E21Q
G9304A	gyrA_Rv0006	nonsyn_G668D
C18064T	pknA_Rv0015c	syn_

```

## 4. Doing it all via the job scheduler

Because sequencing data files are big, the analyses usually takes a lot of time. Because HPCC is a shared resource, they have rules for fair usage. One of those rules is that you _cannot_ run scripts that take a long time (say, more than 30 minutes) from your log-in screen. If your sequencing data has a lot of samples, the qcmap.py script may run for days. 

The alternative (and the better way) is then to run via the "job scheduler". This utlity on the HPCC clusters queues the jobs (running your script = job) that are then run in the background (which means you can log-off and use your computer for other stuff than stare at it for days!).

The folder SNPTB/scheduler has the script run_SNPTB_pipeline.sh that submits your job to the scheduler. It does everything except comparing mutations between two samples (section 3.3).

Added benefit is that now you need to run _one_ script as opposed to three! 

Just like before, you need to load the software modules you'd need (if they are not loaded already, type 'module list' to see which modules are pre-loaded in your environment).
```
$ module load python/2.7.11
$ module load bowtie2/2.2.6
$ module load samtools/1.2
$ module load bcftools/1.2
$ module load java/1.8.0_66
```
Now, you need to _update_ the run_SNPTB_pipeline.sh to provide full path to your data directory, as well as change the time requirements of your job. You can open the script for editing using vi:
```
$ vi run_SNPTB_pipeline.sh
``` 
The script currently asks for 2 hours of computing time:
```
#SBATCH -t 2:00:00
```
Please be mindful of how much time you ask, if you ask way too much your future jobs may wait in the queue for longer. 

The other thing you need to change is that full path of SNPTB/code folder on your system:
```
#SBATCH -D <your_home_directory_path>/SNPTB/code
```
Now to run the script (this is a _shell_ script, note the .sh extension), type the following at the prompt:
```
$ sbatch run_SNPTB_pipeline.sh
```
You will see something like this on your screen:
```
Submitted batch job 5170110
```
This means your job was submitted successfully, and your job id is 5170110.

Now you can log-off from the HPCC, do other stuff, log back in again, and check the status of your job by typing the following:
```
$ squeue -u <netid>
```

It prints a message:
```
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           5170110      main run_SNPT   ag1349  R       0:01      1 node125
```

This shows that your job has been running (R) for 0:01 seconds. If your job is waiting in the queue, it would show the symbol "Q".

If you don't see your job listed when you type squeue command, it means that either your job was completed successfully or your job was terminated. To check that, you can inspect the SNPTB<jobid>.err and the SNPTB<jobid>.out files in the SNPTB/code folder.
   
Note that the SNPTB<jobid>.out file now contains the output of qcmap.py:
   
```
Data organization complete.
0:00:02.440814
Number of samples to be processed: 2

/home/ag1349/SNPTB/test/D23_S15
Trimming adapters and doing quality control ...
Running FastQC ...
Completed QC analysis. Files are in  <your_dir_path>/SNPTB/test/D23_S15/qc
Mapping filtered reads to the reference ...
SNP calling ...
Completed read mapping and SNP calling. Files are in  <your_dir_path>/SNPTB/test/D23_S15/mapped

/home/ag1349/SNPTB/test/A31_S8
Trimming adapters and doing quality control ...
Running FastQC ...
Completed QC analysis. Files are in  <your_dir_path>/SNPTB/test/A31_S8/qc
Mapping filtered reads to the reference ...
SNP calling ...
Completed read mapping and SNP calling. Files are in  <your_dir_path>/SNPTB/test/A31_S8/mapped
Quality control and read mapping complete.
0:42:03.958316
```


## 5. Getting ready to publish
### 5.1 What to include in the “Material and Methods” and how to cite the pipeline and its dependencies

Please see HOW_TO_CITE for this information. https://github.com/aditi9783/SNPTB/blob/master/HOW_TO_CITE.md

### 5.2 Depositing raw and meta-data to databases

Here is an excellent resource: https://github.com/faircloth-lab/home/wiki/Submitting-Read-Data-to-NCBI-SRA

