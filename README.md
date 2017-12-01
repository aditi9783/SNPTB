# SNPTB
_A Bioinformatics Pipeline for Analyses of Mycobacterium tuberculosis Whole Genome Sequencing Data_


**Aditi Gupta**

Center for Emerging Pathogens, Rutgers University, Newark NJ 07103

Email: ag1349 [at] njms.rutgers.edu

October 18, 2017

Please see [docs](docs) for tutorial on how to use SNPTB.

Citation information is in [HOW_TO_CITE](HOW_TO_CITE.md).

SNPTB is publicly available under the GNU GPL License. 

SNPTB contains the following directories:
* **code**: contains all scripts
* **software**: contains Trimmomatic and FastQC software
* **scheduler**: contains script to run the SNPTB pipeline on SLURM job scheduler
* **data**: empty directory, option to save raw sequencing data here
* **output**: empty directory. When SNPTB is run, a new folder is created in "output" that contains a file with SNP annotation
* **test**: contains sample file that shows the format to specify pairs of samples for pairwise SNP identification
* **H37Rv**: all files of the reference _Mycobacterium tuberculosis_ genome H37RV needed for running SNPTB

