#!/usr/bin/env python

import os, subprocess, sys
from datetime import datetime
startTime = datetime.now()

# determine full path of SNPTB
SNPTBroot = sys.path[0].rstrip("/code") 

def qcmap( seqinfo ): # Input: tuple of [seqdir, sd, fq_1, fq_2]
    trim_path = SNPTBroot+"/software/Trimmomatic-0.36/trimmomatic-0.36.jar"
    adaptor_path = SNPTBroot+"/software/Trimmomatic-0.36/adapters/TruSeq3-PE.fa"
    seqdir, sd, fq_1, fq_2 = seqinfo
    # create a new dir for trimmomatic output (qc)
    qcdir = seqdir+"/qc"
    os.makedirs(qcdir)
    os.chdir(qcdir)
    print "Trimming adapters and doing quality control ..."
    trim_cmd = "java -jar "+trim_path+" PE -phred33 "+fq_1+" "+fq_2+" s1_pe.fq.gz s1_se.fq.gz s2_pe.fq.gz s2_se.fq.gz ILLUMINACLIP:"+adaptor_path+":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:20 &> qc.run"
    os.system(trim_cmd)
    # run FastQC to detect qual control
    qcfiles = ["s1_pe.fq.gz", "s1_se.fq.gz", "s2_pe.fq.gz", "s2_se.fq.gz"]
    print "Running FastQC ..."
    for f in qcfiles:
        fastqc_cmd = SNPTBroot+"/software/FastQC/fastqc "+qcdir+"/"+f+" &>fastqc.run" # save STDERR output 
        #print fastqc_cmd
        os.system(fastqc_cmd)
    print "Completed QC analysis. Files are in ",qcdir

    # map the reads to the reference H37Rv using bowtie2 (reference has already been indexed)
    # assume that bowtie2 has been loaded in the environment
    mapdir = seqdir+"/mapped"
    os.makedirs(mapdir)
    os.chdir(mapdir)
    print "Mapping filtered reads to the reference ..."
    bowtie_cmd = "bowtie2 -p 1 -x "+SNPTBroot+"/H37Rv/H37Rv -1 "+qcdir+"/s1_pe.fq.gz -2 "+qcdir+"/s2_pe.fq.gz -U "+qcdir+"/s1_se.fq.gz,"+qcdir+"/s2_se.fq.gz -S "+mapdir+"/aln.sam &> map.run"
    #print bowtie_cmd
    os.system(bowtie_cmd)
            
    # use samtools (assume pre-loaded) to convert sam to bam file, and then to sort the bam file
    os.system("samtools view -bS aln.sam > aln.bam")
    os.system("samtools sort aln.bam aln.sorted")
    os.system("samtools mpileup -f "+SNPTBroot+"/H37Rv/H37Rv.fa aln.sorted.bam > aln.sorted.bam.mpileup")

    print "SNP calling ..."
    # the following commands work for samtools/1.2 and bcftools/1.2, both of which should be loaded in env
    os.system("samtools mpileup -go "+sd+".bcf -f "+SNPTBroot+"/H37Rv/H37Rv.fa aln.sorted.bam")
    os.system("bcftools call -vm -o "+sd+".vcf "+sd+".bcf")
    os.system("rm "+sd+".bcf") # delete bcf file as mpileup output is already stored

    # change sorted bam file name to include strain name (subdir name) as well
    os.system("mv aln.sorted.bam "+sd+".aln.sorted.bam")

    # create index for reading bam file in Tablet
    os.system("samtools index "+sd+".aln.sorted.bam")

    # delete sam file, unsorted bam file, raw fastq files, and quality controlled fastq files
    os.system("rm aln.sam")
    os.system("rm aln.bam")
    os.chdir(qcdir)
    os.system("rm *_*.fq.gz")
    os.system("rm "+fq_1) # fq_1 has full path to the file
    os.system("rm "+fq_2) # fq_2 has full path to the file
    print "Completed read mapping and SNP calling. Files are in ", mapdir

if __name__ == '__main__':
    seqlist = []
    if len(sys.argv) == 2: # data path is provided.
        if sys.argv[1].endswith("/"):
            datapath = sys.argv[1]
        else:
            datapath = sys.argv[1]+"/"
    else:
        datapath = SNPTBroot+"/data/" # default data path
    seqids = os.walk(datapath).next()[1] # only get the subdir in the datapath, not further down subdirs
    if len(seqids) == 0: # no sample directories to process
        print "No sample directories here: ", datapath
        exit()
    for sd in seqids:
        #if len(seqlist) == 2: # for testing
        #    break   
        seqdir = datapath+sd
        #print("\n",seqdir)
        allf = os.listdir(seqdir)
        fq_1 = ""
        fq_2 = ""
        for seqf in allf: # get paired end files
            if "R1_001.fastq.gz" in seqf:
                fq_1 = seqdir+"/"+seqf
            if "R2_001.fastq.gz" in seqf:
                fq_2 = seqdir+"/"+seqf

        seqlist.append([seqdir, sd, fq_1, fq_2])
    
    print "Number of samples to be processed:", len(seqlist)
    for sq in seqlist:
        print "\n", sq[0] 
        qcmap(sq)

print "Quality control and read mapping complete."
print datetime.now() - startTime
