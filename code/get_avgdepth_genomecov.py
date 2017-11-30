#!/usr/bin/python

import os, sys

# determine full path of SNPTB
SNPTBroot = sys.path[0].rstrip("/code")

if len(sys.argv) == 2: # data path is provided.
    if sys.argv[1].endswith("/"):
        datapath = sys.argv[1]
    else:
        datapath = sys.argv[1]+"/"
else:
    datapath = SNPTBroot+"/data/" # default data path

seqids = os.walk(datapath).next()[1] # only get the subdir in the seqpath, not further down subdirs
outfh = open(datapath+"avg_depth_maprate.txt", 'w')
outfh.write("StrainID\tNreads\tMap_rate\tNpos_withreads\tAvgDepth\t%GenomeCov_5x\t%GenomeCov_10x\t%GenomeCov_20x\n")

for sd in seqids:
    mpileupfh = open(datapath+sd+"/mapped/aln.sorted.bam.mpileup", 'r')
    npos_gt5, npos_gt10, npos_gt20 = 0,0,0 # num genome positions with depth >5, >10, and >20
    depth = [] # reads at each position
    npos = 0 # num of genome positions encountered

    for line in mpileupfh:
        content = line.split()
        cov = int(content[3])
        depth.append(cov)
        npos += 1
        if cov >= 5:
            npos_gt5 += 1
        if cov >= 10:
            npos_gt10 += 1
        if cov >= 20:
            npos_gt20 += 1

    avgcov = sum(depth)/float(npos)
    readcov5 = float(npos_gt5)/npos * 100.0
    readcov10 = float(npos_gt10)/npos * 100.0
    readcov20 = float(npos_gt20)/npos * 100.0
    # get mapping rate
    mapfh = open(datapath+sd+"/mapped/map.run", 'r')
    firstline = mapfh.readline()
    contents_1 = firstline.split()
    Nreads = contents_1[0]
    lastline = mapfh.readlines()[-1]
    contents_last = lastline.split()
    maprate = contents_last[0]
    outfh.write(sd+"\t"+Nreads+"\t"+maprate+"\t"+str(npos)+"\t"+str(avgcov)+"\t"+str(readcov5)+"\t"+str(readcov10)+"\t"+str(readcov20)+"\n")
outfh.close()
