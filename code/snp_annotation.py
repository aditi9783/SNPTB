#!/usr/bin/env python

import os, subprocess
import sys
from datetime import datetime
from snpTB_functions import *

QSCORE = 200.0
# determine full path of SNPTB
SNPTBroot = sys.path[0].rstrip("/code") 

# start of getGenePos ####################
def getGenePos( genef ): # get gene start and end positions
    genes = [] # list of tuples of gene start, gene end, gene name such that start < end (can't identify complement genes)
    genefh = open(genef, 'r') # H37Rv genes
    for line in genefh: # genes are already sorted by position in this genefile
        line = line.rstrip("\n")
        content = line.split()
        start = int(content[1])
        end = int(content[2])
        gname = content[3]

        if start < end:
            genes.append([ start, end, gname ])
        else:
            genes.append([ end, start, gname ])
    genefh.close()
    return genes
# end of getGenePos #####################

# start of checkRegion ##########################
def checkRegion( pos, genes ): # check if this position is in genic region
    glist = []
    for tup in genes:
        if pos >= tup[0] and pos <= tup[1]: # pos is within gene
            glist.append(tup[2])

    if len(glist) == 0: # no genes found that contained this mutation
        return 0, glist
    else:
        return 1, glist
# end of checkRegion ##########################

# start of findMutations #################
def findMutations( seqlist ):
    mutmatrix = [] # list of high Qscore mutations in each strain
    strainids = []
    allmut = [] # all unique mutations in all strains
    #fout = open("distmat_strain_ids.out", 'w') # write the strain ids as key for the distance matrix
    counter = 0
    for sq in seqlist:
        seqdir, vcffile = sq
        #print seqdir, vcffile
        strainids.append(vcffile.rstrip(".vcf"))
    #    fout.write(str(counter)+"\t"+sd+"\n")
        mutlist = []
        fh = open(seqdir+"/"+vcffile, 'r')
        #fh = open(seqdir+"/mapped/"+vcffile, 'r')
        for line in fh:
            if line.startswith("#"):
                    continue
            else:
                content = line.split("\t")
                if float(content[5]) > QSCORE: # minimum quality score for calling a mutation
                    if len(content[3]) == 1 and len(content[4]) == 1: # both the ref base and mutated base are of len 1 => no indels, only SNPs 
                        mut = content[3]+content[1]+content[4]
                        mutlist.append( mut )
                        if mut not in allmut:
                            allmut.append( mut )
        fh.close()
        mutmatrix.append( mutlist[:] )
        # print vcffile, len(mutlist)
    #fout.close()
    return mutmatrix, strainids, allmut
# end of findMutations ###################

#start of SNPlist ##################
def SNPlist( mutmatrix, strainids, allmut, genes, geneseq, gposdict, protnames, complementgenes, outfilename ):
    
    allmutlist = []
    for i in range(len(strainids)):
        sid = strainids[i]
        sidmuts = mutmatrix[i] # all mutations in this strain
        mutbinary = []
        for j in range(len(allmut)):
            if allmut[j] in sidmuts:
                mutbinary.append("1")
            else:
                mutbinary.append("0")
        #print sid,"\n",mutbinary
        allmutlist.append(mutbinary[:])

    fout = open(outfilename, 'w')
    fout.write("SNP,"+",".join(strainids)+",gname,mut_effect\n")
    for i in range(len(allmut)):
        outstr = allmut[i]
        snp = allmut[i]
        pos = int( snp[1:-1] ) # remove 1st char (ref base) and last char (mut base)
        flag, glist = checkRegion( pos, genes )
        gname = ""
        muttype = "NA"
        if flag == 1:
            gname = ";".join(glist)
            for g in glist:
                if g in protnames:
                    muteffect = [] 
                    cflag = 0
                    if g in complementgenes:
                        cflag = 1
                    annotation = getMutEffect(g, geneseq[g], snp, gposdict, cflag)
                    muteffect.append(annotation)
            muttype = ";".join(muteffect)
        else:
            gname = "intergenic"

        for j in range(len(allmutlist)):
            outstr += ","+allmutlist[j][i]
        fout.write(outstr+","+gname+","+muttype+"\n")
    fout.close()
        
# end of SNPlist ###################

if __name__ == '__main__':
    seqlist = []
    seqids = []
    pairedsamples = []
    if len(sys.argv) == 2: # data path is provided.
        if sys.argv[1].endswith("/"):
            datapath = sys.argv[1]
        else:
            datapath = sys.argv[1]+"/"
    else:
        datapath = SNPTBroot+"/data/" # default data path
    seqids = os.walk(datapath).next()[1] # only get the subdir in the datapath, not further down subdirs

    # for output file, create a new folder in the SNPTB/output directory using same name as data directory and adding a timestamp
    pathcontents = datapath.split("/")
    datadir = pathcontents[-2]

    FORMAT = '%Y%m%d%H'
    path = datadir
    newpath = '%s_%s' % (datetime.now().strftime(FORMAT), path)
    os.makedirs(SNPTBroot+"/output/"+newpath)
    outfilename = SNPTBroot+"/output/"+newpath+"/SNPs_qscore_"+str(QSCORE)+".txt"

    for sd in seqids:
        seqdir = datapath+sd+"/mapped"
        # since problematic cases were moved to other location to allow running of parallel computing, check if the path exists
        if os.path.isdir(seqdir):
            allf = os.listdir(seqdir)
            for vfile in allf:
                if ".vcf" in vfile: # if the uncompressed .vcf file is present
                    seqlist.append( [seqdir, vfile] )
        else:
            print("Problem case:", sd, "\n")

    if len(seqlist) == 0:
        print "No VCF files were found in samples at ", datapath
        exit()

    genefile = SNPTBroot+"/H37Rv/genes_loci.txt"
    genes = getGenePos( genefile )
    geneseq, gposdict = extractSeq( SNPTBroot+"/H37Rv/H37Rv_genes.txt" ) # supply file name that has all gene seq in fasta format
    protnames, complementgenes = getH37Rv_proteins( SNPTBroot+"/H37Rv/H37Rv_genpept.gp" ) # not all genes code for proteins. Get gene names of protein coding genes only.

    mutmatrix, strainids, allmut = findMutations( seqlist )
    SNPlist( mutmatrix, strainids, allmut, genes, geneseq, gposdict, protnames, complementgenes, outfilename )
