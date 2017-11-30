#!/usr/bin/python

import math
from genetic_code import code

# define base complements
complement = {"A" : "T", "T" : "A", "C" : "G", "G" : "C"}

##########################################
def readTruePairs( fname ): # read the file name that has the true pairs ids and the mutational distance between them
    truepairfh = open(fname, 'r')
    pairids = [] # list of pairs whose SNPs are to be compared to each other
    for line in truepairfh:
        line = line.rstrip("\n")
        pairids.append( line.split() )
    truepairfh.close()
    return pairids
##########################################

# start of extractSeq ##########################
def extractSeq( fname ): # extract gene or protein seq from respective file handles
    fh = open(fname, 'r')
    seqdict = {}
    gposdict = {}
    seq = ""
    for line in fh:
        if ">" in line: # header line
            if seq != "": # there is some sequence to save
                seqdict[gname] = seq.replace('\n','')
            content = line.split()
            info = {'gene' : '', 'locus_tag' : '', 'location' : ''}
            for i in range(1, len(content)):
                key, val = content[i].split('=')
                key = key.lstrip('[')
                val = val.rstrip(']')
                if key in info: # H37Rv_proteins has a lot more fields, only get these three
                    info[key] = val
                else:
                    continue
            gname = info['gene']+"_"+info['locus_tag']
            locstr = info['location']
            if locstr.startswith("complement"): # info[2] is location string [location=12468..13016] or [location=complement(13133..13558)]
                locstr = locstr.lstrip("complement(")
                locstr = locstr.rstrip(")]")
                #gname = gname+"_c" # not all complement genes end with a 'c'. Thus add this distinguishing feature to the gene name
            elif locstr.startswith("order"): # only one gene has this. [locus_tag=Rv3216] [location=order(3593369..3593437,3593439..3593852)]
                # enter this gene manually, encompassing the entire region covered by this gene
                locstr = "3593369..3593852"
            start, end = locstr.split('..')
            start = start.lstrip('<') # some gene start positions have location mentioned as [location=<2550340..2551326]. Remove the "<" sign.
            end = end.lstrip('>') # some gene end positions have location mentioned as 817531..>817866. Remove the '>' sign.
            gposdict[gname] = [int(start), int(end)]
            #print(gname, start, end)
            seq = ""
        else:
            seq += line
    return seqdict, gposdict
# end of extractSeq ##########################

# start of getH37Rv_proteins #################
def getH37Rv_proteins(fname):
    fh = open(fname, 'r')
    protnames = [] # names of protein coding genes only (no tRNA genes)
    compgenes = [] # complement genes
    flag = 1
    gname = ''
    locustag = ''
    for line in fh:
        if line.startswith("LOCUS"):
            flag = 1
        if flag == 1:
            if "/gene=" in line:
                content = line.split('="')
                gname = content[1].rstrip('"\n')
            elif "/locus_tag" in line:
                content = line.split('="')
                locustag = content[1].rstrip('"\n')
            elif "/coded_by" in line:
                if "complement" in line: # complement gene
                    compgenes.append(gname+"_"+locustag)
                protnames.append(gname+"_"+locustag)
                flag = 0
                gname = ''
                locustag = ''
    return protnames, compgenes
# end of getH37Rv_proteins #################

# start of translateGene ##########################
def translateGene( nt ):
    tseq = ["M"] # no matter the first codon, the amino acid seq always starts with M -> NOT TRUE! Few proteins in H37Rv don't start with M
    for i in range(3, len(nt), 3): # start from 2nd codon, and then read gene sequence in steps of 3
        codon = "".join(nt[i:i+3])
        if codon in code:
            tseq.append( code[codon] )
        else: # some indels change frame, so last codon might be incomplete
            tseq.append( "***" )
    trans_aa = "".join(tseq)
    #trans_aa = trans_aa.rstrip('*') # remove any stop codons at the end
    return trans_aa
# end of translateGene ##########################

# start of checkRegion ##########################
def checkRegion( pos, repeats ): # check if this position is in the repeat regions
    rregion = []
    for tup in repeats:
        if pos >= tup[0] and pos <= tup[1]: # pos is within the repeat region 
            rregion = tup
            break
    return rregion
# end of checkRegion ##########################

##########################################
def getMutEffect( gname, gseq, snp, gposdict, compflag ):
    ntseq = list(gseq)
    pos = int(snp[1:-1])
    refbase = snp[:1]
    mutbase = snp[-1:]
    wtprot = translateGene(ntseq)
    gstart, gend = gposdict[gname]
    mutsite = pos - gstart
    mutseq = ntseq[:]
    if compflag == 1: # complement genes
        mutsite = -mutsite-1
        mutbase = complement[mutbase]
        refbase = complement[refbase]
        region = "".join(ntseq[mutsite-2:mutsite+2+1])
        #print (snp, gname, region)
    if ntseq[mutsite] == refbase:
        mutseq[mutsite] = mutbase
    else:
        print("ERROR! Refbase doesn't match!", snp, gname, "refbase_snp:", refbase, "refbase_seq:", ntseq[mutsite], "\n")
        print(ntseq[mutsite-2:mutsite+2+1])
    mutprot = translateGene(mutseq)
    wtseq = list(wtprot)
    mutatedseq = list(mutprot)
    aamut = ""
    for i in range(len(wtseq)):
        if wtseq[i] != mutatedseq[i]: 
            #print (i, wtseq[i], mutatedseq[i])
            aamut = wtseq[i]+str(i+1)+mutatedseq[i]
    if wtprot != mutprot:
        #print( gname, snp, "Non syn mutation:\n", wtprot, "\n", mutprot)
        return "nonsyn_"+aamut
    else:
        return "syn_"+aamut
##########################################

##########################################
def getMutData( infile, pairids ): # read csv file that has mutation pr/ab 1/0 for each SNP for each strain
    outfile = infile.rstrip(".txt")+"_pair_mutations.txt"
    fh = open(infile, 'r')
    fout = open(outfile, 'w')

    header = fh.readline().split(',')
    # get column index matching the index-case and secondary-case ids
    idcolumn = {} # key: sample id (index or secondary), value: column id
    for i in range(1,len(header)-2): # 1st elt is "SNP" and last two cols are "gname" and "mut_effect", so ignore those 
        sid = header[i]
        idcolumn[sid] = i

    # save mutations in each pair
    pairmutsnp = {} # key: pair ids, value: list of mutations that are different in this pair 
    pairmutgene = {} # key: pair ids, value: list of genes that have mutations in this pair 
    snp_annotation = {} # key: snp, value = mutation type (syn, nonsyn, intergenic, etc)
    for pids in pairids: # initialize the above dicts
        idstr = ":".join(pids)
        pairmutsnp[idstr] = []
        pairmutgene[idstr] = []

    for line in fh:
        line = line.rstrip('\n')
        content = line.split(',')
        # get pairs where this snp differs in the index and sec cases
        for pids in pairids: # for each index case
            idstr = ":".join(pids)
            idxcol = idcolumn[pids[0]] # 1st sample in pair
            sidcol = idcolumn[pids[1]] # 2nd sample in pair
            if content[idxcol] == content[sidcol]:
                continue # mutation either present or absent in both samples of a pair
            else:
                snp = content[0]
                gname = content[-2]
                muttype = content[-1] 

                pairmutsnp[idstr].append(snp)
                pairmutgene[idstr].append(gname)
                snp_annotation[snp] = muttype

    # print mutations in each pair
    for pids in pairids:
        idstr = ":".join(pids)
        fout.write("\nPair "+idstr+"\tNum. SNPs different:"+str(len(pairmutsnp[idstr]))+"\n") #, set(pairmutgene[idstr]))
        for i in range(len(pairmutsnp[idstr])):
            fout.write(pairmutsnp[idstr][i]+"\t"+pairmutgene[idstr][i]+"\t"+snp_annotation[pairmutsnp[idstr][i]]+"\n")
    fout.close()
##########################################

