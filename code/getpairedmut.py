#!/usr/bin/python

from snpTB_functions import *
import sys

## MAIN ##
if len(sys.argv) != 3:
    print "Please provide full paths for the SNPs_score file and the paired ids file, in that order."
    exit()
elif len(sys.argv) == 3:
    truepairf = sys.argv[2] 
    pairids = readTruePairs(truepairf) # 

    # read true mutation file and identify the true mutations different in each pair
    infile = sys.argv[1] 
    getMutData( infile, pairids )

