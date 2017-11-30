#!/usr/bin/python

import os, sys
from datetime import datetime
startTime = datetime.now()

if len(sys.argv) == 2: # data path is provided.
    datapath = sys.argv[1]+"/"
else:
    datapath = "../data/" # default data path

allfiles = os.listdir(datapath)
if len(allfiles) == 0:
    print "No files in location ", datapath
    exit()

for f in allfiles:
    if "R1_001" in f:
        dirname, suffix = f.split('_R1')
        fullpath = datapath+dirname
        os.makedirs(fullpath)
        os.system("cp "+fullpath+"_R*.gz "+fullpath)
print "Data organization complete."
print datetime.now() - startTime
