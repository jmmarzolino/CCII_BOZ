#!/usr/bin/env python
import sys
import re
import os
import pysam
import pysamstats
#
def getpos(inchr,inpos,inA1,inA2):
    i=0
    for rec in pysamstats.stat_variation_strand(samfile, chrom=inchr, start=inpos-1, end=inpos, truncate=True, fafile="/rhome/jmarz001/shared/GENOMES/BARLEY/2019_Release_Morex_vers2/Barley_Morex_V2_pseudomolecules.fasta"):
        i+=1
        #return rec['chrom'], rec['pos']+1, rec["A_fwd"], rec["T_fwd"], rec["G_fwd"], rec["C_fwd"], rec["A_rev"], rec["T_rev"], rec["G_rev"], rec["C_rev"]
        if rec[inA1]>0 and rec[inA2]>0:
            return rec['chrom'], rec['pos']+1, rec[inA1], rec[inA2],"0/1"
        elif rec[inA1]>0 and rec[inA2]==0:
            return rec['chrom'], rec['pos']+1, rec[inA1], rec[inA2],"0/0"
        elif rec[inA1]==0 and rec[inA2]>0:
            return rec['chrom'], rec['pos']+1, rec[inA1], rec[inA2],"1/1"
        elif rec[inA1]==0 and rec[inA2]==0:
            return rec['chrom'], rec['pos']+1, rec[inA1], rec[inA2],"./."
        return rec['chrom'], rec['pos']+1, rec[inA1], rec[inA2]
    if i==0: return inchr,inpos,0,0,"./."
samfile = pysam.AlignmentFile(sys.argv[1], "rb")
f = open(sys.argv[2])
while True:
    line = f.readline()
    line = line.rstrip('\n')
    if not line: break
    lsplit = line.split()
    outdata = getpos(lsplit[0],int(lsplit[1]),lsplit[2],lsplit[3])[0:4]
    forprint = "\t".join(map(str,outdata))
    print forprint
f.close()
