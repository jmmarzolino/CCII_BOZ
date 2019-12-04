#!/usr/bin/env python
import sys
import re
import os
import scipy.stats as stats
######################################################
def sumac(a1,a2,pc):
    if a1>=a2:
        pc[0]+=int(a1)
        pc[1]+=int(a2)
    else:
        pc[0]+=int(a2)
        pc[1]+=int(a1)
def fst(a1,n1,a2,n2):
    h1 = (a1*(n1-a1))/(n1*(n1-1))
    h2 = (a2*(n2-a2))/(n2*(n2-1))
    N =(a1/n1 - a2/n2)**2 - h1/n1 - h2/n2
    D = N + h1 + h2
    return([N,D])
######################################################
n18 = [0,0]
n28 = [0,0]
n58 = [0,0]
nP = [0,0]
NDP18 = [0,0]
NDP28 = [0,0]
NDP58 = [0,0]
seg18=0
seg28=0
seg58=0
segP=0
#######################################################
f = open(sys.argv[1])
while True:
    line = f.readline()
    if not line: break
    split = line.split('\t')
    sumac(split[2],split[3],nP)
    sumac(split[4],split[5],n18)
    sumac(split[6],split[7],n28)
    sumac(split[8],split[9],n58)
    if int(split[2])>0 and int(split[3])>0:
        segP+=1
    if int(split[4])>0 and int(split[5])>0:
        seg18+=1
    if int(split[6])>0 and int(split[7])>0:
        seg28+=1
    if int(split[8])>0 and int(split[9])>0:
        seg58+=1
#######################################################
    Nk,Dk = fst(float(split[3]),
    float(split[2])+float(split[3]),
    float(split[5]),
    float(split[4])+float(split[5]))
    NDP18[0]+=Nk
    NDP18[1]+=Dk
########################################################
    Nk,Dk = fst(float(split[3]),
    float(split[2])+float(split[3]),
    float(split[7]),
    float(split[6])+float(split[7]))
    NDP28[0]+=Nk
    NDP28[1]+=Dk
#######################################################
    Nk,Dk = fst(float(split[3]),
    float(split[2])+float(split[3]),
    float(split[9]),
    float(split[8])+float(split[9]))
    NDP58[0]+=Nk
    NDP58[1]+=Dk
f.close
win=sys.argv[2].split(":")
pos=win[1].split("-")
#################################################################
if NDP18[1]==0:
    NDP18out = 'NA'
else:
    NDP18out = NDP18[0]/NDP18[1]
if NDP28[1]==0:
    NDP28out = 'NA'
else:
    NDP28out = NDP28[0]/NDP28[1]
if NDP58[1]==0:
    NDP58out = 'NA'
else:
    NDP58out = NDP58[0]/NDP58[1]
#################################################################
if float(nP[0]+nP[1])==0:
    nPout = 'NA'
else:
    nPout = float(2*nP[0]*nP[1])/float(nP[0]+nP[1])**2
if float(n18[0]+n18[1])==0:
    n18out = 'NA'
else:
    n18out = float(2*n18[0]*n18[1])/float(n18[0]+n18[1])**2
if float(n28[0]+n28[1])==0:
    n28out = 'NA'
else:
    n28out = float(2*n28[0]*n28[1])/float(n28[0]+n28[1])**2
if float(n58[0]+n58[1])==0:
    n58out = 'NA'
else:
    n58out = float(2*n58[0]*n58[1])/float(n58[0]+n58[1])**2
print(
win[0]+"\t"+
pos[0]+"\t"+
pos[1]+"\t"+
str(nPout)+"\t"+
str(n18out)+"\t"+
str(n28out)+"\t"+
str(n58out)+"\t"+
str(segP)+"\t"+
str(seg18)+"\t"+
str(seg28)+"\t"+
str(seg58)+"\t"+
str(NDP18out)+"\t"+
str(NDP28out)+"\t"+
str(NDP58out)
)
