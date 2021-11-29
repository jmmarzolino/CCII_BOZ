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
#def fst(a1,n1,a2,n2):
#    h1 = (a1*(n1-a1))/(n1*(n1-1))
#    h2 = (a2*(n2-a2))/(n2*(n2-1))
#    N =(a1/n1 - a2/n2)**2 - h1/n1 - h2/n2
#    D = N + h1 + h2
#    return([N,D])
######################################################
n = [0,0]
allsite=0
seg=0
#######################################################
f = open(sys.argv[1])
while True:
    line = f.readline()
    if not line: break
    split = line.split('\t')
    sumac(split[4],split[5],n)
    allsite+=1
    if int(split[4])>0 and int(split[5])>0:
        seg+=1
f.close
win=sys.argv[2].split(":")
pos=win[1].split("-")
#################################################################
if float(n[0]+n[1])==0:
    Heout = 'NA'
else:
    Heout = float(2*n[0]*n[1])/float(n[0]+n[1])**2

print(
win[0]+"\t"+
pos[0]+"\t"+
pos[1]+"\t"+
str(allsite)+"\t"+
str(seg)+"\t"+
str(Heout)
)
