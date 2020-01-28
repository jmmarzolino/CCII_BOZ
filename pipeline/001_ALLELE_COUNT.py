#!/usr/bin/env python
import sys
import re
import os
#

f = open(sys.argv[1])
while True:
    line = f.readline()
    line =line.rstrip('\n')
    if not line: break
    if re.match("^\#", line) is not None:
        next
    else:
        split2 = line.split()
        alleleR = 0
        alleleA = 0
        for index, item in enumerate(split2[9:]):
            split3 = item.split(":")
            if split3[0] == './.':
                next
            elif split3[0] == '0/0' or split3[0] == '0|0':
                alleleR += 2
            elif split3[0] == '0/1' or split3[0] == '0|1' or split3[0] == '1/0' or split3[0] == '1|0':
                alleleR += 1
                alleleA += 1
            elif split3[0] == '1/1' or split3[0] == '1|1':
                alleleA += 2
        outsite = [split2[0],split2[1],split2[3],split2[4],alleleR,alleleA]
        outputf = '\t'.join(map(str,outsite))
        print outputf
f.close()
