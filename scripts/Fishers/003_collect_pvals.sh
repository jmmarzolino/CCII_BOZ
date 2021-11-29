#!/bin/bash -l

cd /bigdata/koeniglab/jmarz001/CCII_BOZ/results
LENGTH=$(ls -l filtered_Fishers_pvalues_allpops* | wc -l)

rm pvals
cat filtered_Fishers_pvalues_allpops{1..$LENGTH} >> pvals
