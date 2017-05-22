#!/usr/bin/env python

# Load WExT
import sys
# sys.path.append('wext/')
from wext import re_test, SADDLEPOINT

# Hard-code parameters
N = 779 # number of samples
X = [] # number of mutations per gene
inp = open ("Results_Figures_and_P_Values/genes_top_16_genes_with_VUS.txt","r")
for line in inp.readlines():
	for i in line.split(','):
		X.append(int(i))

T =  168 # 11 overlaps; total number of exclusive mutations
tbl = []
inp = open ("Results_Figures_and_P_Values/vec_top_16_genes_with_VUS.txt","r")
for line in inp.readlines():
	for i in line.split(','):
		tbl.append(int(i))
# Run WExT saddlepoint
print re_test(T, X, tbl, method=SADDLEPOINT)
