# -*- coding: utf-8 -*-
"""
Created on Mon May 16 11:22:04 2022

@author: jingy
"""
import statistics

aList = list(input("Amino acid sequence: "))

hydrovalues = {"I":4.5,
"V":4.2,
"L":3.8,
"F":2.8,
"C":2.5,
"M":1.9,
"A":1.8,
"G":-0.4,
"T":-0.7,
"S":-0.8,
"W":-0.9,
"Y":-1.3,
"P":-1.6,
"H":-3.2,
"E":-3.5,
"Q":-3.5,
"D":-3.5,
"N":-3.5,
"K":-3.9,
"R":-4.5,
}
hList = [hydrovalues.get(item) for item in aList]
hAvg = statistics.mean(hList)

pivalues = {"A":6.11,
"R":10.76,
"N":5.43,
"D":2.98,
"C":5.15,
"E":3.08,
"Q":5.65,
"G":6.06,
"H":7.64,
"I":6.04,
"L":6.04,
"K":9.47,
"M":5.71,
"F":5.76,
"P":6.30,
"S":5.70,
"T":5.60,
"W":5.88,
"Y":5.63,
"V":6.02,
}
pList = [pivalues.get(item) for item in aList]
piAvg = statistics.mean(pList)

mvalues = {"A":89,
"R":174,
"N":132,
"D":133,
"C":121,
"Q":146,
"E":147,
"G":75,
"H":155,
"I":131,
"L":131,
"K":146,
"M":149,
"F":165,
"P":115,
"S":105,
"T":119,
"W":204,
"Y":181,
"V":117,
}
mList = [mvalues.get(item) for item in aList]
mAvg = statistics.mean(mList)

print("Hydropathy: " + str(hAvg))
print("Isoelectric point: " + str(piAvg))
print("Molecular weight (Da): " + str(mAvg))

import re
import os
import pandas as pd
os.chdir(r'C:\Users\jingy\Downloads')

df = pd.read_csv('sr exp distributions.csv')

conserved_yyc = "YYC"
conserved_wgq = "WGQ"
string = input("Amino acid sequence: ")
match_yyc = re.search(f'({conserved_yyc})+',string)
match_wgq = re.search(f'({conserved_wgq})+',string)
if match_yyc == None:
    print('HCDR3 Length na')
else:
    start_yyc, end_yyc = match_yyc.span()
    start_wgq, end_wgq = match_wgq.span()
    cdr3start, cdr3end = end_yyc, start_wgq
    print("Heavy chain CDR3 Length: " + str(start_wgq-end_yyc))
    cdrdist = df['cdr dist'].tolist()
    cdrPerc = ((sum(i<=(start_wgq-end_yyc) for i in cdrdist))/648)*100
    print("VH CDR3 Length Percentile: " + str(cdrPerc) + "%")

pidist = df['pi dist'].tolist()
piPerc = ((sum(i<=piAvg for i in pidist)) / 648)*100
print("pI Percentile: " + str(piPerc) + "%")

hydrodist = df['hydro dist'].tolist()
hydroPerc = ((sum(i<=hAvg for i in hydrodist)) / 648)*100
print("Hydropathy Percentile: " + str(hydroPerc) + "%")

mwdist = df['mw dist'].tolist()
mwPerc = ((sum(i<=mAvg for i in mwdist))/648)*100
print("Molecular Weight Percentile: " + str(mwPerc) + "%")

from matplotlib import pyplot as plt

cdrhist = plt.hist(cdrdist)
cdrhist_x = plt.xlabel('HCDR3 lengths')
plt.axvline(start_wgq-end_yyc, ymin=0,ymax=180,color='deeppink',label = 'User Input Sequence')
plt.legend(loc='upper right')
plt.show()

pihist = plt.hist(pidist)
pihist_x = plt.xlabel('pI values')
plt.axvline(piAvg, ymin=0,ymax=180,color='deeppink',label = 'User Input Sequence')
plt.legend(loc='upper right')
plt.show()

hydrohist = plt.hist(hydrodist)
hydrohist_x = plt.xlabel('Hydropathy values')
plt.axvline(hAvg, ymin=0,ymax=180,color='deeppink',label = 'User Input Sequence')
plt.legend(loc='upper right')
plt.show()

mwhist = plt.hist(mwdist)
mwhist_x = plt.xlabel('Molecular weight (Da)')
plt.axvline(mAvg, ymin=0,ymax=180,color='deeppink',label = 'User Input Sequence')
plt.legend(loc='upper right')
plt.show()


