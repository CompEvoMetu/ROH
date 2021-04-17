#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Dilek Koptekin
"""

from pybedtools import BedTool
import pandas as pd
import numpy as np
import sys

snps=sys.argv[1]
file=sys.argv[2] + ".new.bed"
fasta1 = open((sys.argv[2] + ".1.fa"), 'w')
fasta2 = open((sys.argv[2] + ".2.fa"), 'w')

fasta="/mnt/NEOGENE2/share/ref/genomes/hsa/hs37d5.fa"


# get all sequences
df = pd.read_csv(file, sep=' ', comment='t', header=None)
df.columns = ['chr', 'posStart', 'posEnd','']
df = df[df['posStart'] != df['posEnd']]
allseq = pd.read_csv(BedTool.from_dataframe(df).sequence(fi=fasta, tab=True, name=True).seqfn, sep="\t", names=["snpID","seq"])

# read snp file
snp = pd.read_csv(snps, sep='\t', comment='t', header=None)
snp.columns = ['chr', 'posStart', 'posEnd','allele1', 'allele2']
snp['chr'] = 1
snp['snpID'] = "1_" + snp['posStart'].astype(str) + "_" + snp['posEnd'].astype(str)

# change given position on seqs
#fa1
allseq.loc[allseq.snpID.isin(snp.snpID), ['seq']] = snp[['allele1']].values
fa1 = "".join(allseq['seq'])
fasta1.writelines(">1\n")
fa1_parts = [fa1[i:i+60] for i in range(0, len(fa1), 60)]
fasta1.writelines("\n".join(fa1_parts))

#fa2
allseq.loc[allseq.snpID.isin(snp.snpID), ['seq']] = snp[['allele2']].values
fa2 = "".join(allseq['seq'])
fasta2.writelines(">1\n")
fa2_parts = [fa2[i:i+60] for i in range(0, len(fa2), 60)]
fasta2.writelines("\n".join(fa2_parts))

fasta1.close()
fasta2.close()
