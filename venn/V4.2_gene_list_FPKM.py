# !/usr/bin/env python
# -*- coding:utf-8 -*-

# python Vgene_list_FPKM.py fpkm.txt
# 从[sample_name]_[depth]_fpkm.txt文件中提取gene list

import sys
import re
from decimal import Decimal

file1 = open(sys.argv[1], 'r')
list1 = []
for line in file1:
    line2 = re.sub('\n','',line)
    ll = line2.split('\t')
    fpkm = Decimal(ll[1])
    if (fpkm > 1):
        gene = re.sub('\"', '', ll[0])
        list1.append(gene)
file1.close()

prefix1 = re.sub('_FPKM.txt', '', sys.argv[1])
file2 = open(prefix1+'_genelist_FPKM.txt', 'w')
file2.write('GENE\n')
for i in list1:
    file2.write(i+'\n')
file2.close()

