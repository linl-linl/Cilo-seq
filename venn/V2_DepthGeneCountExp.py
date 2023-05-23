#!/usr/bin/env python
#coding=utf-8

import os
import re
import numpy as np
import codecs
import argparse

parser = argparse.ArgumentParser(description='取同一测序深度，统计基因表达情况，以count表示： \
                                              python2 VDepthGeneCountExp.py \
                                              -s list_sam \
                                              -g url_to_gtf \
                                              -d depth')

parser.add_argument('-s', '--sam', type=str, default='list_sam', help='a list of sam file')
parser.add_argument('-g', '--gtf', type=str, default='/home/songjia/reference/hg19/gencode.v35lift37.annotation.gtf', help='gtf for htseq-count')
parser.add_argument('-d', '--depth', type=str, required=True, help='depth')

args = parser.parse_args()

def DepthGene(sam, fq, d):
    file_fq = open(fq, 'r')
    cont_fq = file_fq.readlines()
    nn = 0
    nn_max = len(cont_fq)
    pattern = r'@(.*?)\s'
    list_id = []
    while (nn<nn_max):
        if nn%4==0:
            idid = re.match(pattern, cont_fq[nn])
            list_id.append(idid.group(1))
        nn = nn+1
    list_ID = np.random.choice(list_id, d, replace = False)
    dict0 = {}
    for i in list_id:
        dict0[i] = 0
    for i in list_ID:
        dict0[i] = 1
    file_fq.close()

    file_align = open(sam, 'r')
    prefix = re.sub('.sam','',sam)
    file_align_id = codecs.open(prefix+'_'+str(d)+'DepthGene.out.sam','w')
    for line in file_align:
        if line.startswith('@'):
            file_align_id.write(line)
        else:
            ll = line.split('\t')
            if dict0[ll[0]]:
                file_align_id.write(line)
    file_align.close()
    file_align_id.close()

list_sam=[]
file=open(args.sam,'r')
for i in file:
    sam_name=re.sub('\n','',i)
    list_sam.append(sam_name)

d = int(args.depth)
gtf = args.gtf
for i in list_sam:
    sam=i
    fq=re.sub('.sam','.fq',i)
    url_fq=fq
    DepthGene(sam, url_fq, d)
    prefix = re.sub('.sam','',sam)
    DepthGene_out=prefix+'_'+str(d)+'DepthGene.out.sam'

    htseq_out=prefix+'_'+str(d)+'_counts.txt'
    os.system("htseq-count -f sam -r name -s no -a 10 -t gene -i gene_id -m intersection-nonempty %s %s > %s" %(DepthGene_out, gtf, htseq_out))
    
    count_adj = prefix+'_'+str(d)+'_count.txt'
    file_counts = open(htseq_out, 'r')
    file_count = open(count_adj, 'w')
    cont_counts = ''
    cont_count = ''
    cont_counts = file_counts.readlines()
    cont_count = cont_counts[:-5]
    for w in cont_count:
        file_count.write(w)
    file_counts.close()
    file_count.close()
    os.system("rm %s" % (htseq_out))
    print(i+' is finished!')

