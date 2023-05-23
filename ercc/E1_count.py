# ！usr/bin/env python
# -*- coding:utf-8 -*-

import os
import re
import argparse

parser = argparse.ArgumentParser(description='ERCC第一步，gene-readcount统计; \
                                              ERCC第二步，删除counts.txt文件末尾的统计信息，生成count.txt: \
                                              python ERCCcount.py \
                                              -g gtf \
                                              -s list_ERCCcount_sam.txt')

parser.add_argument('-g', '--gtf', type=str, default='/home/songjia/reference/ERCC/ERCC_Controls.gtf', help='ERCC_Controls.gtf')
parser.add_argument('-s', '--sam', type=str, default='', help='list_ERCCcount_sam.txt')

args = parser.parse_args()

if args.sam:
    list_sam = []
    file = open(args.sam, 'r')
    for i in file:
        sam_name = re.sub('\n', '', i)
        list_sam.append(sam_name)
else:
    os.system("ls *Aligned.out.sam > list_ERCCcount_sam.txt")
    list_sam = []
    file = open('list_ERCCcount_sam.txt', 'r')
    for i in file:
        sam_url = re.sub('\n', '', i)
        list_sam.append(sam_url)
    os.system("rm list_ERCCcount_sam.txt")

gtf = args.gtf
for i in list_sam:
    sam = i
    prefix = re.sub('Aligned.out.sam', '', sam)
    htseq_out = prefix+'_counts.txt'
    os.system("htseq-count -f sam -r name -s no -a 10 -t gene -i gene_id -m intersection-nonempty %s %s > %s" % (sam, gtf, htseq_out))

    count_out = prefix+'_count.txt'
    file_counts = open(htseq_out, 'r')
    file_count = open(count_out, 'w')
    cont_counts = file_counts.readlines()
    cont_count = cont_counts[:-5]
    for line in cont_count:
        file_count.write(line)
    file_counts.close()
    file_count.close()
