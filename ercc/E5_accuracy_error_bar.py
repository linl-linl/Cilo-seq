# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 22:08:50 2020

@author: chenyujiePC
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import csv 
import pandas as pd
from sklearn import linear_model
from scipy.special import expit
from scipy.stats import pearsonr
from matplotlib.backends.backend_pdf import PdfPages

pdf = PdfPages('accuracy_error_bar3_fpkm.pdf')
plt.figure(1, figsize=(6, 6))

fpkm2=pd.read_csv( '6-0127ZQQ2_L1_FPKM.txt',sep = '\t',encoding='utf-8')["fpkm"].values
fpkm3=pd.read_csv( '7-0302ZQQ2_L1_FPKM.txt',sep = '\t',encoding='utf-8')["fpkm"].values
fpkm4=pd.read_csv( '8-0312ZQQ1_L1_FPKM.txt',sep = '\t',encoding='utf-8')["fpkm"].values

index=[i for i in range(len(fpkm2)) if fpkm2[i] > 0 and fpkm3[i] > 0 and fpkm4[i] > 0]

fpkm2_=[fpkm2[i] for i in index]
fpkm3_=[fpkm3[i] for i in index]
fpkm4_=[fpkm4[i] for i in index]

a=np.array([np.log(fpkm2_),np.log(fpkm3_),np.log(fpkm4_)])
std=np.std(a, axis=0)

b = np.array([fpkm2_, fpkm3_, fpkm4_])
mean=np.mean(b,axis=0)
Y=np.log(mean)

total_mol=2.39*math.pow(10,7)
# ercc_tmp=pd.read_csv('ERCC_Controls_Analysis.txt',sep = '\t',encoding='utf-8')["concentration in Mix 1 (attomoles/ul)"].values
# ercc_tmp_sum=np.sum(ercc_tmp)
# ercc_exp_real=total_mol * ercc_tmp/ercc_tmp_sum
ercc_exp_real=pd.read_csv('ercc_fpkm.txt',sep = '\t',encoding='utf-8')["fpkm"].values
ercc_=[ercc_exp_real[i] for i in index]
X=np.log(ercc_)

lr_R = pearsonr(X,Y)

plt.scatter(X,Y,s=3)
plt.errorbar(X, Y, yerr=std, fmt="o",color="royalblue",ecolor='grey',elinewidth=1,capsize=3)

plt.xlabel("Log input ERCC molecules (FPKM)")
plt.ylabel("Measureed expression (logFPKM)")
plt.title("R = %.2f "%lr_R[0])
ax=plt.axes()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(1.5)
ax.spines['left'].set_linewidth(1.5)
#ax.spines['bottom'].set_position(('data',0))
#ax.spines['left'].set_position(('data',0))
plt.tick_params(labelsize=14)
plt.ylim(-3, 16)
plt.xlim(-3, 16)
plt.xticks([-2,0,2,4,6,8,10,12,14])
plt.yticks([-2,0,2,4,6,8,10,12,14])
pdf.savefig()
plt.close()
pdf.close()

