# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 16:45:13 2020

@author: chenyujiePC
"""
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import math
import csv 
import pandas as pd
from sklearn import linear_model
from scipy.special import expit
from scipy.stats import pearsonr
from matplotlib.backends.backend_pdf import PdfPages
import argparse

parser = argparse.ArgumentParser(description='生成单样本的准确率分析图: \
                                              python accuracy.py \
                                              -i FPKM.txt')

parser.add_argument('-i', '--fpkm_input', type=str, required=True, default='', help='FPKM.txt')

args = parser.parse_args()


total_mol=2.39*math.pow(10,7)
#数据读取
fpkm_input = args.fpkm_input
prefix = re.sub('_FPKM.txt', '', fpkm_input)
pdf_out = prefix+'_accuracy.pdf'
ercc=pd.read_csv( fpkm_input,sep = '\t',encoding='utf-8')["fpkm"].values
y=(np.int64(ercc>0))
ercc_tmp=pd.read_csv('ERCC_Controls_Analysis.txt',sep = '\t',encoding='utf-8')["concentration in Mix 1 (attomoles/ul)"].values
ercc_exp = ercc_tmp[:,np.newaxis]    #scipy模型中需要将X增加一维进行模拟，Y不需要变化
ercc_exp_sum=np.sum(ercc_exp)
ercc_exp_real=total_mol * ercc_exp/ercc_exp_sum
#ercc_exp_real2=total_mol * ercc_tmp/ercc_exp_sum
X=np.log(ercc_exp_real)

pdf = PdfPages(pdf_out)
plt.figure(1, figsize=(6.5, 5))
index=[i for i in range(len(ercc)) if ercc[i] > 0 ]
X=[ercc_exp_real[i] for i in index]
Y=[ercc[i] for i in index]
ln_ercc_exp_real = np.log(X)
ln_ercc = np.log(Y)
lr = linear_model.LinearRegression()
lr.fit(ln_ercc_exp_real,ln_ercc)
#lr_score = lr.score(X,Y)
#print("lr_score:{}".format(lr_score))
lr_R = pearsonr(ln_ercc_exp_real.ravel(),ln_ercc)
#print("lr_R:{}".format(lr_R))
X_prime = np.linspace(0, 15, 300)
plt.scatter(ln_ercc_exp_real,ln_ercc,s=15,color='cornflowerblue')
plt.plot(X_prime,lr.coef_ * X_prime + lr.intercept_ ,color="red",alpha=0.9)
plt.xlabel("Log input ERCC molecules")
plt.ylabel("Measureed expression (log FPKM)")
plt.annotate(r'$R = %.2f$'%lr_R[0],xy=(6,11),fontsize=12,color='blue',alpha=.8)
plt.title("Accuracy = %.2f \n Pearson correlation between input and \n quantified expression values"%lr_R[0])
ax=plt.axes()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(1.5)
ax.spines['left'].set_linewidth(1.5)
plt.tick_params(labelsize=14)
#plt.ylim(-2, 15)
#plt.xlim(-4, 10)
#plt.xlim(-1, 15)
pdf.savefig()
plt.close()
pdf.close()
