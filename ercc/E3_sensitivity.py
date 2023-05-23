#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 22:42:50 2020

@author: rwd
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

parser = argparse.ArgumentParser(description='生成单样本的灵敏度分析图: \
                                              python sensitivity.py \
                                              -i fpkm_input')

parser.add_argument('-i', '--fpkm_input', type=str, required=True, help='[sample_name]_FPKM.txt')

args = parser.parse_args()


total_mol=2.39*math.pow(10,7)
#数据读取
fpkm_input = args.fpkm_input
prefix = re.sub('_FPKM.txt', '', fpkm_input)
pdf_out = prefix+'_sensitivity.pdf'
ercc=pd.read_csv(fpkm_input, sep = '\t',encoding='utf-8')["fpkm"].values
y=(np.int64(ercc>0))
ercc_tmp=pd.read_csv('ERCC_Controls_Analysis.txt',sep = '\t',encoding='utf-8')["concentration in Mix 1 (attomoles/ul)"].values
ercc_exp = ercc_tmp[:,np.newaxis]    #scipy模型中需要将X增加一维进行模拟，Y不需要变化
ercc_exp_sum=np.sum(ercc_exp)
ercc_exp_real=total_mol * ercc_exp/ercc_exp_sum
#ercc_exp_real2=total_mol * ercc_tmp/ercc_exp_sum
X=np.log(ercc_exp_real)
# Fit the classifier
clf = linear_model.LogisticRegression(C=1e5)
clf.fit(X, y)
X_test = np.linspace(-5, 15, 300)
loss = expit(X_test * clf.coef_ + clf.intercept_).ravel()  #画图时不需要多一维因此进行降维处理
#ravel() matrix to vector
#'Logistic Regression Model'
pdf = PdfPages(pdf_out)
plt.figure(1, figsize=(7,6))
y_jittered = [i +np.random.uniform(-0.05,0.05) for i in y ]
index_1=[i for i in range(len(y_jittered)) if y_jittered[i] > 0.5 ]
index_0=[i for i in range(len(y_jittered)) if y_jittered[i] < 0.5 ]
Y_1=[y_jittered[i] for i in index_1 ]
X_1=[X.ravel()[i] for i in index_1 ]
Y_0=[y_jittered[i] for i in index_0 ]
X_0=[X.ravel()[i] for i in index_0 ]
P05=(-clf.intercept_ / clf.coef_)
mol_lim=np.exp(P05)
#plt.scatter(X.ravel(), y_jittered,s=10,color='cornflowerblue')
plt.scatter(X_1,Y_1,s=15,color='cornflowerblue')
plt.scatter(X_0,Y_0,s=15,color='hotpink')
plt.plot(X_test, loss, linewidth=2.3)
plt.vlines(x=P05,ymin=-1,ymax=0.5,color="grey",linestyle="--")
plt.hlines(y=0.5,xmin=-10,xmax=P05,color="grey",linestyle="--")
#plt.axhline(y=0.5,xmax=2,color="grey",linestyle="--")
plt.xlabel("Log input ERCC molecules")
plt.ylabel("Detection probability")
plt.title("Sensitiyity = %.2f\n Input level with detection probability >0.5"%mol_lim)
plt.xticks(range(-4, 17,2))
plt.yticks([0, 0.5, 1])
plt.ylim(-.1, 1.1)
#plt.xlim(-4, 10)
plt.xlim(-6, 17)
ax=plt.axes()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_linewidth(1.5)
ax.spines['left'].set_linewidth(1.5)
plt.tick_params(labelsize=12)
#plt.tight_layout()
#plt.show()
pdf.savefig()
pdf.close()
plt.close()





