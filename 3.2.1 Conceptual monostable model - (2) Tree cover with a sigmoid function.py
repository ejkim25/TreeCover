#!/usr/bin/env python
# coding: utf-8

# (2) Tree cover with a sigmoid function

# ====================================================================================
# #### Figure 3.4 Multimodal tree cover distribution with tree cover function with sigmoid curve

# In[ ]:

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

b11, b12, b13, b14 = 0.88, 5, 16, 0.05
b21, b22, b23, b24 = 0.62, 5, 11, 0.35

def Teq(x1, x2):
    Teq1 = b11/(1.0+np.exp(b12-b13*x1))+b14
    Teq2 = b21/(1.0+np.exp(b22-b23*x2))+b24
    Teqf = Teq1*Teq2  
    return Teqf, Teq1, Teq2
     
def X1X2_scatter(ax, x1, x2, title, fontsize=11):
    CT=ax.contourf(X1v, X2v, TEQ, 20)
    ax.scatter(x1, x2, s=1, color='black')
    ax.set_xlabel('$X_1$', fontsize=fontsize)
    ax.set_ylabel('$X_2$', fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    plt.colorbar(CT, location='right', label='Tree Cover [%]')

# Data -------------------------------------------------------------------------    
# Make data for contour
X1 = np.linspace(0, 1, 100)
X2 = np.linspace(0, 1, 100)               
X1v, X2v = np.meshgrid(X1, X2)
TEQ, TEQ1, TEQ2 = Teq(X1v, X2v)

# Make data for Teq vs variable (normal)
N = 2000
meanX1, stdX1 = 0.5, 0.3
meanX2, stdX2 = 0.5, 0.3
dataX1 = np.random.normal(meanX1, stdX1, N)   
dataX2 = np.random.normal(meanX2, stdX2, N)
dataTEQ, dataTEQ1, dataTEQ2 = Teq(dataX1, dataX2)

# Graph --------------------------------------------------------------------------
# X1 vs. X2 scatter
fig, ((ax1)) = plt.subplots(figsize=(4,3), nrows=1, ncols=1)
X1X2_scatter(ax1, dataX1, dataX2, '$X_1, X_2$~$N$({}, {}\u00b2); $X_1 \perp X_2$'.format(meanX1, stdX1))
plt.show()

# sigmoid function
X1Range=np.linspace(0,1,100)
X2Range=np.linspace(0,1,100)
TCF_fct, TC1_fct, TC2_fct = Teq(X1Range, X2Range)
plt.figure(figsize=(4,3))
plt.plot(X1Range, TC1_fct, color='blue', label='$T_1$')
plt.plot(X2Range, TC2_fct, color='black', label='$T_2$')
plt.xlabel('$X_1$, $X_2$', fontsize=12)
plt.ylabel('Tree cover', fontsize=12)
plt.legend(loc='lower right')
plt.ylim(-0.05, 1.05)
plt.show()

# Tree cover distribution w/ kde
fig = plt.figure(figsize=(4,3))
sns.histplot(dataTEQ, kde=True, color='black', bins=20, fill=True, element="step", edgecolor=None)
plt.xlim(0,1)
plt.ylim(0,300)
plt.xlabel('Tree cover fraction', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.tight_layout()
plt.show()

# ====================================================================================
# #### Figure A.8 Sensitivity analysis

# In[ ]:

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

b11, b12, b14 = 0.88, 5, 0.05
b21, b22, b23, b24 = 0.88, 5, 11, 0.05

def Teq(x1, x2, b13):
    Teq1 = b11/(1.0+np.exp(b12-b13*x1))+b14
    Teq2 = b21/(1.0+np.exp(b22-b23*x2))+b24
    Teqf = Teq1*Teq2  
    return Teqf, Teq1, Teq2

def SigmoidFct(ax, b13, fontsize=12):
    X1Range=np.linspace(0,1,100)
    X2Range=np.linspace(0,1,100)
    teqF, teq1, teq2 = Teq(X1Range, X2Range, b13)
    ax.plot(X1Range, teq1, color='blue', label='$T_1$')
    ax.plot(X2Range, teq2, color='black', label='$T_2$')
    ax.set_xlabel('$X_1$, $X_2$', fontsize=fontsize)
    ax.set_ylabel('Tree cover', fontsize=fontsize)
    ax.set_ylim(-0.05, 1.05)
    ax.legend(loc='lower right')   

def X1X2_scatter(ax, x1, x2, title, fontsize=11):
    CT=ax.contourf(X1v, X2v, TEQ, 20)
    ax.scatter(x1, x2, s=1, color='black')
    ax.set_xlabel('$X_1$', fontsize=fontsize)
    ax.set_ylabel('$X_2$', fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    plt.colorbar(CT, location='right', label='Tree Cover [%]')

def TreecoverDistrib(ax, dataTEQ, fontsize=12):
    sns.histplot(dataTEQ, kde=True, color='black', bins=20, fill=True, element="step", edgecolor=None)
    ax.set_xlim(0,1)
    ax.set_ylim(0,400)
    ax.set_xlabel('Tree cover fraction', fontsize=fontsize)
    ax.set_ylabel('Frequency', fontsize=fontsize)

# Graph ----------------------------------------------------------------------------

b13_array = [26, 16, 10, 6]
for b13 in b13_array:
    # Make data for contour
    X1 = np.linspace(0, 1, 100)
    X2 = np.linspace(0, 1, 100)               
    X1v, X2v = np.meshgrid(X1, X2)
    TEQ, TEQ1, TEQ2 = Teq(X1v, X2v, b13)

    # Make data for Teq vs variable (normal)
    N = 2000
    meanX1, stdX1 = 0.5, 0.3
    meanX2, stdX2 = 0.5, 0.3
    dataX1 = np.random.normal(meanX1, stdX1, N)   
    dataX2 = np.random.normal(meanX2, stdX2, N)   
    dataTEQ, dataTEQ1, dataTEQ2 = Teq(dataX1, dataX2, b13)
    
    # Graphs
    fig, ((ax1, ax2, ax3)) = plt.subplots(figsize=(18,4), nrows=1, ncols=3)
    SigmoidFct(ax1, b13)
    X1X2_scatter(ax2, dataX1, dataX2, '$X_1, X_2$~$N$({}, {}\u00b2); $X_1 \perp X_2$'.format(meanX1, stdX1))
    TreecoverDistrib(ax3, dataTEQ)
    plt.show()

