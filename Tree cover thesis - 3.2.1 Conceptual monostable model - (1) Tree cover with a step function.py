#!/usr/bin/env python
# coding: utf-8

# (1) Tree cover with a step function:

# #### Figure 3.3 Trimodal tree cover distribution with tree cover in step function

# In[ ]:


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

k1 = 0.40 # cut point of T1
k2 = 0.60 # cut point of T2
a1l = 0.05 # lower extreme level of T1
a1h = 0.87 # higher extreme level of T1
a2l = 0.35 # lower extreme level of T2
a2h = 0.93 # higher extreme level of T2

def Teq(x1, x2):
    Teq1 = np.where(x1 < k1, a1l, a1h) # variable x1: arbitrary step function with x-axis btw 0 and 1
    Teq2 = np.where(x2 < k2, a2l, a2h) # variable x2: same as above
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

# Data ===========================================================================
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

# Graph ==========================================================================
# Plot the 3D surface
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(X1v, X2v, TEQ, cmap='viridis', linewidth=0, antialiased=False)
ax.set_zlim(0,1)
fig.colorbar(surf, shrink=0.5, aspect=5, location='left')
ax.set_xlabel('$X_1$')
ax.set_ylabel('$X_2$')
ax.set_zlabel('Tree Cover')
plt.show()

# X1 vs. X2 scatter
fig, ((ax1)) = plt.subplots(figsize=(4,3), nrows=1, ncols=1)
X1X2_scatter(ax1, dataX1, dataX2, '$X_1, X_2$~$N$({}, {}\u00b2); $X_1 \perp X_2$'.format(meanX1, stdX1))
plt.show()

# Step function
X1Range=np.linspace(0,1,100)
X2Range=np.linspace(0,1,100)
teqF, teq1, teq2 = Teq(X1Range, X2Range)
plt.figure(figsize=(4,3))
plt.plot(X1Range, teq1, color='blue', label='$T_1$')
plt.plot(X2Range, teq2, color='black', label='$T_2$')
plt.xlabel('$X_1$, $X_2$', fontsize=12)
plt.ylabel('Tree cover', fontsize=12)
plt.legend(loc='lower right')
plt.ylim(-0.05, 1.05)
plt.show()

# Tree cover distribution w/ kde
fig = plt.figure(figsize=(4,3))
sns.histplot(dataTEQ, kde=True, color='black', bins=20, fill=True, element="step", edgecolor=None)
plt.xlim(0,1)
plt.ylim(0,1000)
plt.xlabel('Tree cover fraction', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.tight_layout()
plt.show()


# #### Figure A.5 Sensitivity analysis 1

# In[ ]:


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def Teq(x1, x2, a1l):
    k1 = 0.40 # cut point of T1
    k2 = 0.60 # cut point of T2
    #a1l = 0.05 # lower extreme level of T1
    a1h = 0.90 # higher extreme level of T1
    a2l = 0.05 # lower extreme level of T2
    a2h = 0.90 # higher extreme level of T2
    Teq1 = np.where(x1 < k1, a1l, a1h) # variable x1: arbitrary step function with x-axis btw 0 and 1
    Teq2 = np.where(x2 < k2, a2l, a2h) # variable x2: same as above
    Teqf = Teq1*Teq2 
    return Teqf, Teq1, Teq2

def Step_function(ax, a1l, fontsize=12):
    X1Range=np.linspace(0,1,100)
    X2Range=np.linspace(0,1,100)
    teqF, teq1, teq2 = Teq(X1Range, X2Range, a1l)
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
    ax.set_ylim(0,1000)
    ax.set_xlabel('Tree cover fraction', fontsize=fontsize)
    ax.set_ylabel('Frequency', fontsize=fontsize)

# Graph ==========================================================================

a1l_array = [0.05, 0.25, 0.45, 0.65]
for a1l in a1l_array:
    # Make data for contour
    X1 = np.linspace(0, 1, 100)
    X2 = np.linspace(0, 1, 100)               
    X1v, X2v = np.meshgrid(X1, X2)
    TEQ, TEQ1, TEQ2 = Teq(X1v, X2v, a1l)

    # Make data for Teq vs variable (normal)
    N = 2000
    meanX1, stdX1 = 0.5, 0.3
    meanX2, stdX2 = 0.5, 0.3
    dataX1 = np.random.normal(meanX1, stdX1, N)   
    dataX2 = np.random.normal(meanX2, stdX2, N)   
    dataTEQ, dataTEQ1, dataTEQ2 = Teq(dataX1, dataX2, a1l)
    
    # Graphs
    fig, ((ax1, ax2, ax3)) = plt.subplots(figsize=(18,4), nrows=1, ncols=3)
    Step_function(ax1, a1l)
    X1X2_scatter(ax2, dataX1, dataX2, '$X_1, X_2$~$N$({}, {}\u00b2); $X_1 \perp X_2$'.format(meanX1, stdX1))
    TreecoverDistrib(ax3, dataTEQ)
    plt.show()


# #### Figure A.6 Sensitivity analysis 2

# In[ ]:


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def Teq(x1, x2, a1h):
    k1 = 0.40 # cut point of T1
    k2 = 0.60 # cut point of T2
    a1l = 0.05 # lower extreme level of T1
    #a1h = 0.90 # higher extreme level of T1
    a2l = 0.05 # lower extreme level of T2
    a2h = 0.90 # higher extreme level of T2
    Teq1 = np.where(x1 < k1, a1l, a1h) # variable x1: arbitrary step function with x-axis btw 0 and 1
    Teq2 = np.where(x2 < k2, a2l, a2h) # variable x2: same as above
    Teqf = Teq1*Teq2 
    return Teqf, Teq1, Teq2

def Step_function(ax, a1h, fontsize=12):
    X1Range=np.linspace(0,1,100)
    X2Range=np.linspace(0,1,100)
    teqF, teq1, teq2 = Teq(X1Range, X2Range, a1h)
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
    ax.set_ylim(0,1000)
    ax.set_xlabel('Tree cover fraction', fontsize=fontsize)
    ax.set_ylabel('Frequency', fontsize=fontsize)

# Graph ==========================================================================

a1h_array = [0.9, 0.7, 0.5, 0.3]
for a1h in a1h_array:
    # Make data for contour
    X1 = np.linspace(0, 1, 100)
    X2 = np.linspace(0, 1, 100)               
    X1v, X2v = np.meshgrid(X1, X2)
    TEQ, TEQ1, TEQ2 = Teq(X1v, X2v, a1h)

    # Make data for Teq vs variable (normal)
    N = 2000
    meanX1, stdX1 = 0.5, 0.3
    meanX2, stdX2 = 0.5, 0.3
    dataX1 = np.random.normal(meanX1, stdX1, N)   
    dataX2 = np.random.normal(meanX2, stdX2, N)   
    dataTEQ, dataTEQ1, dataTEQ2 = Teq(dataX1, dataX2, a1h)
    
    # Graphs
    fig, ((ax1, ax2, ax3)) = plt.subplots(figsize=(18,4), nrows=1, ncols=3)
    Step_function(ax1, a1h)
    X1X2_scatter(ax2, dataX1, dataX2, '$X_1, X_2$~$N$({}, {}\u00b2); $X_1 \perp X_2$'.format(meanX1, stdX1))
    TreecoverDistrib(ax3, dataTEQ)
    plt.show()


# #### Figure A.7 Sensitivity analysis 3

# In[ ]:


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def Teq(x1, x2, a1l):
    k1 = 0.40 # cut point of T1
    k2 = 0.60 # cut point of T2
    #a1l = 0.05 # lower extreme level of T1
    a1h = 0.90 # higher extreme level of T1
    a2l = a1l # lower extreme level of T2
    a2h = 0.90 # higher extreme level of T2
    Teq1 = np.where(x1 < k1, a1l, a1h) # variable x1: arbitrary step function with x-axis btw 0 and 1
    Teq2 = np.where(x2 < k2, a2l, a2h) # variable x2: same as above
    Teqf = Teq1*Teq2 
    return Teqf, Teq1, Teq2

def Step_function(ax, a1l, fontsize=12):
    X1Range=np.linspace(0,1,100)
    X2Range=np.linspace(0,1,100)
    teqF, teq1, teq2 = Teq(X1Range, X2Range, a1l)
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
    ax.set_ylim(0,1000)
    ax.set_xlabel('Tree cover fraction', fontsize=fontsize)
    ax.set_ylabel('Frequency', fontsize=fontsize)

# Graph ==========================================================================

a1l_array = [0.05, 0.1, 0.2, 0.3]
for a1l in a1l_array:
    # Make data for contour
    X1 = np.linspace(0, 1, 100)
    X2 = np.linspace(0, 1, 100)               
    X1v, X2v = np.meshgrid(X1, X2)
    TEQ, TEQ1, TEQ2 = Teq(X1v, X2v, a1l)

    # Make data for Teq vs variable (normal)
    N = 2000
    meanX1, stdX1 = 0.5, 0.3
    meanX2, stdX2 = 0.5, 0.3
    dataX1 = np.random.normal(meanX1, stdX1, N)   
    dataX2 = np.random.normal(meanX2, stdX2, N)   
    dataTEQ, dataTEQ1, dataTEQ2 = Teq(dataX1, dataX2, a1l)
    
    # Graphs
    fig, ((ax1, ax2, ax3)) = plt.subplots(figsize=(18,4), nrows=1, ncols=3)
    Step_function(ax1, a1l)
    X1X2_scatter(ax2, dataX1, dataX2, '$X_1, X_2$~$N$({}, {}\u00b2); $X_1 \perp X_2$'.format(meanX1, stdX1))
    TreecoverDistrib(ax3, dataTEQ)
    plt.show()

