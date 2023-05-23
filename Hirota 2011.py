#!/usr/bin/env python
# coding: utf-8

# ### This is the final picture of simulation of Hirota 2011:
import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline

# Equilibrium =========================================================================
a = 0.1
b1 = 20; b2 = 50; b3 = 80
c = 2000
x_eq = np.array(range(0, 100, 1))
p_eq = a*(x_eq-b1)*(x_eq-b2)*(x_eq-b3) + c  # for x'=0
# derivate to get the critical points:
# derivative=0: a[(x-b2)(x-b3)+(x-b1)(x-b3)+(x-b1)(x-b2)]=0 => x^2-(b2+b3)x+b2*b3+x^2-(b1+b3)x+b1*b3+x^2-(b1+b2)x+b1*b2]=0
# => 3x^2-2(b1+b2+b3)x+b1*b2+b2*b3+b1*b3=0
A = 3; B = -2*(b1+b2+b3); C = b1*b2+b2*b3+b1*b3
x_critical_1 = (-B - (B ** 2 - 4 * A * C) ** 0.5) / (2 * A)
x_critical_2 = (-B + (B ** 2 - 4 * A * C) ** 0.5) / (2 * A)
x_critical_1 = round(x_critical_1)
x_critical_2 = round(x_critical_2)
#=====================================================================================

p_range = np.arange(0, 4000, 10)  # precipitation steps
y0 = np.zeros(len(p_range))       # state initialization to start from zero
z0 = np.full(len(p_range), 100)   # state initialization to start from 100
mu = 120                          # noise level in random generation size for tree cover
nu = 2500                         # noise level in random generation size for precipitation

def dxdt(x, p):
    return -0.2*(x-20)*(x-50)*(x-80) - 2000 + p

time = np.linspace(0, 1, num=1000)

y = np.zeros((len(time), len(p_range)))  # 1000 x 400
precip_y = np.zeros((len(time), len(p_range)))
y[0] = y0
for i, p in enumerate(p_range):
    for n, dt in enumerate(np.diff(time), 0):
        precip_y[n+1,i] = p + nu*np.random.normal(0, np.sqrt(time[-1]/time.size))
        y[n+1,i] = y[n,i] + dxdt(y[n,i], precip_y[n+1,i])*dt + mu*np.random.normal(0, np.sqrt(dt))
        
z = np.zeros((len(time), len(p_range)))
precip_z = np.zeros((len(time), len(p_range)))
z[0] = z0
for i, p in enumerate(p_range):
    for n, dt in enumerate(np.diff(time), 0):
        precip_z[n+1,i] = p + nu*np.random.normal(0, np.sqrt(time[-1]/time.size))
        z[n+1,i] = z[n,i] + dxdt(z[n,i], precip_z[n+1,i])*dt + mu*np.random.normal(0, np.sqrt(dt))
        
# Graph ==============================================================================
plt.figure(figsize=(8,6))
plt.scatter(precip_y[800:,], y[800:,], s=0.1, color='black')
plt.scatter(precip_z[800:,], z[800:,], s=0.1, color='black')
plt.plot(p_eq[:x_critical_1], x_eq[:x_critical_1], 'r')
plt.plot(p_eq[x_critical_1:x_critical_2], x_eq[x_critical_1:x_critical_2], 'r--')
plt.plot(p_eq[x_critical_2:], x_eq[x_critical_2:], 'r')        
plt.xlabel('Mean annual precipitation (mm/yr)', fontweight='bold')
plt.ylabel('Tree cover (%)', fontweight='bold')
plt.xlim(-200, 4200)
plt.ylim(-5, 105)
plt.xticks([0, 1000, 2000, 3000, 4000])
plt.title('Tree cover - Precipitation')
plt.show()


# ### Plot histogram

# In[24]:


# plot histogram
w = np.concatenate((y, z))
plt.hist(w[200:,20], bins=40)
plt.xlabel('Tree Cover')
plt.ylabel('Frequency')
plt.show()


# In[261]:


w.shape


# ### Autocorrelation graph

# In[305]:


from statsmodels.graphics.tsaplots import plot_acf

# plot autocorrelation
plt.rc("figure", figsize=(11,5))
plot_acf(y[:,1])
plt.show()
# it seems there is no autocorrelation for both time series starting from zero and those starting from 100


# ### Augmented Dickey-Fuller test
# - To test if the stochastic process (time series) is stationary
# - Checking this because for the stochastic process to be ergodic, it's implicitly assumed that the process is stationary

# In[310]:


from statsmodels.tsa.stattools import adfuller
result = adfuller(w[:, 0])  # Augmented Dickey-Fuller test
result
# ADF Statistic: -4.311
# p-value: 0.000425 <= 0.05
# => -4.311 < -3.434 at 1%  
# => reject the null hypothesis (process has no unit root => time series is stationary or doesn't have time-dependent structure. No trend?)


# - it seems like that around p=2000mm, the time series combined y and z is not stationary => not ergodic(?)
# - but far from the middle of p=2000mm, for example p=100mm or p=4000mm, the time series is stationary => can be ergodic.
# - if we check y and z independently, their individual times series is stationary => can be ergodic.

# ### Check space for time (ergodicity)

# In[61]:


import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

np.random.seed(2)

p_mid = 2000
y0 = 0       # state initialization to start from zero
z0 = 100     # state initialization to start from 100
mu = 120     # noise level in random generation size for tree cover

def dxdt(x, p):
    return -0.2*(x-20)*(x-50)*(x-80) - 2000 + p

time = np.linspace(0, 1, num=1000)

y = np.zeros(len(time))
y[0] = y0
for n, dt in enumerate(np.diff(time), 0):
    y[n+1] = y[n] + dxdt(y[n], p_mid)*dt + mu*np.random.normal(0, np.sqrt(dt))
        
z = np.zeros(len(time))
z[0] = z0
for n, dt in enumerate(np.diff(time), 0):
    z[n+1] = z[n] + dxdt(z[n], p_mid)*dt + mu*np.random.normal(0, np.sqrt(dt))
    
treecover = np.concatenate((y,z), axis=0)

ax = sns.displot(data=treecover, kde=True, bins=40, height=3, aspect=2) # aspect: height will be 2 times width
ax.set(xlabel='Tree fraction estimate', ylabel='Count', title='Tree Cover PDF for p=2000mm', xlim=(0, 100), ylim=(0, 300))
plt.show()

# Graph ==============================================================================
#plt.figure(figsize=(8,6))
#plt.hist(treecover, bins=20, density=True, alpha=0.5)
#plt.xlabel('Tree fraction estimate', fontweight='bold')
#plt.ylabel('Frequency', fontweight='bold')
#plt.xlim(0, 100)
#plt.show()

