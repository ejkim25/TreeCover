#!/usr/bin/env python
# coding: utf-8

# ## Mortality influenced by low wood density, high precipitation, high or low temperature, tree type

# - Productivity in function of Precipitation, Temperature, Soil Nutrient
#     - Soil Nutrient: high soil nutrients -> higher productivity, but soil nutrients are washed away by strong precipitation in tropical rainforest
#     - Precipitation: sigmoid-shaped function together with uniformly distributed Precipitation
# - Mortality influenced by low wood density (let's use low tree cover), high precipitation, high temperature, and tree type
#     - Mortality(tree cover) by wood density (horizontally reflected sigmoidal function): higher wood density, higher survival probability
#     - Mortality(precipitation): low or high precipitation increases the mortality
#     - Mortality(temperature): low or high temperature increases the mortality (26<T<29.5oC is ideal for the tree)

# ### Uniform scatter of Precipitation and soil nutrient depending on the precipitation

# In[150]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator

# Equilibrium =========================================================================
" f(x) = (b-m)x-bx^2 "
" f(x) = 0 => x=0 or x=1-m/b=1-1/bm (bm: b/m) "
" if b>=m, x=0 unstable and x=1-1/bm stable " 
" df/dx=0 => b-m-2bx = 0 "
bm_eq = np.arange(0.000001, 5, 0.1)
x_eq = np.zeros(len(bm_eq))
for i, BM in enumerate(bm_eq):
    if BM < 1:
        x_eq[i] = 0
    elif BM >= 1:
        x_eq[i] = 1 - 1/BM

# Definition ==========================================================================
# parameters for productivity(precipitation) equation:
alpha = 0.000054           # original from Bathiany 2013: 0.0011
c = 138                    # original from Bathiany 2013: 28
gamma = 0.00017*9100

P1 = c*np.exp(gamma/2) # =~300                                       # original: 60.68548303127268
P2 = c*np.exp(gamma/2) + np.exp(gamma)/np.sqrt(0.03*alpha) # =~4000  # original: 878.3903705049846

def productivityPrec(p): # productivity (beta) in function of temperature "sine" function and of precipitation (P) - "sigmoid" equation                 
    if p < P1:
        betaP = 0
    elif P1 <= p <= P2:
        betaP = 1.0 - 1.0/(1+alpha*((p-P1)/np.exp(gamma))**2)
    else:
        betaP = 1.0
    return betaP

def productivitySN(sn):                            
    #betaSN = sn          # productivity to be equal to soil nutrients  
    betaSN = 1/(1+np.exp(-10*sn+5))  # sigmoidal shape of productivity
    return betaSN

def productivityTemp(temp):       # productivity (beta) in function of temperature "sine" function and of precipitation (P) - "sigmoid" equation                     
    betaT = np.sin((temp/56)*np.pi)  # contribution of temperature to productivity; gives the hill shape function with x-axis 0-70
    return betaT

def mortalityPrec(p):
    if p <= 2000:
        mortP = 0.25
    elif 2000 < p <= 3500:
        mortP = 0.009
    elif p > 3500:
        mortP = 0.25
        #mortP = (0.85/1500)*p-1.26
    return mortP

def mortalityTemp(tp):
    if tp <= 26:
        mortT = 0.35
    elif 26 < tp <= 29.5:
        mortT = 0.009
    elif tp > 29.5:
        mortT = 0.35
    return mortT

def mortalityWoodDensity(x):
    mortWD = 0.01/(1+np.exp(10*x-5))
    return mortWD

def dxdt(init, x, Prec):
    " f(x) = beta*x*(1-x)-m*x "
    u, v, w = 1/3, 1/3, 1/3   # weight                      
    k = -1/4000
    nu = 1                                                  
    sn = 1 + k*Prec + nu*np.random.normal(0, 0.05)
    Temp = 5 + 0.02*Prec + np.random.normal(0, 1)

    prod = u*productivityPrec(Prec) + v*productivitySN(sn) + w*productivityTemp(Temp) + np.random.normal(0, 0.05)
    if init == 1:
        mort = 0.4*mortalityPrec(Prec) + 0.5*mortalityWoodDensity(x) + 0.1*mortalityTemp(Temp) + np.random.normal(0, 0.05)
        # when high tree cover, the precipitation is more important factor for the mortality?
    else:
        mort = 0.05*mortalityPrec(Prec) + 0.05*mortalityWoodDensity(x) + 0.9*mortalityTemp(Temp) + np.random.normal(0, 0.05)
        # when low tree cover, the temperature is more important factor for the mortality?
    return prod*x*(1-x)-mort*x, prod/mort, prod, sn

#======================================================================================
np.random.seed(22)
#======================================================================================
time = np.linspace(0, 1000, num=10000)
mu = 0.05  # noise level in random generation size for tree cover

# Tree cover generation from T0=0.0 and T0=1.0 ========================================
dTdt = np.zeros((len(time),2))
bm = np.zeros((len(time),2))
PREC = np.random.uniform(0, 4000, (len(time),2))               # generate random precipitation with mean=2000, std=600
TEMP = np.zeros((len(time),2))
#TEMP = np.random.uniform(20, 45, (len(time),2))
SN = np.zeros((len(time), 2))

b = np.zeros((len(time),2))
T = np.zeros((len(time),2))

T[0] = [0.0, 1.0]   # state initialization

for i in range(2):
    for n, dt in enumerate(np.diff(time), 0):
        dTdt[n,i], bm[n,i], b[n,i], SN[n,i] = dxdt(i, T[n,i], PREC[n,i])
        T[n+1,i] = T[n,i] + dTdt[n,i]*dt + mu*np.random.normal(0, np.sqrt(dt))

# Graph ===============================================================================    
# Time series
plt.figure(figsize=(4,3))
plt.plot(T[:,0], color='black')
plt.plot(T[:,1], color='blue')
plt.ylim(0,1)
plt.title('Tree cover time series', fontsize=14)
plt.xlabel('Time')
plt.ylabel('Tree Cover') 
plt.show()

# Scatter: Precipitation vs Soil Nutrient
plt.figure(figsize=(4,3))
plt.scatter(SN, PREC, color='black', s=0.1)
plt.xlim(0,1)
plt.ylim(0,4000)
plt.title('Precipitation vs Soil Nutrient', fontsize=14)
plt.xlabel('Soil Nutrient')
plt.ylabel('Precipitation [mm/yr]')
plt.show()

# Graph Tree cover - Productivity/Mortality
plt.figure(figsize=(4,3))
plt.plot(bm_eq, x_eq, color='red')
plt.scatter(bm[1000:,], T[1000:,], s=0.1, color='black')
plt.xlim(0,5)
plt.ylim(0,1)
plt.title('Tree cover - Productivity/Mortality', fontsize=14)
plt.xlabel('\u03B2/m')
plt.ylabel('Tree Cover') 
plt.show()

# Productivity histogram
plt.figure(figsize=(4,3))
plt.hist(b[1000:,], bins=40)
plt.xlim(0,1)
#plt.ylim(0,100)
plt.title('Productivity frequency', fontsize=14)
plt.xlabel('Productivity')
plt.ylabel('Frequency')
plt.show()

# Tree fraction histogram
plt.figure(figsize=(4,3))
plt.hist(T[1000:,], bins=40)
plt.xlim(0,1)
#plt.ylim(0,100)
plt.title('Tree cover frequency', fontsize=14)
plt.xlabel('Tree fraction estimate')
plt.ylabel('Frequency')
plt.show()

# 3D Scatter: Tree cover - Precipitation - Soil Moisture
fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(SN[1000:], PREC[1000:], T[1000:], s=0.1, color='black')
ax.set_xlabel('Soil Nutrient')
ax.set_ylabel('Precipitation [mm/yr]')
ax.set_zlabel('Tree Cover')
ax.set_xlim(0, 1)
ax.set_ylim(0, 4000)
ax.set_zlim(0, 1.0)
plt.show()

# 3D Scatter: Tree cover - Soil Moisture - Precipitation
fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(PREC[1000:], SN[1000:], T[1000:], s=0.1, color='black')
ax.set_xlabel('Precipitation [mm/yr]')
ax.set_ylabel('Soil Nutrient')
ax.set_zlabel('Tree Cover')
ax.set_xlim(0, 4000)
ax.set_ylim(0, 1)
ax.set_zlim(0, 1.0)
plt.show()

