#!/usr/bin/env python
# coding: utf-8

# ============================================================================================
# ============================================================================================
# #### Figure 3.1a-b Simulation result of tree cover vs. precipitation for Hirota’s hypothesis
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Equilibrium =========================================================================
a = 0.08    
b1, b2, b3 = 20, 50, 80
c = 1900
x_eq = np.array(range(0, 100, 1))
p_eq = a*(x_eq-b1)*(x_eq-b2)*(x_eq-b3) + c  # for x'=0
" derivate to get the critical points: "
" derivative=0: "  
" => a[(x-b2)(x-b3)+(x-b1)(x-b3)+(x-b1)(x-b2)]=0 "
" => x^2-(b2+b3)x+b2*b3+x^2-(b1+b3)x+b1*b3+x^2-(b1+b2)x+b1*b2]=0 " 
" => 3x^2-2(b1+b2+b3)x+b1*b2+b2*b3+b1*b3=0 "
A = 3; B = -2*(b1+b2+b3); C = b1*b2+b2*b3+b1*b3
x_critical_1 = (-B - (B ** 2 - 4 * A * C) ** 0.5) / (2 * A)
x_critical_2 = (-B + (B ** 2 - 4 * A * C) ** 0.5) / (2 * A)
x_critical_1 = round(x_critical_1)
x_critical_2 = round(x_critical_2)
#=====================================================================================

p_range = np.arange(0, 3800, 10)  # precipitation steps
y0 = np.zeros(len(p_range))       # state initialization to start from 0%
z0 = np.full(len(p_range), 100)   # state initialization to start from 100%
mu = 75                           # noise level in random generation size for tree cover
nu = 1000                         # noise level in random generation size for precipitation

def dxdt(x, p):
    return -a*(x-b1)*(x-b2)*(x-b3) - c + p

time = np.linspace(0, 1, num=500)

y = np.zeros((len(time), len(p_range)))
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
# (a) tree cover vs. precipitation scatter
plt.figure(figsize=(5,4))
plt.scatter(precip_y[100:,], y[100:,], s=0.05, color='black')
plt.scatter(precip_z[100:,], z[100:,], s=0.05, color='black')
plt.plot(p_eq[:x_critical_1], x_eq[:x_critical_1], 'r')
plt.plot(p_eq[x_critical_1:x_critical_2], x_eq[x_critical_1:x_critical_2], 'r--')
plt.plot(p_eq[x_critical_2:], x_eq[x_critical_2:], 'r')        
plt.xlabel('Mean annual precipitation (mm/yr)', fontsize=12) # fontweight='bold'
plt.ylabel('Tree cover (%)', fontsize=12)
plt.xlim(-200, 4000)
plt.xticks([0, 1000, 2000, 3000], fontsize=10)
plt.show()

# (b) tree cover distribution
w = np.concatenate((y,z), axis=None)
plt.figure(figsize=(5,4))
sns.histplot(w, kde=True, color='black', bins=40, fill=False, stat='probability')
plt.xlabel('Tree cover (%)', fontsize=12)
plt.ylabel('Relative frequency', fontsize=12)
plt.xlim(-2, 102)
plt.show()

# ============================================================================================
# ============================================================================================
# #### Figure 3.1c Simulation result of tree cover vs. precipitation for Hirota’s hypothesis
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

np.random.seed(2)

p_mid = 1900
y0 = 0       # state initialization to start from 0%
z0 = 100     # state initialization to start from 100%
mu = 75      # noise level in random generation size for tree cover

def dxdt(x, p):
    return -0.08*(x-20)*(x-50)*(x-80) - 1900 + p

time = np.linspace(0, 1, num=500)

y = np.zeros(len(time))
y[0] = y0
for n, dt in enumerate(np.diff(time), 0):
    y[n+1] = y[n] + dxdt(y[n], p_mid)*dt + mu*np.random.normal(0, np.sqrt(dt))
        
z = np.zeros(len(time))
z[0] = z0
for n, dt in enumerate(np.diff(time), 0):
    z[n+1] = z[n] + dxdt(z[n], p_mid)*dt + mu*np.random.normal(0, np.sqrt(dt))
    
treecover = np.concatenate((y,z), axis=0)

# Graph ==============================================================================
plt.figure(figsize=(5,4))
sns.histplot(treecover, kde=True, color='black', bins=40, fill=False, stat='count')
plt.xlabel('Tree cover (%)', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.text(50, 100, 'P = {} mm/yr'.format(p_mid), horizontalalignment='center', fontsize=12)
plt.xlim(-2, 102)
plt.show()

# ============================================================================================
# ============================================================================================
# #### Figure 3.2 Simulation result of Good’s model
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import interp1d

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
#======================================================================================
def dxdt(x, b, m):
    " f(x) = beta*x*(1-x)-m*x "
    prod = b + np.random.normal(0, 0.01)
    mort = m + np.random.normal(0, m*0.15)
    return prod*x*(1-x)-mort*x, prod/mort

def score(TreeCoverArray):
    " score = fs*ff/(fdiplo*fdiphi) "
    count_tl = sum(1 for i in TreeCoverArray if i < 0.025)              # treeless
    count_diplo = sum(1 for i in TreeCoverArray if i>=0.025 and i<0.05) # btw treeless and savanna
    count_s = sum(1 for i in TreeCoverArray if i>=0.05 and i<0.3)       # savanna
    count_diphi = sum(1 for i in TreeCoverArray if i>=0.3 and i<0.7)    # btw savanna and forest
    count_f = sum(1 for i in TreeCoverArray if i>=0.7 and i<0.9)        # forest
    count_ff = sum(1 for i in TreeCoverArray if i>=0.9)                 # higher tree fraction above forest

    " score_value = count_s*count_f/(count_diplo*count_diphi) "
    if count_diplo == 0:
        score_value = count_s*count_f/(0.001*count_diphi) # not to divide by zero
    elif count_diphi == 0:
        score_value = count_s*count_f/(count_diplo*0.001) # not to divide by zero
    else:
        score_value = count_s*count_f/(count_diplo*count_diphi)
        
    return score_value
#======================================================================================
N_iteration = 1000  # nr of iterations to get the high score
N_b_range = 200

np.random.seed(22)

b_range = np.linspace(0, 1, num=8)  # beta (productivity) steps
m_range = np.random.uniform(0, 0.5, len(b_range))
interpolate = interp1d(b_range, m_range)

# N_b_range for N_iteration
b_generate = np.random.randint(0, 100, (N_b_range, N_iteration))/100
m_generate = interpolate(b_generate)
b_generate_withnoise = b_generate + np.random.normal(0, 0.01, b_generate.shape)
m_generate_withnoise = m_generate + np.random.normal(0, m_generate*0.15)  
T_generate = np.zeros(b_generate.shape)  # or with m_generate.shape

score_array = np.zeros(N_iteration)

for i in range(N_iteration):  # for each iteration
    for j, b_g in enumerate(b_generate_withnoise[:,i]):
        if b_g < m_generate_withnoise[j, i]:  # should satisfy b >= m
            T_generate[j, i] = 0
        else:
            if b_g == 0.0:
                T_generate[j, i] = 1 - m_generate_withnoise[j, i]/0.001  # avoid division by zero productivity
            else:
                T_generate[j, i] = 1 - m_generate_withnoise[j, i]/b_g
    score_array[i] = score(T_generate[:,i])

max_score_index = np.argmax(score_array)
print(score_array[max_score_index])  # print the maximum score

#======================================================================================
time = np.linspace(0, 1000, num=1000)
mu = 0.001                                               # noise level in random generation size for tree cover

T = np.zeros((len(time), N_b_range))
T0 = np.full(N_b_range, 0.5)                             # state initialization with 0.5
bm = np.zeros((len(time), N_b_range))
dTdt = np.zeros((len(time), N_b_range))
T[0] = T0

for l, b in enumerate(b_generate[:,max_score_index]):
    for n, dt in enumerate(np.diff(time), 0):
        dTdt[n,l], bm[n,l] = dxdt(T[n,l], b, m_generate[l,max_score_index])
        T[n+1,l] = T[n,l] + dTdt[n,l]*dt + mu*np.random.normal(0, np.sqrt(dt))
        
# Graphs ==============================================================================
# (a) Mortality profile
b_line = np.linspace(0, 1, num=5)
plt.figure(figsize=(3.5,3))
plt.plot(b_generate_withnoise[:,max_score_index], m_generate_withnoise[:,max_score_index], 'bo', b_range, m_range, 'b-', markersize=2)
plt.plot(b_line, (1-0.0)*b_line, label='T=0.0', color='black')
plt.plot(b_line, (1-0.2)*b_line, label='T=0.2', color='black')
plt.plot(b_line, (1-0.4)*b_line, label='T=0.4', color='black')
plt.plot(b_line, (1-0.6)*b_line, label='T=0.6', color='black')
plt.plot(b_line, (1-0.8)*b_line, label='T=0.8', color='black')
plt.text(0.8, 0.89, 'T=0.0', fontsize=10, color='red', weight='bold')
plt.text(0.8, 0.7, 'T=0.2', fontsize=10, color='red', weight='bold')
plt.text(0.8, 0.52, 'T=0.4', fontsize=10, color='red', weight='bold')
plt.text(0.8, 0.36, 'T=0.6', fontsize=10, color='red', weight='bold')
plt.text(0.8, 0.2, 'T=0.8', fontsize=10, color='red', weight='bold')
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xlabel('Productivity', fontsize=12)
plt.ylabel('Mortality', fontsize=12)
plt.show()        

# (b) Tree cover fraction estimate with beta random data
plt.figure(figsize=(3.5,3))
plt.hist(T_generate[:,max_score_index], bins=20, color = 'gray')
plt.xlim(0, 1)
plt.xlabel('Tree cover fraction estimate', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.show()

# (c) Scatter tree cover vs beta/m with Euler-Maruyama method
plt.figure(figsize=(3.5,3))
plt.scatter(bm[100:,], T[100:,], s=0.1, color='black')
plt.plot(bm_eq, x_eq, 'r')
plt.xlabel('\u03B2/m', fontsize=12)
plt.ylabel('Tree cover fraction', fontsize=12)
plt.xlim(0, 5)
plt.ylim(0, 1)
plt.show()

# (d) Tree cover distribution from (c)
T_vector=np.concatenate(T)
plt.figure(figsize=(3.5,3))
sns.histplot(np.where(T_vector<0, 0, T_vector), kde=True, color='black', bins=40, fill=False, stat='probability')
plt.xlabel('Tree cover fraction', fontsize=12)
plt.ylabel('Relative frequency', fontsize=12)
plt.xlim(-0.02, 1.02)
plt.show()
