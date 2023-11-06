#!/usr/bin/env python
# coding: utf-8
# ======================================================================================
# Case 2: two variables with unimodal input and one variable with bimodal input
# ======================================================================================

# ======================================================================================
# #### Figure 3.7 Outcome of Case 2 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks

np.random.seed(123)

# Equilibrium -------------------------------------------------------------
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
# -------------------------------------------------------------------------
def productivityTemp(temp):                                  
    #betaT = np.sin((temp/15)*np.pi-4.2)
    betaT = -0.0175*(temp-20)*(temp-35)
    return betaT

def productivityTheta(theta):       
    betaTheta = theta                
    return betaTheta

def productivityPrec(p):                        
    alpha = 0.000054           
    c = 138                    
    gamma = 0.00017*9100

    P1 = c*np.exp(gamma/2) # =~300                                       
    P2 = c*np.exp(gamma/2) + np.exp(gamma)/np.sqrt(0.03*alpha) # =~4000  
    
    betaP = np.where(p < P1, 0.0, np.where(p <= P2, 1.0 - 1.0/(1+alpha*((p-P1)/np.exp(gamma))**2) , 1.0))

    return betaP

def Teq(temp, theta, prec):
    mort = 0.12
    prod = productivityTemp(temp)*productivityTheta(theta)*productivityPrec(prec) #+ np.random.normal(0, 0.05)
    BM = prod/mort
    
    T_eq = np.where(BM < 1, 0, 1-1/BM)  # "if" statement for matrix

    return T_eq, BM

def PrecTemp_scatter(ax, temp, prec, TEQv, title, text, fontsize=10):
    CT = ax.contourf(TEMPv, PRECv, TEQv, 20) 
    ax.scatter(temp, prec, s=1, color='black')
    ax.set_xlabel('Temperature [°C]', fontsize=fontsize)
    ax.set_ylabel('Precipitation [mm/yr]', fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)
    AT=AnchoredText(text, loc=4, prop=dict(size=8))
    AT.patch.set_alpha(0.5)
    ax.add_artist(AT)
    plt.colorbar(CT, location='right', label='Tree Cover [%]')
    ax.set_xlim(20,35)
    ax.set_ylim(0,3000)

def PrecTheta_scatter(ax, theta, prec, TEQv, title, text, fontsize=10):
    CT2 = ax.contourf(TTv2, PRECv2, TEQv, 20) 
    ax.scatter(theta, prec, s=1, color='black')
    ax.set_xlabel(r'$\theta$', fontsize=fontsize)
    ax.set_ylabel('Precipitation [mm/yr]', fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)
    AT=AnchoredText(text, loc=4, prop=dict(size=8))
    AT.patch.set_alpha(0.5)
    ax.add_artist(AT)
    plt.colorbar(CT2, location='right', label='Tree Cover [%]')
    ax.set_xlim(0,1)
    ax.set_ylim(0,3000)
    
def TreeCover_hist(ax, teq_lfc, teq_hfc, title, fontsize=10):
    ax.hist(teq_lfc, bins=20, alpha=0.5, label = r'$\theta_{Low}$', color = 'grey')
    ax.hist(teq_hfc, bins=20, alpha=0.6, label = r'$\theta_{High}$', color = 'lightskyblue')
    ax.set_xlabel('Tree fraction estimate', fontsize=fontsize)
    ax.set_ylabel('Tree Cover frequency', fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)
    ax.set_xlim(0,1)
    ax.set_ylim(0,500)
    ax.legend(loc='upper center')

# Data -------------------------------------------------------------------------
# Make data for contour 1 (X: temperature; Y: precipitation)
TEMP = np.linspace(20, 35, 100)
PREC = np.linspace(0, 4000, 100)               
LTT = np.full((TEMP.size, PREC.size), 0.25) # lower range of theta
HTT = np.full((TEMP.size, PREC.size), 0.75) # higher range of theta
TEMPv, PRECv = np.meshgrid(TEMP, PREC)
TEQ_ltt, BM_ltt = Teq(TEMPv, LTT, PRECv)
TEQ_htt, BM_htt = Teq(TEMPv, HTT, PRECv)
TEColor= np.random.uniform(-0.01, 1.01, (TEMP.size, PREC.size)) # fake data for color scale

# Make data for contour 2 (X: theta; Y: precipitation)
TT2 = np.linspace(0, 1, 100)
PREC2 = np.linspace(0, 4000, 100)
TEMP2 = np.full((TT2.size, PREC2.size), 28) #---------------------------------> can change temperature
TTv2, PRECv2 = np.meshgrid(TT2, PREC2)
TEQ2, BM2 = Teq(TEMP2, TTv2, PRECv2)

# 1) Make random data for Teq vs productivity (Precip: normal; Temp: normal; FC: uniform; P-T indep)
datasize1 = 1500
meanT1, stdT1 = 27.5, 3
meanP1, stdP1 = 1500, 500
LTT_beta1 = np.random.uniform(0.2, 0.3, datasize1) # lower range of theta
HTT_beta1 = np.random.uniform(0.7, 0.8, datasize1) # higher range of theta
TT_beta1 = np.concatenate((LTT_beta1, HTT_beta1), axis=0)
TEMP_beta1 = np.random.normal(meanT1, stdT1, datasize1)  
TEMP_beta1_double = np.random.normal(meanT1, stdT1, datasize1*2)
PREC_beta1 = np.random.normal(meanP1, stdP1, datasize1) 
PREC_beta1_double = np.random.normal(meanP1, stdP1, datasize1*2)
TEQ_LTT1, BM_LTT1 = Teq(TEMP_beta1, LTT_beta1, PREC_beta1)
TEQ_HTT1, BM_HTT1 = Teq(TEMP_beta1, HTT_beta1, PREC_beta1)

data1 = np.concatenate((TEQ_LTT1, TEQ_HTT1), axis=None)
kde1 = gaussian_kde(data1)
x_grid1 = np.linspace(-0.1, 1.1, 1000)
density_estimation1 = kde1(x_grid1)
peaks1, _ = find_peaks(density_estimation1, height=0.5)
num_modes1 = len(peaks1)
print("Number of modes for case 1):", num_modes1)

# 2) Make random data for Teq vs productivity (Precip: normal; Temp: normal; FC: uniform; P-T corr)
datasize2 = 1500
meanT2, stdT2 = 27.5, 3
stdEpsP2 = 300
LTT_beta2 = np.random.uniform(0.2, 0.3, datasize2) # lower range of theta
HTT_beta2 = np.random.uniform(0.7, 0.8, datasize2) # higher range of theta
TT_beta2 = np.concatenate((LTT_beta2, HTT_beta2), axis=0)
TEMP_beta2 = np.random.normal(meanT2, stdT2, datasize2)  
TEMP_beta2_double = np.random.normal(meanT2, stdT2, datasize2*2)
slope2 = (2200-500)/(35-20) # slope btw temp and prec
cst2 = 500-slope2*20
PREC_beta2 = slope2*TEMP_beta2 + cst2 + np.random.normal(0, stdEpsP2, TEMP_beta2.size)
PREC_beta2_double = slope2*TEMP_beta2_double + 500-slope2*20 + np.random.normal(0, stdEpsP2, TEMP_beta2_double.size)
TEQ_LTT2, BM_LTT2 = Teq(TEMP_beta2, LTT_beta2, PREC_beta2)
TEQ_HTT2, BM_HTT2 = Teq(TEMP_beta2, HTT_beta2, PREC_beta2)

data2 = np.concatenate((TEQ_LTT2, TEQ_HTT2), axis=None)
kde2 = gaussian_kde(data2)
x_grid2 = np.linspace(-0.1, 1.1, 1000)
density_estimation2 = kde2(x_grid2)
peaks2, _ = find_peaks(density_estimation2, height=0.5)
num_modes2 = len(peaks2)
print("Number of modes for case 2):", num_modes2)

corr2 = round(np.corrcoef(TEMP_beta2, PREC_beta2)[0,1],2)
print('slope2:', round(slope2,2))
print("corr coeff 2:", corr2)

# Graph ---------------------------------------------------------------------------------
# Precipitation(Temperature) scatter
fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(figsize=(8,6), nrows=2, ncols=2)    
PrecTemp_scatter(ax11, TEMP_beta1, PREC_beta1, TEQ_ltt, r'$\theta_L$~$U$(0.2, 0.3); P-T indep', '\n'.join(('P~$N$({}, {}\u00b2)'.format(meanP1, stdP1), 'T~$N$({}, {}\u00b2)'.format(meanT1, stdT1))))
PrecTemp_scatter(ax12, TEMP_beta2, PREC_beta2, TEQ_ltt, r'$\theta_L$~$U$(0.2, 0.3); P-T corr (r={})'.format(corr2), '\n'.join(('P=a$\cdot$T+b+$\u03B5_{P}$, $\u03B5_{P}$~$N$(0, 300\u00b2)', 'T~$N$({}, {}\u00b2)'.format(meanT2, stdT2))))
PrecTemp_scatter(ax21, TEMP_beta1, PREC_beta1, TEQ_htt, r'$\theta_H$~$U$(0.7, 0.8); P-T indep', '\n'.join(('P~$N$({}, {}\u00b2)'.format(meanP1, stdP1), 'T~$N$({}, {}\u00b2)'.format(meanT1, stdT1))))
PrecTemp_scatter(ax22, TEMP_beta2, PREC_beta2, TEQ_htt, r'$\theta_H$~$U$(0.7, 0.8); P-T corr (r={})'.format(corr2), '\n'.join(('P=a$\cdot$T+b+$\u03B5_{P}$, $\u03B5_{P}$~$N$(0, 300\u00b2)', 'T~$N$({}, {}\u00b2)'.format(meanT2, stdT2))))
plt.tight_layout()
plt.show()

# Precipitation(Theta) scatter
fig, ((ax11, ax12)) = plt.subplots(figsize=(7.5,3), nrows=1, ncols=2)
PrecTheta_scatter(ax11, TT_beta1, PREC_beta1_double, TEQ2, 'P,T~$N$; P-T indep', '\n'.join((r'$\theta_L$~$U$(0.2, 0.3)', r'$\theta_H$~$U$(0.7, 0.8)')))
PrecTheta_scatter(ax12, TT_beta2, PREC_beta2_double, TEQ2, 'P,T~$N$; P-T corr (r={})'.format(corr2), '\n'.join((r'$\theta_L$~$U$(0.2, 0.3)', r'$\theta_H$~$U$(0.7, 0.8)')))
plt.tight_layout()
plt.show()

# Tree fraction histogram
fig, ((ax1, ax2)) = plt.subplots(figsize=(6.1,3), nrows=1, ncols=2)
TreeCover_hist(ax1, TEQ_LTT1, TEQ_HTT1, r'P,T~$N$; $\theta$~$U$; P-T indep') 
TreeCover_hist(ax2, TEQ_LTT2, TEQ_HTT2, r'P,T~$N$; $\theta$~$U$; P-T corr (r={})'.format(corr2))
plt.tight_layout()
plt.show()
