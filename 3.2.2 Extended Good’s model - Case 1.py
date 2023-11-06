# =======================================================================================================================================================
# #### Case 1: two variables with unimodal input
# =======================================================================================================================================================

# =======================================================================================================================================================
# #### Figure 3.5 Outcome of Case 1
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

np.random.seed(123)

# Equilibrium --------------------------------------------------------------------------
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
# --------------------------------------------------------------------------------------
def productivityTemp(temp):               # productivity (beta) as a function of temperature                     
    #betaT = np.sin((temp/15)*np.pi-4.2)  # contribution of temperature to productivity
    betaT = -0.0175*(temp-20)*(temp-35)
    return betaT

def productivityPrec(p):       # productivity (beta) as a function of precipitation                 
    alpha = 0.000054           # original from Bathiany et al. (2013): 0.0011
    c = 138                    # original from Bathiany et al. (2013): 28
    gamma = 0.00017*9100

    P1 = c*np.exp(gamma/2) # =~300                                       # original: 60.68548303127268
    P2 = c*np.exp(gamma/2) + np.exp(gamma)/np.sqrt(0.03*alpha) # =~4000  # original: 878.3903705049846
    
    betaP = np.where(p < P1, 0.0, np.where(p <= P2, 1.0 - 1.0/(1+alpha*((p-P1)/np.exp(gamma))**2) , 1.0)) # "if" statement for matrix

    return betaP

mort = 0.15 # only for the graph Productivity scatter; should be same as below
def Teq(Temp, Prec):
    mort = 0.15
    prod = productivityTemp(Temp)*productivityPrec(Prec)
    BM = prod/mort
    
    T_eq = np.where(BM < 1, 0, 1-1/BM)  # "if" statement for matrix

    return T_eq, BM

def PrecTemp_scatter(ax, temp, prec, title, text, fontsize=10):
    CT=ax.contourf(TEMPv, PRECv, TEQ, 20)
    ax.scatter(temp, prec, s=1, color='black')
    ax.set_xlabel('Temperature [°C]', fontsize=fontsize)
    ax.set_ylabel('Precipitation [mm/yr]', fontsize=fontsize)
    ax.set_title(title, fontsize=12, loc='left', weight='bold')
    AT=AnchoredText(text, loc=1, prop=dict(size=8))
    AT.patch.set_alpha(0.5)
    ax.add_artist(AT)
    ax.set_xlim(20,35)
    ax.set_ylim(0,3000)
    plt.colorbar(CT, location='right', label='Tree Cover [%]')
     
def TreeCover_hist(ax, treecover, title, fontsize=10):
    ax.hist(treecover, bins=20, alpha=0.5)  # density=True
    ax.set_xlabel('Tree fraction estimate', fontsize=fontsize)
    ax.set_ylabel('Frequency', fontsize=fontsize)
    ax.set_title(title, fontsize=12, loc='left', weight='bold')
    ax.set_xlim(0,1)
    ax.set_ylim(0,500)
    
def Productivity_hist(ax, prod, title, fontsize=10):
    ax.hist(prod, bins=20, alpha=0.5, color='green')
    ax.set_xlabel('Productivity fraction', fontsize=fontsize)
    ax.set_ylabel('Productivity frequency', fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)
    ax.set_xlim(0,1)
    ax.set_ylim(0,500)

def TreeCoverEQ_scatter(ax, bm, teq, title, fontsize=10):
    ax.plot(bm_eq, x_eq, color='red') 
    ax.scatter(bm, teq, s=10, color='black')
    ax.set_xlabel('\u03B2/m', fontsize=fontsize)
    ax.set_ylabel('Tree Cover', fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)
    ax.set_xlim(0,5)
    ax.set_ylim(0,1)   

# Data -------------------------------------------------------------------------------    
# Make data for contour
TEMP = np.linspace(20, 35, 100)
PREC = np.linspace(0, 4000, 100)               
TEMPv, PRECv = np.meshgrid(TEMP, PREC)
TEQ, BM = Teq(TEMPv, PRECv)

# Make data for Teq vs productivity (normal)
datasize1 = 2000
meanT1, stdT1 = 27.5, 3
meanP1, stdP1 = 1000, 500  # peak the center where the gradient is sharper
TEMP_beta1 = np.random.normal(meanT1, stdT1, datasize1)   
PREC_beta1 = np.random.normal(meanP1, stdP1, datasize1)   
TEQ_beta1, BM_beta1 = Teq(TEMP_beta1, PREC_beta1)

# Make data for Teq vs productivity (normal-corr)
datasize2 = 2000
meanT2, stdT2 = 27.5, 3
slope2 = (1500-500)/(35-20) # =1000/15=200/3=~66.67
cst2 = 500-slope2*20        # from slope=(500-cst2)*20; cst2=500-(200/3)*20=500-4000/3=-833.33        200/3*20-833.33=500
stdEpsP2 = 200
TEMP_beta2 = np.random.normal(meanT2, stdT2, datasize2)
PREC_beta2 = slope2*TEMP_beta2 + cst2 + np.random.normal(0, stdEpsP2, TEMP_beta2.size)  
   # PREC as linear transformation of a normally distributed random variable TEMP~N(mu_X, sigma_X^2):
   # E[Y]=E[aX+b]=aE[X]+b=a*mu_X+b
   # Var[Y]=Var[aX+b]=a^2*Var[X]=a^2*sigma_X^2
   # And sum of two normal distrib (TEMP and error term) is normal (PREC)
TEQ_beta2, BM_beta2 = Teq(TEMP_beta2, PREC_beta2)
corr2 = round(np.corrcoef(TEMP_beta2, PREC_beta2)[0,1],2)
print("corr coeff 2:", corr2)
print("a:", round(slope2,2))
print("b:", round(cst2,2))

# Graph ------------------------------------------------------------------------------
# Precipitation(Temperature) scatter
fig, ((ax1, ax2)) = plt.subplots(figsize=(7.5,3), nrows=1, ncols=2)
PrecTemp_scatter(ax1, TEMP_beta1, PREC_beta1, 'a', '\n'.join(('P~$N$({}, {}\u00b2)'.format(meanP1, stdP1), 'T~$N$({}, {}\u00b2)'.format(meanT1, stdT1))))
PrecTemp_scatter(ax2, TEMP_beta2, PREC_beta2, 'b', '\n'.join(('P=a$\cdot$T+b+$\u03B5_P$, $\u03B5_P$~$N$(0, {}\u00b2)'.format(stdEpsP2), 'T~$N$({}, {}\u00b2)'.format(meanT2, stdT2))))
plt.tight_layout()
plt.show()

# Tree fraction histogram
fig, ((ax1, ax2)) = plt.subplots(figsize=(6,3), nrows=1, ncols=2)
TreeCover_hist(ax1, TEQ_beta1, 'a')
TreeCover_hist(ax2, TEQ_beta2, 'b')
plt.tight_layout()
plt.show()

# =======================================================================================================================================================
# #### Figure 3.6a-b Sensitivity analysis for Case 1
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

np.random.seed(123)

def productivityTemp(temp):                                   
    #betaT = np.sin((temp/15)*np.pi-4.2)
    betaT = -0.0175*(temp-20)*(temp-35)
    return betaT

def productivityPrec(p):                        
    alpha = 0.000054           
    c = 138                    
    gamma = 0.00017*9100

    P1 = c*np.exp(gamma/2) # =~300                                     
    P2 = c*np.exp(gamma/2) + np.exp(gamma)/np.sqrt(0.03*alpha) # =~4000
    
    betaP = np.where(p < P1, 0.0, np.where(p <= P2, 1.0 - 1.0/(1+alpha*((p-P1)/np.exp(gamma))**2) , 1.0))

    return betaP

def Teq(temp, prec):
    mort = 0.15
    prod = productivityTemp(temp)*productivityPrec(prec)
    BM = prod/mort
    
    T_eq = np.where(BM < 1, 0, 1-1/BM) 

    return T_eq, BM, prod

# -----------------------------------------------------------------------------------------
# 1) Make random data for Teq vs productivity (Precip: normal; Temp: normal; P-T indep)
datasize1 = 2000
meanT1_rng = [20, 25, 30, 35]        
meanP1 = 1000                      
stdT1 = 3
stdP1 = 500

param_values = ['20°C', '25°C', '30°C', '35°C']
line_colors = ['red', 'green', 'blue', 'aqua']

plt.figure(figsize=(4,3))
for i, meanT1 in enumerate(meanT1_rng):
    TEMP_beta1 = np.random.normal(meanT1, stdT1, datasize1)  # peak the center where the gradient is sharper
    PREC_beta1 = np.random.normal(meanP1, stdP1, datasize1)  # peak the center where the gradient is sharper
    TEQ_1, BM_1, PROD_1 = Teq(TEMP_beta1, PREC_beta1)
    sns.kdeplot(TEQ_1, label=param_values[i], color=line_colors[i])

plt.legend(loc='upper center')
plt.xlim(-0.02, 1.02)
plt.ylim(-0.2, 5.2)
plt.title('Sensitivity of $\mu_T$[°C] for $\mu_P$={}mm/yr'.format(meanP1))
plt.xlabel('Tree fraction estimate', fontsize=12)
plt.ylabel('Density', fontsize=12)
plt.show()

# -------------------------------------------------------------------------------------
# 2) Make random data for Teq vs productivity (Precip: normal; Temp: normal; P-T corr)
datasize2 = 2000
meanT2_rng = [20, 25, 30, 35]    
slope2 = (1500-500)/(35-20)

stdT2 = 3        # 3,5,7     # std of temperature 
stdEpsP2 = 200   # 300, 500  # std of noise in precipitation

param_values = ['20°C', '25°C', '30°C', '35°C']
line_colors = ['red', 'green', 'blue', 'aqua']

plt.figure(figsize=(4,3))
for i, meanT2 in enumerate(meanT2_rng):
    TEMP_beta2 = np.random.normal(meanT2, stdT2, datasize2)
    cst2 = 1000-meanT2*slope2  # from slope2 = (meanP2-cst2)/meanT2
    PREC_beta2 = slope2*TEMP_beta2 + cst2 + np.random.normal(0, stdEpsP2, TEMP_beta2.size)    
    TEQ_2, BM_2, PROD_2 = Teq(TEMP_beta2, PREC_beta2)
    sns.kdeplot(TEQ_2, label=param_values[i], color=line_colors[i])
    
plt.legend(loc='upper center')
plt.xlim(-0.02, 1.02)
plt.ylim(-0.2, 5.2)
plt.title('Sensitivity of $\mu_T$[°C] for slope $a$={}'.format(round(slope2,1)))
plt.xlabel('Tree fraction estimate', fontsize=12)
plt.ylabel('Density', fontsize=12)
plt.show()

# =======================================================================================================================================================
# #### Figure 3.6c-d Sensitivity analysis for Case 1
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

np.random.seed(123)

def productivityTemp(temp):                                   
    #betaT = np.sin((temp/15)*np.pi-4.2)
    betaT = -0.0175*(temp-20)*(temp-35)
    return betaT

def productivityPrec(p):                        
    alpha = 0.000054          
    c = 138                    
    gamma = 0.00017*9100

    P1 = c*np.exp(gamma/2) # =~300                                       
    P2 = c*np.exp(gamma/2) + np.exp(gamma)/np.sqrt(0.03*alpha) # =~4000
    
    betaP = np.where(p < P1, 0.0, np.where(p <= P2, 1.0 - 1.0/(1+alpha*((p-P1)/np.exp(gamma))**2) , 1.0))

    return betaP

def Teq(temp, prec):
    mort = 0.15
    prod = productivityTemp(temp)*productivityPrec(prec)
    BM = prod/mort
    
    T_eq = np.where(BM < 1, 0, 1-1/BM)

    return T_eq, BM, prod

# ------------------------------------------------------------------------------------------
# 1) Make random data for Teq vs productivity (Precip: normal; Temp: normal; FC: uniform; P-T indep)
datasize1 = 2000
meanP1_rng = [500, 1000, 1500, 2000]        
meanT1 = 27.5                     
stdT1 = 3   
stdP1 = 500

param_values = ['500mm/yr', '1000mm/yr', '1500mm/yr', '2000mm/yr']
line_colors = ['red', 'green', 'blue', 'aqua']

plt.figure(figsize=(4,3))
for i, meanP1 in enumerate(meanP1_rng):
    TEMP_beta1 = np.random.normal(meanT1, stdT1, datasize1)
    PREC_beta1 = np.random.normal(meanP1, stdP1, datasize1)
    TEQ_1, BM_1, PROD_1 = Teq(TEMP_beta1, PREC_beta1)
    sns.kdeplot(TEQ_1, label=param_values[i], color=line_colors[i])

plt.legend(loc='upper center', fontsize=9)
plt.xlim(-0.02, 1.02)
plt.ylim(-0.2, 5.2)
plt.title('Sensitivity of $\mu_P$[mm/yr] for $\mu_T$={}°C'.format(meanT1))
plt.xlabel('Tree fraction estimate', fontsize=12)
plt.ylabel('Density', fontsize=12)
plt.show()

# ----------------------------------------------------------------------------------------
# 2) Make random data for Teq vs productivity (Precip: normal; Temp: normal; FC: uniform; P-T corr)
datasize2 = 2000
meanP2_rng = [500, 1000, 1500, 2000]         
slope2 = (1500-500)/(35-20)  # =~66.67

meanT2 = 27.5
stdT2 = 3        # 3,5,7     # std of temperature   
stdEpsP2 = 200   # 300, 500  # std of noise in precipitation

param_values = ['500mm/yr', '1000mm/yr', '1500mm/yr', '2000mm/yr']
line_colors = ['red', 'green', 'blue', 'aqua']

plt.figure(figsize=(4,3))
for i, meanP2 in enumerate(meanP2_rng):
    cst2 = meanP2-27.5*slope2  # from slope2 = (meanP2-cst2)/(27.5-0)
    TEMP_beta2 = np.random.normal(meanT2, stdT2, datasize2)
    PREC_beta2 = slope2*TEMP_beta2 + cst2 + np.random.normal(0, stdEpsP2, TEMP_beta2.size)    
    TEQ_2, BM_2, PROD_2 = Teq(TEMP_beta2, PREC_beta2)
    sns.kdeplot(TEQ_2, label=param_values[i], color=line_colors[i])

plt.legend(loc='upper center', fontsize=9)
plt.xlim(-0.02, 1.02)
plt.ylim(-0.2, 5.2)
plt.title('Sensitivity of $\mu_P$[mm/yr] for slope $a$={}'.format(round(slope2,1)))
plt.xlabel('Tree fraction estimate', fontsize=12)
plt.ylabel('Density', fontsize=12)
plt.show()
