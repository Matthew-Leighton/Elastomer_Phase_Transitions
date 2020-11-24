import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from matplotlib import rc
import matplotlib.gridspec as gridspec
from pyswarm import pso
from scipy.optimize import minimize


rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def Psi(psi0,strain,zeta):
	#strain is lambda
    x = ((zeta+1)*(strain**3 -1) + (zeta-1)*(strain**3 +1) * np.cos(2*psi0))/( 2*strain**(3/2)*(zeta-1)*np.sin(2*psi0) )
    return (1/2) * np.arctan(1/x)

def ExactPsi(psi0,strainlist,zeta):
	#Returns correct branch of the solution
	psi = Psi(psi0,strainlist,zeta)
	for i in range(len(strainlist)):
		if psi[i]<0:
			psi[i]=psi[i]+ (np.pi/2)
	return psi

def ExactPsi(psi0,strain,zeta):
    x = ((zeta+1)*(strain**3 -1) + (zeta-1)*(strain**3 +1) * np.cos(2*psi0))/( 2*strain**(3/2)*(zeta-1)*np.sin(2*psi0) )
    psi = (1/2) * np.arctan(1/x)
    if psi<0:
        psi=psi+ (np.pi/2)
    return psi

def f(strain,zeta,psi0):
    psi = ExactPsi(psi0,strain,zeta)
    fe = 1/strain + (1/strain)*(1+(zeta-1)*np.sin(psi0)**2)*(1+(1/zeta -1)*np.sin(psi)**2)
    fe += (strain**2)*(1+(zeta-1)*np.cos(psi0)**2)*(1+(1/zeta -1)*np.cos(psi)**2)
    fe += np.sqrt(strain)*(2 - zeta - 1/zeta)*np.sin(2*psi0)*np.sin(2*psi)/2
    return fe/2


N=1000
zeta=1.3
strainlist = np.linspace(0.8,1.1,num=N)
lslist = ['-','-.','--',':']

def functiontominimize(criticalstrains):
    strain1 = criticalstrains[0]
    strain2 = criticalstrains[1]
    f1 = f(strain1,zeta,psi0)
    f2 = f(strain2,zeta,psi0)
    fprime1 = (f1 - f(strain1+0.0001,zeta,psi0))/(-0.0001)
    fprime2 = (f2 - f(strain2+0.0001,zeta,psi0))/(-0.0001)
    return abs( (f2 - f1)/(strain2-strain1)  - fprime1) + abs( (f2 - f1)/(strain2-strain1)  - fprime2) + abs(fprime1-fprime2)

psi0list = [10**(-6),0.1,0.3,0.5] #[10**(-6),0.1,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4]
labellist= ['$\psi_0=0$','$\psi_0=0.1$','$\psi_0=0.3$','$\psi_0=0.5$'] #['0','0.1','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4']

fig, ax = plt.subplots()

for j in range(len(psi0list)):
    psi0 = psi0list[j]
    
    
    commont = minimize(functiontominimize,[0.85,0.99],bounds = ((0.8,0.9),(0.95,1)),tol=10**(-8)).x

    newf = np.zeros(N)
    for i in range(N):
        if strainlist[i]<commont[0]:
            newf[i] = f(strainlist[i],zeta,psi0)
        elif strainlist[i]>commont[1]:
            newf[i] = f(strainlist[i],zeta,psi0)
        else:
            x = (strainlist[i]-commont[0])/(commont[1]-commont[0])
            f1 = f(commont[0],zeta,psi0)
            f2 = f(commont[1],zeta,psi0)
            newf[i] = x*(f2-f1) + f1
            
        if psi0list[j]>=0.4: #psi0list[j]==0.4 or 
            newf[i] = f(strainlist[i],zeta,psi0)

    stress=np.zeros(N)        
    for i in range(N):
        if i>0:
            stress[i] = (newf[i] - newf[i-1])/(strainlist[i]-strainlist[i-1])
    stress[0] = stress[1] - (stress[2]-stress[1])
    
    ax.plot((strainlist-1)*100,stress,label=labellist[j],lw=3,ls=lslist[j])


ax.set_ylabel('$\\frac{\sigma}{\mu}$',fontsize=20,rotation=0,labelpad=15x)
ax.set_xlabel('Fibril Strain',fontsize=16)
ax.set_ylim(-0.6,0.3)
ax.set_xlim(-20,10)
ax.legend(loc='best',fontsize=16)

ax.set_xticks([-20,-15,-10,-5,0,5,10])
ax.set_xticklabels(['$-20\%$','$-15\%$', '$-10\%$','$-5\%$', '$0\%$','$5\%$','$10\%$'],fontsize=14)
ax.tick_params(axis='y', labelsize=14)
ax.minorticks_on()

plt.show()