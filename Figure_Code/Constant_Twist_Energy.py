import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from matplotlib import rc
import matplotlib.gridspec as gridspec


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

def f(strain,zeta,psi0):
    psi = ExactPsi(psi0,strain,zeta)
    fe = 1/strain + (1/strain)*(1+(zeta-1)*np.sin(psi0)**2)*(1+(1/zeta -1)*np.sin(psi)**2)
    fe += (strain**2)*(1+(zeta-1)*np.cos(psi0)**2)*(1+(1/zeta -1)*np.cos(psi)**2)
    fe += np.sqrt(strain)*(2 - zeta - 1/zeta)*np.sin(2*psi0)*np.sin(2*psi)/2
    return fe/2



zeta = 1.3
N=1000
lambdalist = np.linspace(0.8,1.1,num=N)
epsilonlist = np.linspace(-0.2,0.1,num=N)

fig=plt.figure()
gs=gridspec.GridSpec(2,1,width_ratios=[1],height_ratios=[1,1])
ax1=plt.subplot(gs[0])
ax2=plt.subplot(gs[1])
ax1.minorticks_on()
ax2.minorticks_on()

ax1.plot(100*epsilonlist, ExactPsi(0.000001,lambdalist,zeta),label = '$\psi_0 = 0$',lw=3,ls='-')
ax1.plot(100*epsilonlist, ExactPsi(0.1,lambdalist,zeta),label = '$\psi_0 = 0.1$',lw=3,ls='-.')
ax1.plot(100*epsilonlist, ExactPsi(0.3,lambdalist,zeta),label = '$\psi_0 = 0.3$',lw=3,ls='--')
ax1.plot(100*epsilonlist, ExactPsi(0.5,lambdalist,zeta),label = '$\psi_0 = 0.5$',lw=3,ls=':')

ax2.plot(100*epsilonlist, f(lambdalist,zeta,0.000001),label = '$\psi_0 = 0$',lw=3,ls='-')
ax2.plot(100*epsilonlist, f(lambdalist,zeta,0.1),label = '$\psi_0 = 0.1$',lw=3,ls='-.')
ax2.plot(100*epsilonlist, f(lambdalist,zeta,0.3),label = '$\psi_0 = 0.3$',lw=3,ls='--')
ax2.plot(100*epsilonlist, f(lambdalist,zeta,0.5),label = '$\psi_0 = 0.5$',lw=3,ls=':')


ax1.set_title('A)',loc='left',fontsize=20)
ax1.set_ylabel('$\psi$',fontsize=20,rotation=0,labelpad=15)
#ax1.set_xlabel('Fibril Strain',fontsize=16)
ax1.legend(loc='best',fontsize=16)
ax1.set_ylim(-0.1,1.7)
ax1.set_xlim(-20,10)

ax1.set_xticks([-20,-15,-10,-5,0,5,10])
ax1.set_xticklabels(['$-20\%$','$-15\%$', '$-10\%$','$-5\%$', '$0\%$','$5\%$','$10\%$'],fontsize=14)



ax2.set_title('B)',loc='left',fontsize=20)
ax2.set_ylabel('$\\frac{f}{\mu}$',fontsize=30,rotation=0,labelpad=15)
ax2.set_xlabel('Fibril Strain',fontsize=16)
#ax2.legend(loc='best',fontsize=16)
#ax1.set_ylim(-0.1,1.7)
ax2.set_xlim(-20,10)

ax2.set_xticks([-20,-15,-10,-5,0,5,10])
ax2.set_xticklabels(['$-20\%$','$-15\%$', '$-10\%$','$-5\%$', '$0\%$','$5\%$','$10\%$'],fontsize=14)

#ax1.tick_params(axis='x', labelsize=14)
ax1.tick_params(axis='y', labelsize=14)
#ax2.tick_params(axis='x', labelsize=16)
ax2.tick_params(axis='y', labelsize=14)

plt.tight_layout(pad=0.5)
plt.show()


