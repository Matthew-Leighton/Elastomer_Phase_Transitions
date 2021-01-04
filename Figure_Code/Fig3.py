import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from matplotlib import rc
import matplotlib.gridspec as gridspec
from pyswarm import pso
import scipy.optimize as opt
from scipy.integrate import quad
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)



def ExactPsi(psi0,strain,zeta):
    x = ((zeta+1)*(strain**3 -1) + (zeta-1)*(strain**3 +1) * np.cos(2*psi0))/( 2*strain**(3/2)*(zeta-1)*np.sin(2*psi0) )
    psi = (1/2) * np.arctan(1/x)
    if psi<0:
        psi=psi+ (np.pi/2)
    return psi

def integrand(r,strain,zeta,psi0_surf):
    psi0 = psi0_surf * r
    psi = ExactPsi(psi0*r,strain,zeta)
    fe = 1/strain + (1/strain)*(1+(zeta-1)*np.sin(psi0)**2)*(1+(1/zeta -1)*np.sin(psi)**2)
    fe += (strain**2)*(1+(zeta-1)*np.cos(psi0)**2)*(1+(1/zeta -1)*np.cos(psi)**2)
    fe += np.sqrt(strain)*(2 - zeta - 1/zeta)*np.sin(2*psi0)*np.sin(2*psi)/2
    return fe * r

def volume_averaged_f(strain,zeta,psi0_surf):
    fe = quad(integrand,0,1,args=(strain,zeta,psi0_surf))[0]
    
    return fe



# Figure Setup
fig=plt.figure()
gs=gridspec.GridSpec(2,1,width_ratios=[1],height_ratios=[1,1])
ax1=plt.subplot(gs[0])
ax2=plt.subplot(gs[1])
ax1.minorticks_on()
ax2.minorticks_on()



# A plot: (fe landscape)

N=100
strainlist = np.linspace(0.8,1.1,num=N)
zeta = 1.3


psi0_surf = 0.1
psi_surf = []
f = []
for strain in strainlist:
    f.append(volume_averaged_f(strain,zeta,psi0_surf))
ax1.plot(100*(strainlist-1),f,label='$\psi_0^{surf}=0.1$',linewidth=2,color='blue')

psi0_surf = 0.3
psi_surf = []
f = []
for strain in strainlist:
    f.append(volume_averaged_f(strain,zeta,psi0_surf))
ax1.plot(100*(strainlist-1),f,label='$\psi_0^{surf}=0.3$',linewidth=2.5,linestyle='-.',color='green')

psi0_surf = 0.5
psi_surf = []
f = []
for strain in strainlist:
    f.append(volume_averaged_f(strain,zeta,psi0_surf))
ax1.plot(100*(strainlist-1),f,label='$\psi_0^{surf}=0.5$',linewidth=3,linestyle='--',color='xkcd:orange')


ax1.set_xlabel('Fibril Strain',fontsize=16)
ax1.set_ylabel(r'$\frac{\langle f \rangle}{\mu}$',fontsize=30,rotation=0,labelpad=15)
ax1.legend(loc='best',fontsize=16)
ax1.set_xlim(-20,10)
ax1.set_xticks([-20,-15,-10,-5,0,5,10])
ax1.set_xticklabels(['$-20\%$','$-15\%$', '$-10\%$','$-5\%$', '$0\%$','$5\%$','$10\%$'],fontsize=14)


# B plot (psi examples):

N=1000
radius = np.linspace(0,1,num=N)
zeta = 1.3

psi0 = radius*0.3
ax2.plot(radius,psi0,label='$\psi_0$',lw=1.5,linestyle='-',color='black')


strain = 1.05
psi=[]
for j in range(N):
	psi.append(ExactPsi(psi0[j],strain,zeta))
ax2.plot(radius,psi,label='$5\%$ Strain',lw=2,ls=':',color='green')

strain = 0.98
psi=[]
for j in range(N):
	psi.append(ExactPsi(psi0[j],strain,zeta))
ax2.plot(radius,psi,label='$-2\%$ Strain',lw=2.5,ls='-.',color='blue')

strain = 0.8
psi=[]
for j in range(N):
	psi.append(ExactPsi(psi0[j],strain,zeta))
psi[0]=np.pi/2
ax2.plot(radius,psi,label='$-20\%$ Strain',lw=3,ls='--',color='xkcd:orange')


ax2.set_xlabel('$r/R$',fontsize=16)
ax2.set_ylabel('$\psi$',fontsize=30)
ax2.set_xlim(0,1)
ax2.set_ylim(0,1.6)
ax2.legend(loc='best',fontsize=16)


ax1.set_title('A)',loc='left',fontsize=20)
ax2.set_title('B)',loc='left',fontsize=20)

#ax1.tick_params(axis='x', labelsize=14)
ax1.tick_params(axis='y', labelsize=14)
ax2.tick_params(axis='x', labelsize=14)
ax2.tick_params(axis='y', labelsize=14)
plt.tight_layout(pad=0.5)

plt.show()
