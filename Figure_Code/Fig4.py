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


# Post-strain twist angle:
def ExactPsi(psi0,strain,zeta):
    x = ((zeta+1)*(strain**3 -1) + (zeta-1)*(strain**3 +1) * np.cos(2*psi0))/( 2*strain**(3/2)*(zeta-1)*np.sin(2*psi0) )
    psi = (1/2) * np.arctan(1/x)
    if psi<0:
        psi=psi+ (np.pi/2)
    return psi

# Integrand for the next function:
def integrand(r,strain,zeta,psi0_surf):
    psi0 = psi0_surf * r
    psi = ExactPsi(psi0*r,strain,zeta)
    fe = 1/strain + (1/strain)*(1+(zeta-1)*np.sin(psi0)**2)*(1+(1/zeta -1)*np.sin(psi)**2)
    fe += (strain**2)*(1+(zeta-1)*np.cos(psi0)**2)*(1+(1/zeta -1)*np.cos(psi)**2)
    fe += np.sqrt(strain)*(2 - zeta - 1/zeta)*np.sin(2*psi0)*np.sin(2*psi)/2
    return fe * r

# Computes the volume-averaged free energy
def volume_averaged_f(strain,zeta,psi0_surf):
    fe = quad(integrand,0,1,args=(strain,zeta,psi0_surf))[0]
    
    return fe


### Functions needed for Coexist()

# First Derivative Calculation
def f1(strain,psi0_surf,zeta):
    h = 0.00001
    f_1 = (volume_averaged_f(strain+h,zeta,psi0_surf) - volume_averaged_f(strain-h,zeta,psi0_surf))/(h*2)
    return f_1

# Second derivative calculation
def f2(strain,psi0_surf,zeta):
    h = 0.00001
    f_2 = (volume_averaged_f(strain+h,zeta,psi0_surf) - 2*volume_averaged_f(strain,zeta,psi0_surf) + volume_averaged_f(strain-h,zeta,psi0_surf))/(h**2)
    return f_2

#Checks if point is in spinodal region
def is_spinodal(strain,psi0_surf,zeta):
    if f2(strain,psi0_surf,zeta)>=0:
        return 0
    elif f2(strain,psi0_surf,zeta)<0:
        return 1

# Function to minimize for coexistence calculation
def functiontominimize(criticalstrains,psi0_surf,zeta):
    strain1 = criticalstrains[0]
    strain2 = criticalstrains[1]
    f1 = volume_averaged_f(strain1,zeta,psi0_surf)
    f2 = volume_averaged_f(strain2,zeta,psi0_surf)
    fprime1 = (f1 - volume_averaged_f(strain1+0.00001,zeta,psi0_surf))/(-0.00001)
    fprime2 = (f2 - volume_averaged_f(strain2+0.00001,zeta,psi0_surf))/(-0.00001)
    return abs( (f2 - f1)/(strain2-strain1)  - fprime1) + abs( (f2 - f1)/(strain2-strain1)  - fprime2) + abs(fprime1-fprime2)

def Spinodal(psi0_surf,zeta):
    lambdalist = np.linspace(0.5,1,num=1000)
    count=0
    for k in range(1000):
        if count==0:
            if is_spinodal(lambdalist[k],psi0_surf,zeta)==1:
                start_spinodal=lambdalist[k]
                count+=1
        elif count==1:
            if is_spinodal(lambdalist[k],psi0_surf,zeta)==0:
                end_spinodal = lambdalist[k]
                count+=1

    if count==0:
        return np.array([np.NaN,np.NaN])
    else:
        return np.array([start_spinodal,end_spinodal])

#Calculates coexistence region for general psi_0 (still must be constant):
def Coexist(psi0_surf,zeta):
    spinodalpoints = Spinodal(psi0_surf,zeta)

    if np.isnan(spinodalpoints[0]):
        return np.NaN,np.NaN
    elif psi0_surf<0.4:
        common_t = pso(functiontominimize,[0.5,spinodalpoints[1]-0.000001],[spinodalpoints[0],1],args=(psi0_surf,zeta))[0]    
        return common_t[0],common_t[1]
    else:
        common_t = pso(functiontominimize,[0.91,spinodalpoints[1]],[spinodalpoints[0],1],args=(psi0_surf,zeta))[0]    
        return common_t[0],common_t[1]



# Figure Setup
fig=plt.figure()
gs=gridspec.GridSpec(2,1,width_ratios=[1],height_ratios=[1,1])
ax1=plt.subplot(gs[0])
ax2=plt.subplot(gs[1])
ax1.minorticks_on()
ax2.minorticks_on()



# A plot: (fe landscape)

N=1000
strainlist = np.linspace(0.8,1.1,num=N)
zeta = 1.3


psi0_surf = 0.1
psi_surf = []
f = []
for strain in strainlist:
    f.append(volume_averaged_f(strain,zeta,psi0_surf))
ax1.plot(100*(strainlist-1),f,label='$\psi_0^{surf}=0.1$',ls='-.',linewidth=2,color='blue')

psi0_surf = 0.3
psi_surf = []
f = []
for strain in strainlist:
    f.append(volume_averaged_f(strain,zeta,psi0_surf))
ax1.plot(100*(strainlist-1),f,label='$\psi_0^{surf}=0.3$',linewidth=2.5,linestyle='--',color='xkcd:orange')

psi0_surf = 0.5
psi_surf = []
f = []
for strain in strainlist:
    f.append(volume_averaged_f(strain,zeta,psi0_surf))
ax1.plot(100*(strainlist-1),f,label='$\psi_0^{surf}=0.5$',linewidth=3,linestyle=':',color='xkcd:red')


ax1.set_xlabel('Fibril Strain',fontsize=16)
ax1.set_ylabel(r'$\frac{\langle f \rangle}{\mu}$',fontsize=30,rotation=0,labelpad=15)
ax1.legend(loc='best',fontsize=14)
ax1.set_xlim(-20,10)
ax1.set_xticks([-20,-15,-10,-5,0,5,10])
ax1.set_xticklabels(['$-20\%$','$-15\%$', '$-10\%$','$-5\%$', '$0\%$','$5\%$','$10\%$'],fontsize=14)


# B plot (psi examples):

N=1000
radius = np.linspace(0,1,num=N)
zeta = 1.3
psi0 = radius*0.3

Lambda_High,Lambda_Low = Coexist(0.3,1.3)
print(Lambda_Low,Lambda_High)

strain = Lambda_High
psi=[]
for j in range(N):
	psi.append(ExactPsi(psi0[j],strain,zeta))
psi[0]=np.pi/2
ax2.plot(radius,psi,label='$\epsilon_H = -12.4\%$',lw=3,ls='--',color='xkcd:purple')

strain = Lambda_Low
psi=[]
for j in range(N):
    psi.append(ExactPsi(psi0[j],strain,zeta))
ax2.plot(radius,psi,label='$\epsilon_L = -2.6\%$',lw=2.5,ls='-.',color='xkcd:green')


ax2.plot(radius,psi0,label='$\psi_0(r)$',lw=1.5,linestyle='--',color='xkcd:orange')


strain = 1.05
psi=[]
for j in range(N):
    psi.append(ExactPsi(psi0[j],strain,zeta))
ax2.plot(radius,psi,label='$\epsilon = +5\%$',lw=2,ls=':',color='black')




ax2.set_xlabel('$r/R$',fontsize=18)
ax2.set_ylabel('$\psi$',fontsize=30)
ax2.set_xlim(0,1)
ax2.set_ylim(0,1.6)
ax2.legend(loc='best',fontsize=14)


ax1.set_title('A)',loc='left',fontsize=20)
ax2.set_title('B)',loc='left',fontsize=20)

#ax1.tick_params(axis='x', labelsize=14)
ax1.tick_params(axis='y', labelsize=14)
ax2.tick_params(axis='x', labelsize=14)
ax2.tick_params(axis='y', labelsize=14)
plt.tight_layout(pad=0.5)

plt.show()
