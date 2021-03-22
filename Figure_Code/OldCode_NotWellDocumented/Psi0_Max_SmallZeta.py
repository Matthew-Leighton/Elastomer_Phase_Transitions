import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from pyswarm import pso

#%matplotlib notebook

from matplotlib import rc

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)


def lambdac_exact(psi0,zeta):
    return ( (zeta-1) * (1 - np.sin(psi0)**2 + zeta*np.sin(psi0)**2)/ (zeta - 1 + np.cos(psi0)**2  - 2*zeta*np.cos(psi0)**2 + zeta**2 * np.cos(psi0)**2))**(1/3)


def f(strain,zeta,psi0):
    psi = ExactPsi(psi0,strain,zeta)
    fe = 1/strain + (1/strain)*(1+(zeta-1)*np.sin(psi0)**2)*(1+(1/zeta -1)*np.sin(psi)**2)
    fe += (strain**2)*(1+(zeta-1)*np.cos(psi0)**2)*(1+(1/zeta -1)*np.cos(psi)**2)
    fe += np.sqrt(strain)*(2 - zeta - 1/zeta)*np.sin(2*psi0)*np.sin(2*psi)/2
    return fe

def ExactPsi(psi0,strainlist,zetalist):
	psilist = np.zeros(len(zetalist))
	for i in range(len(zetalist)):
		strain = strainlist[i]
		zeta=zetalist[i]
		x = ((zeta+1)*(strain**3 -1) + (zeta-1)*(strain**3 +1) * np.cos(2*psi0))/( 2*strain**(3/2)*(zeta-1)*np.sin(2*psi0) )
		psi = (1/2) * np.arctan(1/x)
		if psi<0:
			psi=psi+ (np.pi/2)
		psilist[i]=psi
	return psilist



def f2(strain,psi0,zeta):
    h = 0.00001
    f_2 = (f(strain+h,zeta,psi0) - 2*f(strain,zeta,psi0) + f(strain-h,zeta,psi0))/(h**2)
    return f_2


#Getting psi_0_max to 4 sigfigs:
'''
zetalist=np.linspace(1,1.1,num=1000)



psi0=0.423
plt.plot(zetalist,f2(lambdac_exact(psi0,zetalist),psi0,zetalist),label='$\psi_0=$'+str(psi0))

psi0=0.4235
plt.plot(zetalist,f2(lambdac_exact(psi0,zetalist),psi0,zetalist),label='$\psi_0=$'+str(psi0))

psi0=0.4237
plt.plot(zetalist,f2(lambdac_exact(psi0,zetalist),psi0,zetalist),label='$\psi_0=$'+str(psi0))


psi0=0.4238
plt.plot(zetalist,f2(lambdac_exact(psi0,zetalist),psi0,zetalist),label='$\psi_0=$'+str(psi0))

psi0=0.4239
plt.plot(zetalist,f2(lambdac_exact(psi0,zetalist),psi0,zetalist),label='$\psi_0=$'+str(psi0))



plt.hlines(0,1,1.1,color='black',ls=':')

plt.xlabel('$\zeta$',fontsize=20)
plt.ylabel(r"$f''(\lambda_c)$",fontsize=20)
plt.ylim(-0.01,0.02)
plt.xlim(1,1.05)

plt.legend(loc='best',fontsize=14)
plt.tick_params(axis='x', labelsize=14)
	#ax2.tick_params(axis='x', labelsize=16)
plt.tick_params(axis='y', labelsize=14)

plt.tight_layout(pad=0.5)
plt.show()'''


zetalist=np.linspace(1.0015,1.0024,num=1000)

psi0=0.423848
plt.plot(zetalist,f2(lambdac_exact(psi0,zetalist),psi0,zetalist),label='$\psi_0=$'+str(psi0))

psi0=0.423849
plt.plot(zetalist,f2(lambdac_exact(psi0,zetalist),psi0,zetalist),label='$\psi_0=$'+str(psi0))

psi0=0.42385
plt.plot(zetalist,f2(lambdac_exact(psi0,zetalist),psi0,zetalist),label='$\psi_0=$'+str(psi0))


psi0=0.423851
plt.plot(zetalist,f2(lambdac_exact(psi0,zetalist),psi0,zetalist),label='$\psi_0=$'+str(psi0))

psi0=0.423852
plt.plot(zetalist,f2(lambdac_exact(psi0,zetalist),psi0,zetalist),label='$\psi_0=$'+str(psi0))



plt.hlines(0,1,1.1,color='black',ls=':')

plt.xlabel('$\zeta$',fontsize=20)
plt.ylabel(r"$f''(\lambda_c)$",fontsize=20)
plt.ylim(-0.00003,0.00007)
plt.xlim(1.0015,1.0024)

plt.legend(loc='best',fontsize=14)
plt.tick_params(axis='x', labelsize=14)
	#ax2.tick_params(axis='x', labelsize=16)
plt.tick_params(axis='y', labelsize=14)

plt.tight_layout(pad=0.5)
plt.show()

