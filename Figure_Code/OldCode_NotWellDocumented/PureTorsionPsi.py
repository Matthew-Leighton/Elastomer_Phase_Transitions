import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from matplotlib import rc
import matplotlib.gridspec as gridspec
from pyswarm import pso
from scipy.optimize import minimize

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)


def Psi(psi0,tau,zeta,r):
	#strain is lambda
    x = (1 + zeta*(r**2 * tau**2 - 1))/(2*zeta*tau*r)
    return (1/2) * np.arctan(1/x)

def ExactPsi(psi0,strain,zeta,r):
    #returns correct branch of the solution
    x = (1 + zeta*(r**2 * tau**2 - 1))/(2*zeta*tau*r)
    psi = (1/2) * np.arctan(1/x)
    if psi<0:
        psi=psi+ (np.pi/2)
    return psi

#Params and setup
N=1000
zeta=1.3
psi0=0
rlist=np.linspace(0,1,num=N)
taulist = [-0.01,-0.1,-1,-10]
lslist = ['-','-.','--',':']
colorlist = ['black','blue','xkcd:orange','xkcd:red']
lwlist = [2,2.5,3,2]


labellist= [r'$\tau/L=-0.01$',r'$\tau/L=-0.1$',r'$\tau/L=-1$',r'$\tau/L=-10$'] #['0','0.1','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4']



# Compute and plot the stress-strain curves:

fig, ax = plt.subplots()

for j in range(len(taulist)):
	tau = taulist[j]

	psilist=np.zeros(N)
	for i in range(len(rlist)):
		psilist[i] = ExactPsi(psi0,tau,zeta,rlist[i])

	ax.plot(rlist,psilist,label=labellist[j],ls=lslist[j],color=colorlist[j],lw=lwlist[j])


ax.set_ylabel('$\psi$',fontsize=20,rotation=0,labelpad=15)
ax.set_xlabel('$r/R$',fontsize=20)
ax.set_xlim(0,1)
ax.set_ylim(0,1.7)
ax.legend(loc='best',fontsize=16)

ax.tick_params(axis='x', labelsize=14)
ax.tick_params(axis='y', labelsize=14)
ax.minorticks_on()

plt.show()