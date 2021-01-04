import math as m
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar,minimize
import scipy.integrate as integrate
import matplotlib.gridspec as gridspec
from scipy.misc import derivative
from scipy.optimize import fsolve
import matplotlib.cm as cmap
from matplotlib.widgets import Slider, Button, RadioButtons
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)



def F(psi,gr,K,Lambda,delta,eta):
    return np.abs(gr - (1/2)*np.sin(2*psi) - K*np.tan(2*psi)*np.sin(psi)**2 - Lambda*((4*np.pi**2) - (eta**2) * np.cos(psi)**2)*np.tan(2*psi)*(delta*eta*gr)**2)

# This function returns psi(r)
#    inputs are the constants K and Lambda, a trial value of psi, an r value (gr), and the quantities delta(r) and eta(r)
def psifunction(gr,K,deltalist,etalist,Lambda):
    psilist=np.zeros(len(gr))
    for i in range(len(gr)):
        psilist[i] = minimize(F,psilist[i-1],args=(gr[i],K,Lambda,deltalist[i],etalist[i])).x
    return psilist



def CalculateStructure_2(gr,K,Lambda,omega):
	N = len(gr)

	psi=np.zeros(N)
	etalist=np.zeros(N)
	deltalist=np.ones(N)

	for i in range(N):
		if i==0:
			psi[i] = gr[i]
			quadint2 = (gr[i])*(gr[i]*np.cos(psi[i])**2)/2
			quadint4 = (gr[i])*(gr[i]*np.cos(psi[i])**4)/2
		else:
			psi[i] = minimize(F,psi[i-1],args=(gr[i],K,Lambda,deltalist[i-1],etalist[i-1])).x
			quadint2 += (gr[i]-gr[i-1]) * (gr[i]*np.cos(psi[i])**2 + gr[i-1]*np.cos(psi[i-1])**2)/2
			quadint4 += (gr[i]-gr[i-1]) * (gr[i]*np.cos(psi[i])**4 + gr[i-1]*np.cos(psi[i-1])**4)/2
			etalist[i] = np.sqrt(4*(np.pi**2)*quadint2/quadint4)
			deltalist[i] =  ( 1 - (8*np.pi**4 * Lambda/omega) + (4*np.pi**2 * (Lambda/omega)*(1/(gr[i]**2))*etalist[i]**2 * quadint2) )
			deltalist[i] = np.sqrt(max(deltalist[i],0))

	return psi,etalist,deltalist



K=100
Lambda=6
omega=0.1
gr = np.logspace(-2,2,num=1000)


psi,etalist,deltalist = CalculateStructure_2(gr,K,Lambda,omega)




##### Plot Molecular Strain Figure (Figure 9):
def PlotMolecularStrain(gr,psi,etalist,deltalist):

	r0loc = np.where(deltalist == min(deltalist))[0][0]
	vsloc = 150

	molecularstrain = ((2*np.pi/etalist[-1] - np.cos(psi))/np.cos(psi))
	molecularstrainsmall = ((2*np.pi/etalist[r0loc] - np.cos(psi[:r0loc]))/np.cos(psi[:r0loc]))
	molecularstrainverysmall = ((2*np.pi/etalist[vsloc] - np.cos(psi[:vsloc]))/np.cos(psi[:vsloc]))

	plt.plot(gr,molecularstrain*100,label='$R >R_0$',lw=3,color='blue')
	plt.plot(gr[:r0loc],molecularstrainsmall*100,label='$R = R_0$',lw=3,color='xkcd:orange',ls='-.')
	plt.plot(gr[:vsloc],molecularstrainverysmall*100,label='$R < R_0$',lw=3,color='xkcd:red',ls='--')

	plt.title('$K=100,\Lambda=6, \omega=0.1$',fontsize=16)
	plt.xlabel('$r$',fontsize=26)
	plt.ylabel('Molecular Strain (\%)',fontsize=18)
	plt.xscale('log')
	plt.xlim(0.01,130)
	plt.ylim(-0.5,0.2)
	plt.scatter(100,molecularstrain[-1]*100,marker='o',s=200,color='blue')
	plt.scatter(gr[r0loc],molecularstrainsmall[-1]*100,marker='o',s=200,color='xkcd:orange')
	plt.scatter(gr[vsloc],molecularstrainverysmall[-1]*100,marker='o',s=200,color='xkcd:red')
	plt.scatter(gr[np.where(deltalist == min(deltalist))], -0.47,marker='v',color='red',s=200)

	plt.minorticks_on()
	plt.tick_params(axis='x', labelsize=17)
	plt.tick_params(axis='y', labelsize=17)
	plt.tight_layout(pad=0.5)
	plt.legend(loc='best',fontsize=20)

	plt.show()



PlotMolecularStrain(gr,psi,etalist,deltalist)

