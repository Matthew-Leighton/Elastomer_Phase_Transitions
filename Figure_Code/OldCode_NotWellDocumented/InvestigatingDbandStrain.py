import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

bellepsilondata=np.array([0,0.1,0.25,0.3])
belltwistdata = np.array([16,14,12,11])*np.pi/180


def epsilon_Dband(psi,psi0,mean_psi0_squared):
	return (1/2)*(1- ((psi)/psi0)**2) * mean_psi0_squared

def psi(epsilon,psi0,mean_psi0_squared):
	return np.sqrt(1-(2*epsilon/mean_psi0_squared)) * psi0

def psifinal(strain,psi0,zeta):
	return psi0*(zeta-1)/(zeta*strain**(3/2) - strain**(-3/2))



plt.scatter(bellepsilondata,belltwistdata,color='k',marker='D')

#epsilonlist =np.linspace(0,0.5,num=100)

'''
plt.plot(epsilonlist, psi(epsilonlist,16*np.pi/180,60*np.pi/180),lw=3)

plt.xlim(0,0.5)
plt.ylim(0,17*np.pi/180)
plt.xlabel(r'$\epsilon_{D-Band}$',fontsize=16)
plt.ylabel(r'$\langle \psi\rangle$',fontsize=16)

#plt.legend(loc='best',fancybox=True,fontsize='x-large')

plt.minorticks_on()
plt.tick_params(axis='x', labelsize=14)
plt.tick_params(axis='y', labelsize=14)

plt.tight_layout(pad=0.5)

plt.show()'''


zetalist = [1.05,1.1,1.15,1.2]
psi0=16*np.pi/180
mean_psi0_squared =0.01

lambdalist = np.linspace(1,1.05,num=100)



for zeta in zetalist:
	psilist = psifinal(lambdalist,psi0,zeta)
	epsilonlist = epsilon_Dband(psilist,psi0,mean_psi0_squared)

	plt.plot(epsilonlist*100,psilist,lw=3,label='$\zeta=$'+str(zeta))

plt.xlim(0,0.5)
plt.ylim(0,17*np.pi/180)
plt.xlabel('D-band strain (\%)',fontsize=16)
plt.ylabel(r'$\langle \psi\rangle$',fontsize=16)

plt.legend(loc='best',fontsize=16)

plt.minorticks_on()
plt.tick_params(axis='x', labelsize=14)
plt.tick_params(axis='y', labelsize=14)

plt.tight_layout(pad=0.5)

plt.show()


