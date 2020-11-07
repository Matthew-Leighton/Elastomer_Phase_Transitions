import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from matplotlib import rc

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

bellepsilondata=np.array([1-1,1.014-1,1.028-1,1.05-1])
bellfdata = np.array([1,14/16,12/16,11/16])
bellferror = np.array([np.sqrt(2)/16,np.sqrt((14/16)**2 + 1)/16,np.sqrt((12/16)**2 + 1)/16,np.sqrt((11/16)**2 + 1)/16])


def Psi_over_Psi0(strain,zeta):
	#strain is lambda
    return (zeta-1)/ (zeta*strain**(3/2) - strain**(-3/2))


fig, ax = plt.subplots()
ax.errorbar(100*bellepsilondata,bellfdata,bellferror,label='Corneal Fibrils (Bell et al, 2018)',fmt='D:k',ls='none',ms=10)


epsilonlist = np.linspace(0,0.1,num=1000)
lambdalist = epsilonlist + np.ones(1000)
ax.plot(100*(epsilonlist),Psi_over_Psi0(lambdalist,1.2),label='$\zeta=1.2$',linestyle='-',lw=3)
ax.plot(100*(epsilonlist),Psi_over_Psi0(lambdalist,1.3),label='$\zeta=1.3$',linestyle='--',lw=3)
ax.plot(100*(epsilonlist),Psi_over_Psi0(lambdalist,1.4),label='$\zeta=1.4$',linestyle='-.',lw=3)



ax.set_xlabel('Fibril Strain',fontsize=20,labelpad=10)
ax.set_ylabel(r'$\frac{\langle \psi\rangle}{\langle\psi_0\rangle}$',fontsize=30,rotation=0,labelpad=30)
ax.set_xlim(0,10)
ax.set_xticks([0,2,4,6,8,10])
ax.set_xticklabels(['0\%','2\%', '4\%','6\%', '8\%','10\%'],fontsize=16)
ax.set_yticks([0.4,0.5,0.6,0.7,0.8,0.9,1.0])
ax.set_yticklabels(['0.4','0.5','0.6','0.7','0.8','0.9','1.0'],fontsize=16)
ax.set_ylim(0.5,1.05)
ax.legend(loc='best',fancybox=True,fontsize='x-large')
ax.xaxis.set_minor_locator(MultipleLocator(0.2))
ax.yaxis.set_minor_locator(MultipleLocator(0.02))

#plt.subplots_adjust(left=0.16,bottom=0.18,top=0.97,right=0.96)
plt.tight_layout(pad=0.5)
plt.show()
