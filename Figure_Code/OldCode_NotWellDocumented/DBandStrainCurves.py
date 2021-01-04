import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc


rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


def periodstrain(strain,zeta,psi0_R):
	a = psi0_R*(zeta-1)/(zeta*strain**(3/2) - strain**(-3/2))

	period = 1 - a**2/4 + (19/288)*a**4

	return (period/ (1-psi0_R**2 /4 + (19/288)*psi0_R**4) - 1)*100



zetalist=[1.2,1.3,1.4,1.5]
strainlist = np.linspace(1,1.1,num=100)

fig=plt.figure()
gs=gridspec.GridSpec(1,3,width_ratios=[1,1,1],height_ratios=[1])

ax1=plt.subplot(gs[0])
ax2=plt.subplot(gs[1])
ax3=plt.subplot(gs[2])

for zeta in zetalist:
	ax1.plot(100*(strainlist-1),periodstrain(strainlist,zeta,0.2),label='$\zeta=$'+str(zeta),lw=3)
ax1.set_xlim(0,10)
ax1.set_ylim(0,)
ax1.set_xlabel('Fibril Strain (\%)',fontsize=16)
ax1.set_ylabel('D-Band Strain (\%)',fontsize=16)
ax1.legend(loc='best',fontsize=16)
ax1.minorticks_on()
ax1.tick_params(axis='x', labelsize=16)
ax1.tick_params(axis='y', labelsize=16)
ax1.set_title('A) $\psi_R=0.2$',fontsize=16)

for zeta in zetalist:
	ax2.plot(100*(strainlist-1),periodstrain(strainlist,zeta,0.3),label='$\zeta=$'+str(zeta),lw=3)
ax2.set_xlim(0,10)
ax2.set_ylim(0,)
ax2.set_xlabel('Fibril Strain (\%)',fontsize=16)
ax2.set_ylabel('D-Band Strain (\%)',fontsize=16)
ax2.legend(loc='best',fontsize=16)
ax2.minorticks_on()
ax2.tick_params(axis='x', labelsize=16)
ax2.tick_params(axis='y', labelsize=16)
ax2.set_title('B) $\psi_R=0.3$',fontsize=16)

for zeta in zetalist:
	ax3.plot(100*(strainlist-1),periodstrain(strainlist,zeta,0.4),label='$\zeta=$'+str(zeta),lw=3)
ax3.set_xlim(0,10)
ax3.set_ylim(0,)
ax3.set_xlabel('Fibril Strain (\%)',fontsize=16)
ax3.set_ylabel('D-Band Strain (\%)',fontsize=16)
ax3.legend(loc='best',fontsize=16)
ax3.minorticks_on()
ax3.tick_params(axis='x', labelsize=16)
ax3.tick_params(axis='y', labelsize=16)
ax3.set_title('C) $\psi_R=0.4$',fontsize=16)


plt.tight_layout(pad=0.5)
plt.show()