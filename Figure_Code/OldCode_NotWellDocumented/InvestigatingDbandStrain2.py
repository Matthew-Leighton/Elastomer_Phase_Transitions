import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

bellepsilondata=np.array([0,0.1,0.24,0.3])
belltwistdata = np.array([16,14,12,11])*np.pi/180


def psi(epsilon_D,psi0,psi0squared):
	return psi0 * np.sqrt(1 - 2*epsilon_D/psi0squared)

'''def AmplitudeChange(epsilon_D,L_over_W,mean_psi0_squared,psi0):
	psi_notnot = psi(epsilon_D,psi0,mean_psi0_squared)
	return -8*(np.pi**4)*L_over_W * (psi_notnot/psi0)**4 * epsilon_D'''

def AmplitudeChange(epsilon_D,L_over_W,mean_psi0_squared):
	return  -8*(np.pi**4)*L_over_W * epsilon_D * (1- 2*epsilon_D/mean_psi0_squared)**2


psierror = [np.pi/180,np.pi/180,np.pi/180,np.pi/180]
epsilonerror = [0,0.065,0.035,0.055]



epsilon_D_list = np.linspace(0,0.01,num=1000)
psi0=16*np.pi/180
mean_psi0_squaredlist = [0.01,0.02,0.04,0.08]#,0.12,0.16]



for i in mean_psi0_squaredlist:
	plt.plot(100*epsilon_D_list,psi(epsilon_D_list,psi0,i),label=r'$\langle\psi_0^2\rangle =$'+str(i),lw=3)

plt.errorbar(bellepsilondata,belltwistdata,yerr = psierror,xerr = epsilonerror,color='k')


plt.xlim(0,1)
plt.ylim(0,17*np.pi/180)
plt.xlabel('D-band strain (\%)',fontsize=16)
plt.ylabel(r'$\langle \psi\rangle$',fontsize=16)

plt.legend(loc='best',fontsize=16)

plt.minorticks_on()
plt.tick_params(axis='x', labelsize=14)
plt.tick_params(axis='y', labelsize=14)

plt.tight_layout(pad=0.5)

plt.show()
'''

L_over_W_list = [0.1,1,10,100]
mean_psi0_squared = 0.01#(5*np.pi/180)**2
psi0=16*np.pi/180

#plt.scatter(bellepsilondata,-8*(np.pi**4)*1*(belltwistdata/psi0)**4 * bellepsilondata)

for i in L_over_W_list:
	#plt.plot(100*epsilon_D_list,AmplitudeChange(epsilon_D_list,i,mean_psi0_squared,psi0),label=r'$\Lambda/\omega =$'+str(i),lw=3)
	plt.plot(100*epsilon_D_list,AmplitudeChange(epsilon_D_list,i,mean_psi0_squared),label=r'$\Lambda/\omega =$'+str(i),lw=3)


plt.xlim(0,0.8)
plt.ylim(-50,5)
plt.xlabel('D-band strain (\%)',fontsize=16)
plt.ylabel('Change in D-Band Amplitude(\%)',fontsize=16)
#plt.yscale('log')

plt.legend(loc='best',fontsize=16)

plt.minorticks_on()
plt.tick_params(axis='x', labelsize=14)
plt.tick_params(axis='y', labelsize=14)

plt.tight_layout(pad=0.5)

plt.show()'''



