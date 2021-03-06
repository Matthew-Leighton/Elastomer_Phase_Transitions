import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from matplotlib import rc
import matplotlib.gridspec as gridspec
from pyswarm import pso


rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def Psi(psi0,strain,zeta):
	#strain is lambda
    x = ((zeta+1)*(strain**3 -1) + (zeta-1)*(strain**3 +1) * np.cos(2*psi0))/( 2*strain**(3/2)*(zeta-1)*np.sin(2*psi0) )
    return (1/2) * np.arctan(1/x)

def ExactPsi_other(psi0,strainlist,zeta):
	#Returns correct branch of the solution
	psi = Psi(psi0,strainlist,zeta)
	for i in range(len(strainlist)):
		if psi[i]<0:
			psi[i]=psi[i]+ (np.pi/2)

	return psi

def ExactPsi(psi0,strain,zeta):
    x = ((zeta+1)*(strain**3 -1) + (zeta-1)*(strain**3 +1) * np.cos(2*psi0))/( 2*strain**(3/2)*(zeta-1)*np.sin(2*psi0) )
    psi = (1/2) * np.arctan(1/x)
    if psi<0:
        psi=psi+ (np.pi/2)
    return psi

def f(strain,zeta,psi0):
    psi = ExactPsi(psi0,strain,zeta)
    fe = 1/strain + (1/strain)*(1+(zeta-1)*np.sin(psi0)**2)*(1+(1/zeta -1)*np.sin(psi)**2)
    fe += (strain**2)*(1+(zeta-1)*np.cos(psi0)**2)*(1+(1/zeta -1)*np.cos(psi)**2)
    fe += np.sqrt(strain)*(2 - zeta - 1/zeta)*np.sin(2*psi0)*np.sin(2*psi)/2
    return fe/2

def f_other(strain,zeta,psi0):
    psi = ExactPsi_other(psi0,strain,zeta)
    fe = 1/strain + (1/strain)*(1+(zeta-1)*np.sin(psi0)**2)*(1+(1/zeta -1)*np.sin(psi)**2)
    fe += (strain**2)*(1+(zeta-1)*np.cos(psi0)**2)*(1+(1/zeta -1)*np.cos(psi)**2)
    fe += np.sqrt(strain)*(2 - zeta - 1/zeta)*np.sin(2*psi0)*np.sin(2*psi)/2
    return fe/2

def f1(strain,psi0,zeta):
    h = 0.00001
    f_1 = (f(strain+h,zeta,psi0) - f(strain-h,zeta,psi0))/(h*2)
    return f_1

def f2(strain,psi0,zeta):
    h = 0.00001
    f_2 = (f(strain+h,zeta,psi0) - 2*f(strain,zeta,psi0) + f(strain-h,zeta,psi0))/(h**2)
    return f_2

def is_spinodal(strain,psi0,zeta):
    if f2(strain,psi0,zeta)>=0:
        return 0
    elif f2(strain,psi0,zeta)<0:
        return 1


def functiontominimize(criticalstrains,psi0,zeta):
    strain1 = criticalstrains[0]
    strain2 = criticalstrains[1]
    f1 = f(strain1,zeta,psi0)
    f2 = f(strain2,zeta,psi0)
    fprime1 = (f1 - f(strain1+0.00001,zeta,psi0))/(-0.00001)
    fprime2 = (f2 - f(strain2+0.00001,zeta,psi0))/(-0.00001)
    return abs( (f2 - f1)/(strain2-strain1)  - fprime1) + abs( (f2 - f1)/(strain2-strain1)  - fprime2) + abs(fprime1-fprime2)


def commontangentpoints(psi0,zeta,spinodalpoints):
	print([0.7,spinodalpoints[1]],[spinodalpoints[0],1])
	#common_t = minimize(functiontominimize,[np.mean([0,spinodalpoints[0]-0.001]),np.mean([spinodalpoints[1]+0.001,1])],bounds = ((0.7,spinodalpoints[0]),(spinodalpoints[1],1)),args=(psi0,zeta),tol=10**(-16)).x
	common_t = pso(functiontominimize,[0.7,spinodalpoints[1]-0.000001],[spinodalpoints[0],1],args=(psi0,zeta))[0]    
	return common_t


zeta = 1.3
N=1000
psi0=0.1
lambdalist = np.linspace(0.8,1.1,num=N)
epsilonlist = np.linspace(-0.2,0.1,num=N)


count=0
for k in range(N):
	if count==0:
		if is_spinodal(lambdalist[k],psi0,zeta)==1:
			start_spinodal=lambdalist[k]
			count+=1
	elif count==1:
		if is_spinodal(lambdalist[k],psi0,zeta)==0:
			end_spinodal = lambdalist[k]
			count+=1

if start_spinodal>0 and end_spinodal>0:
	commont = commontangentpoints(psi0,zeta,[start_spinodal,end_spinodal])
	start_coexist = commont[0]
	end_coexist = commont[1]



fig, ax = plt.subplots()

plt.minorticks_on()


plt.plot(100*epsilonlist, f_other(lambdalist,zeta,0.1),label = '$\psi_0 = 0.1$',lw=3,ls='-',color = 'tab:orange')

plt.ylabel('$\\frac{f}{\mu}$',fontsize=30,rotation=0,labelpad=15)
plt.xlabel('Fibril Strain',fontsize=16)
plt.xlim(-18,8)
'''plt.scatter(100*(commont[0]-1),f(commont[0],zeta,psi0),marker='o',color='black',s=300)
plt.scatter(100*(commont[1]-1),f(commont[1],zeta,psi0),marker='o',color='black',s=300)

linlist=np.linspace(commont[0],commont[1],num=1000)
plt.plot(100*(linlist-1),f(commont[0],zeta,psi0) -(linlist-commont[0])/(commont[1]-commont[0])*(f(commont[0],zeta,psi0)-f(commont[1],zeta,psi0)),ls='--',color='red',lw=3)

plt.vlines(100*(commont[0]-1),0,2,ls=':')
plt.vlines(100*(commont[1]-1),0,2,ls=':')'''

plt.ylim(1.498,1.520)

ax.set_xticks([-15,-10,-5,0,5])
ax.set_xticklabels(['$-15\%$', '$-10\%$','$-5\%$', '$0\%$','$5\%$'],fontsize=14)
plt.tick_params(axis='y', labelsize=14)
plt.legend(loc='best',fontsize=16)

plt.tight_layout(pad=0.5)
plt.show()


