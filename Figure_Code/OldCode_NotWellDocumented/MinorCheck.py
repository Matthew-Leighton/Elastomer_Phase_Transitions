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

def ExactPsi(psi0,strainlist,zeta):
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
	print([0.6,spinodalpoints[1]],[spinodalpoints[0],1])
	#common_t = minimize(functiontominimize,[np.mean([0,spinodalpoints[0]-0.001]),np.mean([spinodalpoints[1]+0.001,1])],bounds = ((0.7,spinodalpoints[0]),(spinodalpoints[1],1)),args=(psi0,zeta),tol=10**(-16)).x
	common_t = pso(functiontominimize,[0.7,spinodalpoints[1]-0.000001],[spinodalpoints[0],1],args=(psi0,zeta))[0]    
	return common_t


def PlotCoexistenceRegion():
	Coexist_Start = []
	Coexist_End = []

	for j in range(len(psi0list)):
		start_spinodal=np.zeros(N)
		end_spinodal=np.zeros(N)
		start_spinodal[:]=np.NaN
		end_spinodal[:]=np.NaN

		start_coexist=np.zeros(N)
		end_coexist=np.zeros(N)
		start_coexist[:]=np.NaN
		end_coexist[:]=np.NaN

		for i in range(N):
			psi0 = psi0list[j]
			zeta = zetalist[i]

			count=0
			for k in range(N):
				if count==0:
					if is_spinodal(lambdalist[k],psi0,zeta)==1:
						start_spinodal[i]=lambdalist[k]
						count+=1
				elif count==1:
					if is_spinodal(lambdalist[k],psi0,zeta)==0:
						end_spinodal[i] = lambdalist[k]
						count+=1

			if start_spinodal[i]>0 and end_spinodal[i]>0:
				commont = commontangentpoints(psi0,zeta,[start_spinodal[i],end_spinodal[i]])
				start_coexist[i] = commont[0]
				end_coexist[i] = commont[1]

		Coexist_Start.append(start_coexist)
		Coexist_End.append(end_coexist)
		ax1.plot(100*(start_coexist-1),zetalist,color=colorlist[j],label='$\psi_0 = $'+str(psi0list[j]),lw=3)
		ax1.plot(100*(end_coexist-1),zetalist,color=colorlist[j],lw=3)

	return Coexist_Start,Coexist_End

def Plot_Stress_Line(zetalist,Coexist_Start,Coexist_End,psi0list):
	N = len(zetalist)
	sigmalist=np.zeros(N)
	sigmalist[:] = np.NaN

	for j in range(len(psi0list)):
		psi0 = psi0list[j]
		start_coexist = Coexist_Start[j]
		end_coexist = Coexist_End[j]

		for i in range(N):
			zeta = zetalist[i]
			fstart = f(start_coexist[i],zeta,psi0)
			fend = f(end_coexist[i],zeta,psi0)
			sigmalist[i] = (fstart-fend)/(start_coexist[i]-end_coexist[i])

		ax2.plot(sigmalist,zetalist,color=colorlist[j],label='$\psi_0 = $'+str(psi0list[j]),lw=3)







N=50
lambdalist=np.linspace(0.8,1.0,num=N)
zetalist=np.linspace(1.01,2,num=N)
psi0list = [0.1]
colorlist = ['tab:orange']


fig=plt.figure()
gs=gridspec.GridSpec(2,1,width_ratios=[1],height_ratios=[1,1])
ax1=plt.subplot(gs[0])
ax2=plt.subplot(gs[1])
ax1.minorticks_on()
ax2.minorticks_on()

Coexist_Start,Coexist_End = PlotCoexistenceRegion()

ax1.fill_between(lambdalist,Coexist_Start[0],Coexist_End[0],color='orange')

#Plot_Stress_Line(zetalist,Coexist_Start,Coexist_End,psi0list)



ax1.set_title('A)',loc='left',fontsize=20)
ax1.set_ylabel('$\zeta$',fontsize=20,rotation=0,labelpad=15)
ax1.set_xlabel('Fibril Strain',fontsize=16)
ax1.legend(loc='best',fontsize=16)
ax1.set_ylim(1,1.95)
ax1.set_xlim(-20,0)

ax1.set_xticks([-30,-25,-20,-15,-10,-5,0])
ax1.set_xticklabels(['$-30\%$','$-25\%$','$-20\%$','$-15\%$', '$-10\%$','$-5\%$', '$0\%$'],fontsize=14)



ax2.set_title('B)',loc='left',fontsize=20)
ax2.set_ylabel('$\zeta$',fontsize=20,rotation=0,labelpad=15)
ax2.set_xlabel('Stress $\sigma/\mu$',fontsize=16)
#ax2.legend(loc='best',fontsize=16)
ax2.set_ylim(1,1.8)
ax2.set_xlim(-0.25,0)

ax2.set_xticks([-0.25,-0.2,-0.15,-0.1,-0.05,0])
#ax2.set_xticklabels(['$-20\%$','$-15\%$', '$-10\%$','$-5\%$', '$0\%$'],fontsize=14)

#ax1.tick_params(axis='x', labelsize=14)
ax1.tick_params(axis='y', labelsize=14)
ax2.tick_params(axis='x', labelsize=16)
ax2.tick_params(axis='y', labelsize=14)

plt.tight_layout(pad=0.5)
plt.show()