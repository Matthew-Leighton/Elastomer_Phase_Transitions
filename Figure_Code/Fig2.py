import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from matplotlib import rc
import matplotlib.gridspec as gridspec
from pyswarm import pso
import scipy.optimize as opt
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

# This function is needed for Coexist_0()
def coexist(x,zeta): 
    return [-2/x[0]**2+2*x[0]+(1+1/zeta)/x[1]**2-2*zeta*x[1], 4/x[0]-x[0]**2-2*(1+1/zeta)/x[1]+zeta*x[1]**2]

#Computes lambda_L and lambda_H for psi_0=0 and a given value of zeta
def Coexist_0(zeta):
	guess=[.5,.9]
	sol=opt.root(coexist,guess,args=(zeta))
	lambda_low = sol.x[0]
	lambda_high = sol.x[1]

	return lambda_low,lambda_high


### Functions needed for Coexist()

# First Derivative Calculation
def f1(strain,psi0,zeta):
    h = 0.00001
    f_1 = (f_coexist(strain+h,zeta,psi0) - f_coexist(strain-h,zeta,psi0))/(h*2)
    return f_1

# Second derivative calculation
def f2(strain,psi0,zeta):
    h = 0.00001
    f_2 = (f_coexist(strain+h,zeta,psi0) - 2*f_coexist(strain,zeta,psi0) + f_coexist(strain-h,zeta,psi0))/(h**2)
    return f_2

def ExactPsi2(psi0,strain,zeta):
    x = ((zeta+1)*(strain**3 -1) + (zeta-1)*(strain**3 +1) * np.cos(2*psi0))/( 2*strain**(3/2)*(zeta-1)*np.sin(2*psi0) )
    psi = (1/2) * np.arctan(1/x)
    if psi<0:
        psi=psi+ (np.pi/2)
    return psi

def f_coexist(strain,zeta,psi0):
    psi = ExactPsi2(psi0,strain,zeta)
    fe = 1/strain + (1/strain)*(1+(zeta-1)*np.sin(psi0)**2)*(1+(1/zeta -1)*np.sin(psi)**2)
    fe += (strain**2)*(1+(zeta-1)*np.cos(psi0)**2)*(1+(1/zeta -1)*np.cos(psi)**2)
    fe += np.sqrt(strain)*(2 - zeta - 1/zeta)*np.sin(2*psi0)*np.sin(2*psi)/2
    return fe/2

#Checks if point is in spinodal region
def is_spinodal(strain,psi0,zeta):
    if f2(strain,psi0,zeta)>=0:
        return 0
    elif f2(strain,psi0,zeta)<0:
        return 1

# Function to minimize for coexistence calculation
def functiontominimize(criticalstrains,psi0,zeta):
    strain1 = criticalstrains[0]
    strain2 = criticalstrains[1]
    f1 = f_coexist(strain1,zeta,psi0)
    f2 = f_coexist(strain2,zeta,psi0)
    fprime1 = (f1 - f_coexist(strain1+0.00001,zeta,psi0))/(-0.00001)
    fprime2 = (f2 - f_coexist(strain2+0.00001,zeta,psi0))/(-0.00001)
    return abs( (f2 - f1)/(strain2-strain1)  - fprime1) + abs( (f2 - f1)/(strain2-strain1)  - fprime2) + abs(fprime1-fprime2)

def Spinodal(psi0,zeta):
	lambdalist = np.linspace(0.1,1,num=1000)
	count=0
	for k in range(1000):
		if count==0:
			if is_spinodal(lambdalist[k],psi0,zeta)==1:
				start_spinodal=lambdalist[k]
				count+=1
		elif count==1:
			if is_spinodal(lambdalist[k],psi0,zeta)==0:
				end_spinodal = lambdalist[k]
				count+=1
	return np.array([start_spinodal,end_spinodal])

#Calculates coexistence region for general psi_0 (still must be constant):
def Coexist(psi0,zeta):
	spinodalpoints = Spinodal(psi0,zeta)
	common_t = pso(functiontominimize,[0.5,spinodalpoints[1]-0.000001],[spinodalpoints[0],1],args=(psi0,zeta))[0]    
	#common_t = pso(functiontominimize,[0.5,zeta**(-1/3)-0.00001],[zeta**(-1/3),1],args=(psi0,zeta))[0]    

	return common_t[0],common_t[1]


def CoexistenceRegion(psi0,zetalist):
	Lambda_Low=[]
	Lambda_High=[]

	for zeta in zetalist:
		low,high = Coexist(psi0,zeta)
		Lambda_Low.append(low)
		Lambda_High.append(high)

	return np.array(Lambda_Low),np.array(Lambda_High)

def CoexistenceRegion0(zetalist):
	Lambda_Low=[]
	Lambda_High=[]

	for zeta in zetalist:
		low,high = Coexist_0(zeta)
		Lambda_Low.append(low)
		Lambda_High.append(high)

	return np.array(Lambda_Low),np.array(Lambda_High)


def Stress(zetalist,Coexist_Start,Coexist_End,psi0):
	N = len(zetalist)
	sigmalist=np.zeros(N)
	sigmalist[:] = np.NaN

	for i in range(N):
		zeta = zetalist[i]
		fstart = f_coexist(Coexist_Start[i],zeta,psi0)
		fend = f_coexist(Coexist_End[i],zeta,psi0)
		sigmalist[i] = (fstart-fend)/(Coexist_Start[i]-Coexist_End[i])

	return sigmalist

def Stress0(zetalist,Coexist_Start,Coexist_End):
	N = len(zetalist)
	sigmalist=np.zeros(N)
	sigmalist[:] = np.NaN

	for i in range(N):
		zeta = zetalist[i]
		fstart = f_coexist(Coexist_Start[i],zeta,0)
		fend = f_coexist(Coexist_End[i],zeta,psi0)
		sigmalist[i] = (fstart-fend)/(Coexist_Start[i]-Coexist_End[i])

	return sigmalist

def Stress0(zetalist,Coexist_Start,Coexist_End):
	N = len(zetalist)
	sigmalist=np.zeros(N)
	sigmalist[:] = np.NaN

	for i in range(N):
		zeta = zetalist[i]
		sigma = -1/Coexist_End[i]**2 + Coexist_End[i]
		sigmalist[i] = sigma

	return sigmalist


#Parameters:
N=10
lambdalist=np.linspace(0.8,1.0,num=N)
zetalist=np.linspace(1.01,1.8,num=N)
#psi0list = [0.1,0.2,0.3]
#colorlist = ['blue','green','xkcd:orange']


# Figure Setup:
fig=plt.figure()
gs=gridspec.GridSpec(2,1,width_ratios=[1],height_ratios=[1,1])
ax1=plt.subplot(gs[0])
ax2=plt.subplot(gs[1])
ax1.minorticks_on()
ax2.minorticks_on()

# Plot Coexistence Regions:
Lambda_Low,Lambda_High = CoexistenceRegion0(zetalist)
ax1.plot(100*(Lambda_High-1),zetalist,color='black',lw=1.5,ls='-')
ax1.plot(100*(Lambda_Low-1),zetalist,color='black',lw=1.5,ls='-')

sigmalist = Stress0(zetalist,Lambda_High,Lambda_Low)
ax2.plot(sigmalist,zetalist,label='$\psi_0 = 0$',color='black',lw=1.5,ls='-')


psi0 = 0.1
Lambda_Low,Lambda_High = CoexistenceRegion(psi0,zetalist)
ax1.plot(100*(Lambda_High-1),zetalist,color='blue',lw=2,ls='-.')
ax1.plot(100*(Lambda_Low-1),zetalist,color='blue',lw=2,ls='-.')

sigmalist = Stress(zetalist,Lambda_High,Lambda_Low,psi0)
ax2.plot(sigmalist,zetalist,label='$\psi_0 = 0.1$',color='blue',lw=2,ls='-.')


psi0 = 0.2
Lambda_Low,Lambda_High = CoexistenceRegion(psi0,zetalist)
ax1.plot(100*(Lambda_High-1),zetalist,color='green',lw=2.5,ls=':')
ax1.plot(100*(Lambda_Low-1),zetalist,color='green',lw=2.5,ls=':')

sigmalist = Stress(zetalist,Lambda_High,Lambda_Low,psi0)
ax2.plot(sigmalist,zetalist,label='$\psi_0 = 0.2$',color='green',lw=2.5,ls=':')


psi0 = 0.3
Lambda_Low,Lambda_High = CoexistenceRegion(psi0,zetalist)
ax1.plot(100*(Lambda_High-1),zetalist,color='xkcd:orange',lw=3,ls='--')
ax1.plot(100*(Lambda_Low-1),zetalist,color='xkcd:orange',lw=3,ls='--')

sigmalist = Stress(zetalist,Lambda_High,Lambda_Low,psi0)
ax2.plot(sigmalist,zetalist,label='$\psi_0 = 0.3$',color='xkcd:orange',lw=3,ls='--')

'''
psi0 = 0.4
Lambda_Low,Lambda_High = CoexistenceRegion(psi0,zetalist)
ax1.plot(100*(Lambda_High-1),zetalist,color='red',lw=3.5,ls='-')
ax1.plot(100*(Lambda_Low-1),zetalist,color='red',lw=3.5,ls='-')

sigmalist = Stress(zetalist,Lambda_High,Lambda_Low,psi0)
ax2.plot(sigmalist,zetalist,label='$\psi_0 = 0.4$',color='red',lw=3.5,ls='-')
'''



# Axes, Labels, and Aesthetic Stuff:

ax1.text(-13,1.5,'Coexistence',fontsize=16)
ax1.text(-3,1.6,r'$\epsilon_L$',fontsize=20)
ax1.text(-15,1.2,r'$\epsilon_H$',fontsize=20)

ax1.set_title('A)',loc='left',fontsize=20)
ax1.set_ylabel('$\zeta$',fontsize=20,rotation=0,labelpad=15)
ax1.set_xlabel('Fibril Strain',fontsize=16)
#ax1.legend(loc='best',fontsize=16)
ax1.set_ylim(1,1.8)
ax1.set_xlim(-20,0)

ax1.set_xticks([-20,-15,-10,-5,0])
ax1.set_xticklabels(['$-20\%$','$-15\%$', '$-10\%$','$-5\%$', '$0\%$'],fontsize=14)

ax2.set_title('B)',loc='left',fontsize=20)
ax2.set_ylabel('$\zeta$',fontsize=20,rotation=0,labelpad=15)
ax2.set_xlabel('Stress $\sigma/\mu$',fontsize=16)
ax2.legend(loc='best',fontsize=16)
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


