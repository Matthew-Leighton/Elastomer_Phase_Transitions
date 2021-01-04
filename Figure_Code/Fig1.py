import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from matplotlib import rc
import matplotlib.gridspec as gridspec
import scipy.optimize as opt
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

#Calculates free energy density
def f(strain,zeta,psi0):
    psi = ExactPsi(psi0,strain,zeta)
    fe = 1/strain + (1/strain)*(1+(zeta-1)*np.sin(psi0)**2)*(1+(1/zeta -1)*np.sin(psi)**2)
    fe += (strain**2)*(1+(zeta-1)*np.cos(psi0)**2)*(1+(1/zeta -1)*np.cos(psi)**2)
    fe += np.sqrt(strain)*(2 - zeta - 1/zeta)*np.sin(2*psi0)*np.sin(2*psi)/2
    return fe/2

#Analytic result for psi when psi_0=0:
def psi_00(strainlist,zeta):
	psi_list=[]
	for strain in strainlist:
		if strain>= zeta**(-1/3):
			psi_list.append(0)
		else:
			psi_list.append(np.pi/2)
	return psi_list

#Analytic result for f when psi_0=0:
def f0(strainlist,zeta):
	f0_list=[]
	for strain in strainlist:
		if strain>= zeta**(-1/3):
			f0_list.append(strain**2 /2 + 1/strain)
		else:
			f0_list.append(((1+1/zeta)/strain + zeta*strain**2)/2)
	return f0_list



##### Coexistence Region Calculations:

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
	lambdalist = np.linspace(0.5,1,num=1000)
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
	common_t = pso(functiontominimize,[0.7,spinodalpoints[1]-0.000001],[spinodalpoints[0],1],args=(psi0,zeta))[0]    
	#common_t = pso(functiontominimize,[0.5,zeta**(-1/3)-0.00001],[zeta**(-1/3),1],args=(psi0,zeta))[0]    

	return common_t[0],common_t[1]






# Plots Figure 1 in the Manuscript
def PlotFig1():

	#Parameters
	zeta = 1.3
	N=1000
	lambdalist = np.linspace(0.8,1.1,num=N)
	epsilonlist = lambdalist-1

	# Figure Setup
	fig=plt.figure()
	gs=gridspec.GridSpec(2,1,width_ratios=[1],height_ratios=[1,1])
	ax1=plt.subplot(gs[0])
	ax2=plt.subplot(gs[1])
	ax1.minorticks_on()
	ax2.minorticks_on()

	# Plot Twist Functions:
	ax1.plot(100*epsilonlist, psi_00(lambdalist,zeta),label = '$\psi_0 = 0$',lw=2,ls='-',color='black')
	ax1.plot(100*epsilonlist, ExactPsi(0.1,lambdalist,zeta),label = '$\psi_0 = 0.1$',lw=2.5,ls='-.',color='blue')
	ax1.plot(100*epsilonlist, ExactPsi(0.3,lambdalist,zeta),label = '$\psi_0 = 0.3$',lw=3,ls='--',color='xkcd:orange')

	# Plot free energy densities:
	ax2.plot(100*epsilonlist, f0(lambdalist,zeta),label = '$\psi_0 = 0$',lw=2,ls='-',color='black')
	ax2.plot(100*epsilonlist, f(lambdalist,zeta,0.1),label = '$\psi_0 = 0.1$',lw=2.5,ls='-.',color='blue')
	ax2.plot(100*epsilonlist, f(lambdalist,zeta,0.3),label = '$\psi_0 = 0.3$',lw=3,ls='--',color='xkcd:orange')

	# Plot Coexistence Regions:
	lambda_low_0,lambda_high_0 = Coexist_0(zeta)
	ax2.hlines(1.525,100*(lambda_low_0-1),100*(lambda_high_0-1),lw=2,ls='-',color='black')

	lambda_low_0,lambda_high_0 = Coexist(0.1,zeta)
	ax2.hlines(1.524,100*(lambda_low_0-1),100*(lambda_high_0-1),lw=2.5,ls='-.',color='blue')

	lambda_low_0,lambda_high_0 = Coexist(0.3,zeta)
	ax2.hlines(1.523,100*(lambda_low_0-1),100*(lambda_high_0-1),lw=3,ls='--',color='xkcd:orange')



	# Axes, Labels, and Aesthetic Stuff:

	ax1.set_title('A)',loc='left',fontsize=20)
	ax1.set_ylabel('$\psi$',fontsize=20,rotation=0,labelpad=15)
	ax1.legend(loc='best',fontsize=16)
	ax1.set_ylim(-0.1,1.7)
	ax1.set_xlim(-20,10)
	ax1.set_xticks([-20,-15,-10,-5,0,5,10])
	ax1.set_xticklabels(['$-20\%$','$-15\%$', '$-10\%$','$-5\%$', '$0\%$','$5\%$','$10\%$'],fontsize=14)


	ax2.set_title('B)',loc='left',fontsize=20)
	ax2.set_ylabel('$\\frac{f}{\mu}$',fontsize=30,rotation=0,labelpad=15)
	ax2.set_xlabel('Fibril Strain',fontsize=16)
	ax2.set_xlim(-20,10)
	ax2.set_xticks([-20,-15,-10,-5,0,5,10])
	ax2.set_xticklabels(['$-20\%$','$-15\%$', '$-10\%$','$-5\%$', '$0\%$','$5\%$','$10\%$'],fontsize=14)

	#ax1.tick_params(axis='x', labelsize=14)
	ax1.tick_params(axis='y', labelsize=14)
	#ax2.tick_params(axis='x', labelsize=16)
	ax2.tick_params(axis='y', labelsize=14)

	plt.tight_layout(pad=0.5)
	plt.show()


PlotFig1()
