import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.stats import rsample
from lib.myPlots import *
from lib.myclasses import *
from scipy import stats
import getpass
user=getpass.getuser()
pi=np.pi
plt.close('all')
figs=[];axs=[];hdls=[];lbls=[];lgds=[];nfig=-1;fs=18
#%%
save=False;
average=True
M_inf=0.3 #Mach number
AoA=0 #Angle of attack
nwave=1 #Number of WLE wavelengths
spath='/home/rpt1g12/Documents/thesis/data/lowAoA/pint{}/'.format(AoA)
#%%Create first figure
fname='{:1d}pintCl'.format(nwave)
f,a=getFig(fname);nfig+=1 # Get figure and axis objects
figs.append(f);axs.append(a) # Append them to the figures and axes arrays
fname='{:1d}pintCd'.format(nwave)
f,a=getFig(fname);nfig+=1 # Get figure and axis objects
figs.append(f);axs.append(a) # Append them to the figures and axes arrays
for A in [0,15]:
    if A>0:
        wavy=True
        prandtl=True #Apply Prandtl-Glauert correction?
        sName='WLE'
        ls='--'
        factor=1
    else:
        wavy=False
        prandtl=True #Apply Prandtl-Glauert correction?
        sName='SLE'
        ls='-'
        factor=(-1)
    #%%Load piecewise integration from Paraview
    folder='piecewiseIntegration/AoA{:02d}/'.format(AoA)
    sim="{:1d}{}_pint.dat".format(nwave,sName)
    #sim="{:1d}{}_pint0.dat".format(nwave,sName)
    dataset='/home/'+user+'/anaconda3/pyTandem/'+folder+sim
    #dataset='/home/'+user+'/Desktop/paraviewExports/'+sim
    xm,xl,cl,cd=np.loadtxt(dataset,skiprows=1,unpack=True)
    
    nx=len(xm)
    #%%Scale coefficients with prandtl rule
    if prandtl:
        beta=np.sqrt(1-M_inf**2)
    else:
        beta=1
    cl*=beta;cd*=beta
#%%Average every 7 steps
    if average:
        for i in range(3,nx,7):
            cl[i-3:i+4]=cl[i-3:i+4].mean()
            cd[i-3:i+4]=cd[i-3:i+4].mean()
    
    cl/=xl;cd/=xl #And divide by Delta_X
    
    #%% Compute deltas
    if A>0:
        deltaCl+=cl
        deltaCd+=cd
    else:
        deltaCl=-cl
        deltaCd=-cd
    #%%Plotting        
    axs[0].plot(xm,cl,lw=2,color='blue',linestyle=ls,label='{}cl'.format(sName))
    axs[1].plot(xm,cd,lw=2,color='red',linestyle=ls,label='{}cd'.format(sName))
#%%Plot delta
fname='{:1d}pintClDelta'.format(nwave)
f,a=getFig(fname);nfig+=1 # Get figure and axis objects
figs.append(f);axs.append(a) # Append them to the figures and axes arrays
fname='{:1d}pintCdDelta'.format(nwave)
f,a=getFig(fname);nfig+=1 # Get figure and axis objects
figs.append(f);axs.append(a) # Append them to the figures and axes arrays

axs[2].plot(xm,deltaCl,lw=2,color='blue',label='deltaCl')
axs[3].plot(xm,deltaCd,lw=2,color='red',label='deltaCd')
#%%
for i in range(nfig+1):    
    hdl,lbl,lgd=getLabels(axs[i],fontsize=fs,ncol=1)
    hdls.append(hdl)
    lbls.append(lbl)
    lgds.append(lgd)
    fit(axs[i])
    #%% Saving
    if save:
        savePlotFile(path=spath,ax=axs[i])

#%%Clean variables
del deltaCl,deltaCd