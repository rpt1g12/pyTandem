import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
from scipy import stats
import getpass
user=getpass.getuser()
pi=np.pi
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')
figs=[];axs=[];hdls=[];lbls=[];lgds=[];nfig=-1;fs=18
#%%
fmode = 1;
save=False
prandtl=True #Apply Prandtl-Glauert correction?
mode='high'
AoA=6 #Angle of attack
A=0 #WLE amplitude
nwave=1 #Number of WLE wavelengths
if A>0:
    wavy=True
else:
    wavy=False
    
if AoA==0:
    if wavy:
        ts,te=30,45 #Time interval considered
    else:
        ts,te=60,90
elif AoA==6:
    if wavy:
        ts,te=100,135 #Time interval considered
    else:
        ts,te=160,210
elif AoA==10:

    if wavy:
        ts,te=120,213 #Time interval considered
    else:
        if nwave==8:        
            ts,te=355,368
        elif nwave==1:
            ts,te=270,330
#%%Simulation time history
folder='clData/lowAoA/'
sim="{:1d}A{:02d}W11AoA{:02d}".format(nwave,A,AoA)
dataset='/home/'+user+'/anaconda3/pyTandem/'+folder+sim+'.dat'
spath='/home/rpt1g12/Documents/thesis/data/lowAoA/ClCdPSD/'
n,t0,clin,cdin,taoa,tmach=np.loadtxt(dataset,skiprows=1,unpack=True)
M_inf=tmach[-1]
if prandtl:
    beta=np.sqrt(1-M_inf**2)
else:
    beta=1
#%%
f,a=getFig('ClCd_pv');figs.append(f),axs.append(a);nfig+=1
folderPV=folder+'/pvClCd/'

datasetPV='/home/'+user+'/anaconda3/pyTandem/'+folderPV+sim+'.dat';
n,t0,clpin,clvin,cdpin,cdvin=np.loadtxt(datasetPV,skiprows=1,unpack=True)

if ts<0: ts=min(t0) 
if te<0: te=max(t0) 
    
ns=np.where((t0>=ts))[0][0]
ne=np.where((t0<=te))[0][-1]
t=t0[ns:ne]
clp=clpin[ns:ne]*beta;cdp=cdpin[ns:ne]*beta
clv=clvin[ns:ne]*beta;cdv=cdvin[ns:ne]*beta
#%% Plot histories
cl=clp+clv;cd=cdp+cdv

if fmode==0:
    print('lift mode!')
    axs[nfig].plot(t*0.3,clp,label='r$C_{lp}$')
    axs[nfig].plot(t*0.3,clv,label='r$C_{lv}$')
    axs[nfig].set_ylabel(r'$C_l$',fontsize=18)
    axs[nfig].plot(t*0.3,clp+clv,label=r'$C_l$',lw=2,linestyle='--')
elif fmode==1:
    print('drag mode!')
    axs[nfig].plot(t*0.3,cdp,lw=2,label=r'$C_{dp}$')
    axs[nfig].plot(t*0.3,cdv,lw=2,label=r'$C_{dv}$')
    axs[nfig].set_ylabel(r'$C_d$',fontsize=18)
    axs[nfig].plot(t*0.3,cdp+cdv,label=r'$C_d$',lw=2,linestyle='--')
    
axs[nfig].set_xlabel(r'$t$',fontsize=18)

#%%
nfigs=len(figs)
for i in range(nfigs):
    fit(axs[i])
    hdl,lbl,lgd=getLabels(ax=axs[i],ncol=3,fontsize=18)
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
#%%
if save:
    savePlotFile(ax=axs[i],path=spath,sameX=True)