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
prandtl=True #Apply Prandtl-Glauert correction?
AoA=6 #Angle of attack
A=0 #WLE amplitude
nwave=1 #Number of WLE wavelengths
if A>0:
    wavy=True
else:
    wavy=False
ts,te='min','max' #Time interval considered

#%%Simulation time history
folder='clData/lowAoA/'
sim="{:1d}A{:02d}W11AoA{:02d}.dat".format(nwave,A,AoA)
dataset='/home/'+user+'/anaconda3/pyTandem/'+folder+sim
n,t0,clin,cdin,taoa,tmach=np.loadtxt(dataset,skiprows=1,unpack=True)
M_inf=tmach[-1]
if prandtl:
    beta=np.sqrt(1-M_inf**2)
else:
    beta=1
#%%Expected values
dataset='/home/'+user+'/anaconda3/pyTandem/clData/HansenClCd.dat';
aoah,clh,cdh,clh0,cdh0=np.loadtxt(dataset,skiprows=1,unpack=True)
fcl,fcd,fcl0,fcd0=interf(aoah,clh),interf(aoah,cdh),interf(aoah,clh0),interf(aoah,cdh0)
if wavy:
    clh,cdh=fcl.f(AoA),fcd.f(AoA)
else:
    clh,cdh=fcl0.f(AoA),fcd0.f(AoA)
clh,cdh=float(clh),float(cdh)
print(r'Expected Values: (cl,cd)=({:1.3f},{:1.3f})'.format(clh,cdh))
#%%

if ts=='min': ts=min(t0) 
if te=='max': te=max(t0) 
    
ns=np.where((t0>=ts))[0][0]
ne=np.where((t0<=te))[0][-1]
t=t0[ns:ne]
cl=clin[ns:ne]*beta;cd=cdin[ns:ne]*beta
aoa=taoa[ns:ne]

clmean,cdmean,clvar,cdvar=np.mean(cl),np.mean(cd),np.var(cl),np.var(cd)


#%%
f,a=getFig(sim);figs.append(f),axs.append(a);nfig+=1
axs[nfig].plot(t,cl,lw=2,color='blue',label=r'$C_L$')
axs[nfig].plot(t,cd,lw=2,color='red',label=r'$C_D$')
axs[nfig].axhline(y=clh,color='blue',linewidth=2,linestyle='--',label=r'$C_{L_{Hansen}}$') 
axs[nfig].axhline(y=cdh,color='red',linewidth=2,linestyle='--',label=r'$C_{D_{Hansen}}$')
axs[nfig].set_xlabel(r'$t^*$',fontsize=fs)
axs[nfig].set_ylabel(r'$C_L\, &\, C_D$',fontsize=fs)

#%%
nfigs=len(figs)
for i in range(nfigs):
    fit(axs[i])
    hdl,lbl,lgd=getLabels(ax=axs[i],ncol=2,fontsize=15)
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
#%%
axs[nfig].set_xlim(ts,te)

print('\nCl_mean={:2.4f} Cd_mean={:2.4f} Cl_var={:1.4e} Cd_var={:1.4e}'.format(clmean,cdmean,clvar,cdvar))