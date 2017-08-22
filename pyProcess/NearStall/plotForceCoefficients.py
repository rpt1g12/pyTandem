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
save=True
prandtl=True #Apply Prandtl-Glauert correction?
AoA=10 #Angle of attack
A=15 #WLE amplitude
nwave=8 #Number of WLE wavelengths
if A>0:
    wavy=True
else:
    wavy=False
ts,te='min','max' #Time interval considered
#%%
folder='clData/aoa10/'
sim="{:1d}A{:02d}W11AoA{:02d}".format(nwave,A,AoA)
dataset='/home/'+user+'/anaconda3/pyTandem/'+folder+sim+'.dat';
spath='/home/rpt1g12/Documents/thesis/data/nearStall/clcdHistory/'
n,t0,clin,cdin,taoa,tmach=np.loadtxt(dataset,skiprows=1,unpack=True)
if ts=='min': ts=min(t0) 
if te=='max': te=max(t0) 
M_inf=tmach[-1]
if prandtl:
    beta=np.sqrt(1-M_inf**2)
else:
    beta=1
#%%
ns=np.where((t0>=ts))[0][0]
ne=np.where((t0<=te))[0][-1]
cl,t,nsam,fsam=rsample(clin[ns:ne]*beta,t0[ns:ne],verbose=False,force=True,nsample=2048)
cd,t,nsam,fsam=rsample(cdin[ns:ne]*beta,t0[ns:ne],verbose=False,force=True,nsample=2048)
aoa,t,nsam,fsam=rsample(taoa[ns:ne],t0[ns:ne],verbose=False,force=True,nsample=2048)
t-=t[0]
#%%
if wavy:
    fname='{:1d}WLE{:02d}'.format(nwave,AoA)
else:
    fname='{:1d}SLE{:02d}'.format(nwave,AoA)
f,a=getFig('t_cl_'+fname);figs.append(f),axs.append(a);nfig+=1
axs[nfig].plot(t,cl,lw=2,color='blue',label=r'$C_L$')
axs[nfig].plot(t,cd,lw=2,color='red',label=r'$C_D$')
axs[nfig].plot(t,cl/cd,lw=2,color='green',label=r'$C_L/C_D$')
axs[nfig].axvline(x=t[0],color='black',linewidth=2,linestyle='--') 
axs[nfig].set_xlabel(r'$t^*$',fontsize=fs)
axs[nfig].set_ylabel(r'$C_L\, &\, C_D$',fontsize=fs)
#%%
f,a=getFig('t_aoa_'+fname);figs.append(f),axs.append(a);nfig+=1
axs[nfig].plot(t,aoa,lw=2,color='blue',label=r'$\alpha$')
#%%
if len(np.where(aoa>19.99999)[0])>1:
    n0=np.where(aoa>19.99999)[0][0]
    n1=np.where(aoa>19.99999)[0][-1]
else:
    n0=-1
    n1=0
f,a=getFig('aoa_cl_'+fname);figs.append(f),axs.append(a);nfig+=1
axs[nfig].plot(aoa[:n0],cl[:n0],lw=2,color='blue',label='C_Lup')
axs[nfig].plot(aoa[:n0],cd[:n0],lw=2,color='red',label='C_Dup')
axs[nfig].plot(aoa[:n0],cl[:n0]/cd[:n0],lw=2,color='green',label='C_L/C_Dup')
axs[nfig].plot(aoa[n1:],cl[n1:],lw=2,color='blue',linestyle='--',label='C_Ldown')
axs[nfig].plot(aoa[n1:],cd[n1:],lw=2,color='red',linestyle='--',label='C_Ddown')
axs[nfig].plot(aoa[n1:],cl[n1:]/cd[n1:],lw=2,color='green',linestyle='--',label='C_L/C_Ddown')
#%%
nfigs=len(figs)
for i in range(nfigs):
    fit(axs[i])
    hdl,lbl,lgd=getLabels(ax=axs[i],ncol=2,fontsize=15)
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
#%% Saving
    if save:
        savePlotFile(path=spath,ax=axs[i])
        