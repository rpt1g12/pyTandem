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
save=False
prandtl=True #Apply Prandtl-Glauert correction?
mode='high'
AoA=10 #Angle of attack
A=15 #WLE amplitude
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
    if mode=='low':
        nw=1;ovlp=0.5
    else:
        nw=3;ovlp=0.5
elif AoA==6:
    if wavy:
        ts,te=100,135 #Time interval considered
    else:
        ts,te=160,210
    if mode=='low':
        nw=1;ovlp=0.5
    else:
        nw=8;ovlp=0.5
elif AoA==10:
    if mode=='low':
        nw=2;ovlp=0.5
    else:
        nw=16;ovlp=0.5
    if wavy:
        ts,te=120,213 #Time interval considered
    else:
        if nwave==8:        
            ts,te=355,368
            if mode!='low':
                nw=5;ovlp=0.5
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

if ts<0: ts=min(t0) 
if te<0: te=max(t0) 
    
ns=np.where((t0>=ts))[0][0]
ne=np.where((t0<=te))[0][-1]
t=t0[ns:ne]
cl=clin[ns:ne]*beta;cd=cdin[ns:ne]*beta
aoa=taoa[ns:ne]

clmean,cdmean,clvar,cdvar=np.mean(cl),np.mean(cd),np.var(cl),np.var(cd)


#%%
f,a=getFig('ClCd');figs.append(f),axs.append(a);nfig+=1

axs[nfig].plot(t*0.3,cl,label='cl')
axs[nfig].plot(t*0.3,cd,label='cd')

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
#%% Resample signals and Windowing
scale=False; step=False; sclg='density'

clp,tn,nsam,fsam=rsample(clp,t,verbose=True,rmAvg=True)
clv,tn,nsam,fsam=rsample(clv,t,verbose=True,rmAvg=True)
cdp,tn,nsam,fsam=rsample(cdp,t,verbose=True,rmAvg=True)
cdv,tn,nsam,fsam=rsample(cdv,t,verbose=True,rmAvg=True)
t=tn
nseg,novlp,ntt,fmax,fmin=defWin(t,clp,nw,ovlp,verbose=True)

#%% Plot histories

axs[nfig].plot(t*0.3,clp,label='Clp')
axs[nfig].plot(t*0.3,cdp,label='Cdp')
axs[nfig].plot(t*0.3,clv,label='Clv')
axs[nfig].plot(t*0.3,cdv,label='Cdv')
cl=clp+clv;cd=cdp+cdv

axs[nfig].plot(t*0.3,clp+clv,label=r'$C_l$',linestyle='--')
axs[nfig].plot(t*0.3,cdp+cdv,label=r'$C_d$',linestyle='--')

#%%
ff,Pclp=psdw(clp,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
ff,Pclv=psdw(clv,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
ff,Pcdp=psdw(cdp,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
ff,Pcdv=psdw(cdv,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
ff,Pcl=psdw(cl,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
ff,Pcd=psdw(cd,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)

#%% Plot frequencies
#st=ff*np.sin(np.deg2rad(AoA))/0.3
st=ff/0.3
if mode=='low':
    nmax=np.where(st>10)[0][0]
else:
    nmax=len(st)
f,a=getFig(sim+mode);figs.append(f),axs.append(a);nfig+=1
axs[nfig].loglog(st[:nmax],Pclp[:nmax],label='Clp')
axs[nfig].loglog(st[:nmax],Pcdp[:nmax],label='Cdp')
axs[nfig].loglog(st[:nmax],Pclv[:nmax],label='Clv')
axs[nfig].loglog(st[:nmax],Pcdv[:nmax],label='Cdv')


axs[nfig].loglog(st[:nmax],Pcl[:nmax],label=r'$C_l$',linestyle='--')
axs[nfig].loglog(st[:nmax],Pcd[:nmax],label=r'$C_d$',linestyle='--')
#%%
nfigs=len(figs)
for i in range(nfigs):
    fit(axs[i])
    hdl,lbl,lgd=getLabels(ax=axs[i],ncol=2,fontsize=15)
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
#%%
if save:
    savePlotFile(ax=axs[i],path=spath,sameX=True)