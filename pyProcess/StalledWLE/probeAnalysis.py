# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import welch as psdw
from scipy.signal import hilbert
#import pywt #Wavelets
import matplotlib.pyplot as plt;plt.close('all')
from lib.stats import *
from lib.myPlots import *
from lib.myPlot3dOperator import *
import lib.Plot3DClasses as p3d
from scipy import stats
from scipy.interpolate import griddata
pi=np.pi
import getpass
user=getpass.getuser()
#from lib.matplotlib2tikz import save as tikz_save
figs=[];axs=[];hdls=[];lbls=[];lgds=[];nfig=-1;fs=18
plt.close('all')
import importlib
#%%
importlib.reload(p3d)
#%% Options
wavelets=False;plot=False;save=True
#%%Simulation Name and set-up
#theta=1;Ue=1
#hName='8SLE';probR=range(0,17,1)
#theta=0.000639896428216028;Ue=0.4058528675001842
#tpath='/media/{:}/082F-63FE/phd/thesis/data/StalledWLE/fig24_8SLE/'.format(user)
#hName='8SLE_150';probR=range(0,17,1)
#theta=0.000639896428216028;Ue=0.4058528675001842
#tpath='/media/{:}/082F-63FE/phd/thesis/data/StalledWLE/fig24_8SLE/'.format(user)
#hName='8WLE_Central';probR=range(0,17,1)
#theta=0.0006705586212473073;Ue=0.4565491406912271
#tpath='/media/{:}/082F-63FE/phd/thesis/data/StalledWLE/fig25_8WLE/'.format(user)
hName='8WLE_SSL';probR=range(0,17,1)
theta=0.0006922426498962798;Ue=0.3807581073978518
tpath='/media/{:}/082F-63FE/phd/thesis/data/StalledWLE/fig25_8WLE/'.format(user)
#%%
theta=1*np.sin(np.deg2rad(20));Ue=0.3
#%%Probes time history

var='p'
dataset='/home/'+user+'/anaconda3/pyTandem/pgfPlots/'+hName+'_Hist_{}.dat'.format(var);
dummy=np.loadtxt(dataset,skiprows=1,unpack=True)
nprob=dummy.shape[0]-1

#%%
_,t,nsam,fsam=rsample(dummy[1,:],dummy[0,:],verbose=True,force=True)
data=np.zeros((nprob,nsam))
for pb in range(nprob):
    data[pb,:],t,nsam,fsam=rsample(dummy[pb+1,:],dummy[0,:],verbose=True,force=True)
    data[pb,:]=rmvLS(t,data[pb,:],1)

del dummy
#%%
scale=False;sclg='density'
nw=32;ovlp=0.5;
nseg,novlp,ntt,fmax,fmin=defWin(t,data[0,:],nw,ovlp,verbose=False)
pdata=np.zeros((len(probR),int(nseg/2+1)))
for i in range(len(probR)):
    pb=probR[i]
    ff,pdata[i,:]=psdw(data[pb,:],fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
    if plot:
        f,a=getFig('{}Probe{:d}'.format(var,pb));nfig+=1
        figs.append(f);axs.append(a)
        axs[nfig].plot(t*0.3/np.sin(np.deg2rad(20)),data[pb,:])
ff=ff*theta/Ue    
#%%
f,a=getFig('{}Freq'.format(var));nfig+=1
figs.append(f);axs.append(a)
for i in range(len(probR)):
    pb=probR[i]
    axs[nfig].loglog(ff,pdata[i,:],label='p{:02d}'.format(pb),lw=2)

hdl,lbl,lgd=getLabels(ax=axs[nfig],ncol=2,fontsize=fs,loc='lower left')
#%%
if save:
    #savePlotFile(ax=axs[nfig],path=tpath,name='BLscaled_{}.dat'.format(hName))
    savePlotFile(ax=axs[nfig],path=tpath,name='{}.dat'.format(hName))
#%%
if wavelets:
    #plt.close('all')
    nlv=4;wavetype='db4'
    sgnl=data[0,:]
    c=pywt.wavedec(sgnl,wavetype,level=nlv)
    f,a=getFig('Wavelets');nfig+=1
    figs.append(f);axs.append(a)
    A=np.zeros((nlv+1,nsam))
    ad='a'
    ar=[nlv]+[i for i in range(nlv,0,-1)]
    for i in range(nlv+1):
        ii=ar[i]
        A[i]=pywt.upcoef(ad,c[i],wavetype,level=ii,take=nsam);ad='d'
        

    #%%
    D=A[1,:]
    for i in range(2,nlv+1):
        D+=A[i,:]
    
    axs[nfig].plot(t,sgnl,lw=1,label='S')
    axs[nfig].plot(t,A[0,:],lw=2,label='A')
    #axs[nfig].plot(t,D,lw=2,label='D')
    hdl,lbl,lgd=getLabels(ax=axs[nfig],ncol=2,fontsize=fs,loc='best')
    
#    #%% Apply hilbert transform
#    f,a=getFig('freq(t)');nfig+=1
#    figs.append(f);axs.append(a)
#    
#    hD=hilbert(sgnl)
#    E=np.abs(hD)
#    iph = np.unwrap(np.angle(hD))
#    ifreq = (np.diff(iph)/(2.0*np.pi)*fsam)
#    
#    axs[nfig].plot(t[1:],ifreq)


