# -*- coding: utf-8 -*-
import os
import math
import numpy as np
from scipy.signal import spectrogram 
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
from lib.myPlot3dOperator import *
import lib.Plot3DClasses as p3d
from scipy import stats
from scipy.interpolate import griddata
from scipy.signal import argrelextrema as locmaxmin
from matplotlib import colors,ticker
pi=np.pi
import getpass
user=getpass.getuser()
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')
figs=[];axs=[];hdls=[];lbls=[];lgds=[];nfig=-1;fs=18
import importlib
styles=['-','--',':','-.']
colList=['blue','red','cyan','magenta']
importlib.reload(p3d)

#%%
save=False
cutoff=1.5
#%% Sub-set settings
ssName='upSurface'
bkXYZ=[1,1,1]
blocks=[4]
xRanges=[range(110)]
yRanges=[[1]]
zRanges=[None]
vnames=['r','u','v','w','p']; # Variable names
sPattern='solT*.*.q'
var='w'    
#%% Paths set-up
path="/media/{}/dellHDD/post/4A15W11AoA10/heaving/{}/".format(user,ssName)
print('Reading data from:\n {}'.format(path))
#%%
gfile='grid.xyz' # Grid file name
files=p3d.getFileNames(path=path,pattern=sPattern)
sfile=files[0]#'solTavgCf+tw+Cp.q'# Solution file name

fl=p3d.flow(path,"grid.xyz",sfile) # Check-out flow object
fl.rdHdr() # Read header from grid file and set-up flow object properties
fl.rdGrid() # Read Grid
fl.rdSol(vnames=vnames) # Read solution file
flInfo=fl.rdFlowInfo() # Read solution file info, i.e. Mach, AoA, Re and time
#%%
#fl.contourf(varname='w',vmin=-0.05,vmax=0.05,plane=1,nlvl=21,cmap=plt.cm.bwr)
nt=len(files)
t=np.zeros(nt)
i=0
data=np.zeros((4,nt))
for sfile in files:
    fl.rdSol(vnames=vnames,sfile=sfile)
    flInfo=fl.rdFlowInfo()
    data[:,i]=fl.blk[0].var[var].getValues()[75,0,range(12,12+48*3+1,48)]
    t[i]=flInfo[3]
    i+=1
    
#%% Resample
deltaT=50
fsam=64
t0=t[0]
t=t-t0
t1=deltaT
n=np.where(t>t1)[0][0]
tnew=np.linspace(0,deltaT,deltaT*fsam+1)
nsam=len(tnew)
Data=np.zeros((4,nsam))
Data_filt=Data.copy()
for i in range(4):
    Data[i,:],tt,nsam,fsam=rsample(data[i,:n],t[:n],verbose=False,tnew=tnew)
    Data_filt[i,:]=fourierFilter(Data[i,:],'low',fsam,cutoff*0.3,15,False)
#%%
f,a=getFig('{}History'.format(var));figs.append(f),axs.append(a);nfig+=1
for i in [0,1]:
    a.plot(0.3*tt,Data[i,:],color=colList[i],lw=2,label='T{:1d}'.format(i+1))
#%%
f,a=getFig('{}HistoryFilt'.format(var));figs.append(f),axs.append(a);nfig+=1
for i in [0,1]:
    a.plot(0.3*tt,Data_filt[i,:],color=colList[i],lw=2,label='T{:1d}_filt'.format(i+1))

#%% Define windows
vmin=1e-9;vmax=1e-5
nw=24;ovlp=0.80;sclg='density'

sgnl=Data_filt[0,:]
sgnl=Data[0,:]

nseg,novlp,ntt,fmax,fmin=defWin(t,sgnl,nw,ovlp,verbose=False)


#%% Perform STFT
ff2,tt2,psd=spectrogram(sgnl,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)

#%% Perform FFT
ff3,PP=psdw(sgnl,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
#%% Plot frequencies
st=ff3/0.3

f,a=getFig('PSD');figs.append(f),axs.append(a);nfig+=1
axs[nfig].loglog(st[:],PP[:],label='PSD')

#%% Plot Spectogram
TT,FF=np.meshgrid(tt2,ff2/0.3)
f,a=getFig('Spectogram_'+var);figs.append(f),axs.append(a);nfig+=1
im=a.pcolor(TT,FF,psd,norm=colors.LogNorm(vmin=psd.min(), vmax=psd.max()),cmap='jet',vmin=vmin,vmax=vmax)
cb=f.colorbar(im,orientation='vertical')
cb.set_label('PSD',fontsize=fs)
a.set_ylabel('f',fontsize=fs)
a.set_xlabel('t',fontsize=fs)
a.set_ylim(fmin/0.3,fmax/0.3)
a.set_yscale('log')
#%%
nfigs=len(figs)-1
for i in range(nfigs):
    fit(axs[i])
    hdl,lbl,lgd=getLabels(ax=axs[i],ncol=3,fontsize=15)
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)