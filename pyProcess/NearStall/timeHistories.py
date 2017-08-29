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
readSol=False
save=True
bar=False
cutoff=7
filt=False
detrend=True
dtrndString='_original'
#%% Options
sPattern='solT*.q'
vnames=['r','u','v','w','p']; # Variable names
vname='w'
#%% Simulation Options
A=15 #WLE Amplitude, if SLE A=0
AoA=10 #Angle of Attack
iRange=np.asarray(range(30,91,10))
kRange=np.asarray(range(6,6+24*7+1,24))
nwave=8 #Number of LE wavelengths
#%% Paths set-up
if A>0:
    wavy=True
    sfolder='{}WLE'.format(nwave)
    subpath='heaving/ss001/'
    block=1 #Block to look at
else:
    wavy=False
    sfolder='{}SLE'.format(nwave)
    subpath='heaving/ss003/'
    block=4 #Block to look at


simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/dellHDD/post/{}/{}".format(user,simfolder,subpath)
spath='/home/rpt1g12/Documents/thesis/data/nearStall/{}{}History/'.format(sfolder,vname)
sfpath='/home/rpt1g12/Documents/thesis/figures/nearStall/{}{}History/'.format(sfolder,vname)

dataPath=spath+'data.dat'
tPath=spath+'t.dat'

if not os.path.exists(spath):
    os.makedirs(spath)
if not os.path.exists(sfpath):
    os.makedirs(sfpath)
print('Reading data from:\n {}'.format(path))
#%%
gfile='grid.xyz' # Grid file name
files=p3d.getFileNames(path=path,pattern=sPattern);len(files)
sfile=files[0]#'solTavgCf+tw+Cp.q'# Solution file name

fl=p3d.flow(path,gfile,sfile) # Check-out flow object
fl.rdHdr() # Read header from grid file and set-up flow object properties
fl.rdGrid() # Read Grid
fl.rdSol(vnames=vnames) # Read solution file
flInfo=fl.rdFlowInfo() # Read solution file info, i.e. Mach, AoA, Re and time
mach=flInfo[0]

nx,nz=len(iRange),len(kRange)

#%% Read data from Plot3d files if readSol=True
if readSol:
    data=np.zeros((nx,nz,nt))
    t=np.zeros(nt)
    i=0
    for sfile in files:
        fl.rdSol(vnames=vnames,sfile=sfile)
        flInfo=fl.rdFlowInfo()
        for ii in range(nx):
            iii=iRange[ii]
            data[ii,:,i]=fl.blk[block].var[vname].getValues()[iii,1,kRange]
        t[i]=flInfo[3]
        i+=1
    
    data.tofile(dataPath)
    t.tofile(tPath)
else:
    t=np.fromfile(tPath)
    nt=len(t)
    data=np.reshape(np.fromfile(dataPath),(nx,nz,nt))

#%% Resample
deltaT=t[-2]-t[0]
fsam=64
t0=t[0]
t=t-t0

tnew=np.linspace(0,deltaT,deltaT*fsam+1)
nsam=len(tnew)
#%%
n=np.where(t>deltaT)[0][0]
Data=np.zeros((nx,nz,nsam))
Data_filt=Data.copy()
for i in range(nx):
    for k in range(nz):
        Data[i,k,:],tt,nsam,fsam=rsample(data[i,k,:n],t[:n],verbose=False)
        if filt:        
            Data_filt[i,k,:]=fourierFilter(Data[i,k,:],'high',fsam,cutoff*0.3,15,False)
        if detrend:
            Data[i,k,:]=rmvLS(tt,Data[i,k,:],2)
            dtrndString='_detrend'
#%%
tt0=tt+t0-230
xplotRange=[0,1]
zplotRange=range(nz)#[4,5]#,3,4,5]
plt.close('all')
for i in xplotRange:
    f,a=getFig('{}History{:d}{}'.format(vname,i,dtrndString));figs.append(f),axs.append(a);nfig+=1
    for k in zplotRange:
        a.plot(tt0,Data[i,k,:],lw=2,label='T{:1d}'.format(k+1))
    hdl,lbl,lgd=getLabels(ax=a,ncol=2,fontsize=15)
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
    if save:
        savePlotFile(path=spath,ax=a,sameX=True)
#%%
if filt:
    for i in range(nx):
        f,a=getFig('{}HistoryFilt{:d}{}'.format(vname,i,dtrndString));figs.append(f),axs.append(a);nfig+=1
        for k in range(nz):
            a.plot(tt0,Data_filt[i,k,:],lw=2,label='T{:1d}'.format(k+1))
    hdl,lbl,lgd=getLabels(ax=a,ncol=2,fontsize=15)
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
    if save:
        savePlotFile(path=spath,ax=a,sameX=True)
#%% Define windows
vmin=1e-9;vmax=1e-6
nw=16;ovlp=0.9;sclg='density'
nw=64;ovlp=0.9;sclg='density'
nseg,novlp,ntt,fmax,fmin=defWin(tt,Data[0,0,:],nw,ovlp,verbose=False)
print('fmax={}\tfmin={}'.format(fmax/0.3,fmin/0.3))
#%% PSD
plt.close('all')
for i in xplotRange:
    f,a=getFig('{}PSD{:d}'.format(vname,i));figs.append(f),axs.append(a);nfig+=1
    for k in zplotRange:
        if filt:
            sgnl=Data_filt[i,k,:]
        else:
            sgnl=Data[i,k,:]       
        # Perform FFT
        ff3,PP=psdw(sgnl,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
        # Plot frequencies      
        a.loglog(ff3/0.3,PP,label='T{:1d}'.format(k+1))
    a.set_xlim(1,100)
    hdl,lbl,lgd=getLabels(ax=a,ncol=3,fontsize=15)
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
    if save:
        savePlotFile(path=spath,ax=a,sameX=True)
#%% Spectogram
plt.close('all')
for i in xplotRange:
    for k in zplotRange:
        f,a=getFig('{}Spectogram{:d}_{:d}'.format(vname,i,k+1));figs.append(f),axs.append(a);nfig+=1
        if filt:
            sgnl=Data_filt[i,k,:]
        else:
            sgnl=Data[i,k,:]
        # Perform STFT
        ff2,tt2,psd=spectrogram(sgnl,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)#,window='boxcar')
        # Plot Spectogram
        TT,FF=np.meshgrid((tt2+t0-230),ff2/0.3)
        im=a.pcolor(TT,FF,psd,norm=colors.LogNorm(vmin=vmin, vmax=vmax),cmap=kbw,vmin=vmin,vmax=vmax)
        if bar and not save:
            cb=f.colorbar(im,orientation='vertical')
            cb.set_label('PSD',fontsize=fs)
        a.set_ylabel('f',fontsize=fs)
        a.set_xlabel('t',fontsize=fs)
        a.set_ylim(7,100)
        a.set_xlim((tt2[0]+t0-230),(tt2[-1]+t0-230))
        a.set_yscale('log')
        if save:
            figName='{}Spectogram{:d}_{:d}'.format(vname,i,k+1)
            saveFigOnly(path=sfpath,fig=f,ax=a,name=figName,ext='.pdf')

