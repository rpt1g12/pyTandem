# -*- coding: utf-8 -*-
import os
import math
import numpy as np
import scipy
from scipy.signal import spectrogram 
from scipy.signal import welch as psdw
from scipy.signal import blackman
from scipy.signal import correlate
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
bar=True
cutoff=2
filt=False
detrend=True
rmAvg=True
dtrndString='_original'
media='sonyHDD'
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
path="/media/{}/{}/post/{}/{}".format(user,media,simfolder,subpath)
sdpath='/home/{}/Documents/thesis/data/nearStall/{}{}History/'.format(user,sfolder,vname)
spath='/home/{}/Documents/thesis/data/nearStall/{}{}Phase/'.format(user,sfolder,vname)
sfpath='/home/{}/Documents/thesis/figures/nearStall/{}{}History/'.format(user,sfolder,vname)

dataPath=sdpath+'data.dat'
tPath=sdpath+'t.dat'

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
    
t0=t[0]
t=t-t0
#%% Resample
deltaT=t[-1]-t[0]
fsam=64


tnew=np.linspace(0,deltaT,deltaT*fsam)
nsam=len(tnew)
#%%
n=np.where(t<deltaT)[0][-1]+1
Data=np.zeros((nx,nz,nsam))
Data_filt=Data.copy()
detrend=False
for i in range(nx):
    for k in range(nz):
        Data[i,k,:],tt,nsam,fsam=rsample(data[i,k,:n],t[:n],verbose=False,rmAvg=False)

        if filt:        
            #Data_filt[i,k,:]=fourierFilter(Data[i,k,:],'low',fsam,cutoff*0.3,15,False)
            Data[i,k,:]=smooth(Data[i,k,:],int(64/(cutoff*0.3)))
        if detrend:
            Data[i,k,:]=rmvLS(tt,Data[i,k,:],2)
            dtrndString='_detrend'

tt=t[:n]
#%% Fourier decomposition
plt.close('all')
DD=Data
pks=[3,4]
axs=[];figs=[]
nTot=len(DD[0,0,:])
ntop=int(nTot/2)

phase=np.zeros((2,ntop))
tRange=range(ntop,nTot)

if tRange[0]==0:
    suffix='_prestall'
else:
    suffix='_poststall'
    
prefix='P{:d}P{:d}_'.format(pks[0]+1,pks[1]+1)

ff=np.linspace(0,fsam/2,ntop)

#y0=np.sin(2*pi*0.25*tt)
y0=DD[0,pks[0],tRange]#*blackman(ntop)
#y1=3*np.sin(2*pi*0.25*tt+0.25*pi)
y1=DD[0,pks[1],tRange]#*blackman(ntop)

f,a=getFig(prefix+'signals'+suffix);figs.append(f);axs.append(a);nfig+=1
a.plot(tt[tRange]+60,y0,lw=2,color='blue',label='P{:d}'.format(pks[0]+1))
a.plot(tt[tRange]+60,y1,lw=2,color='red',label='P{:d}'.format(pks[1]+1))
hdl,lbl,lgd=getLabels(ax=a,ncol=2,fontsize=15)
hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
if save:
    savePlotFile(path=spath,ax=a,sameX=True)

nw=30;ovlp=0.0;sclg='density'
nseg,novlp,ntt,fmax,fmin=defWin(tt[:ntop],Data[0,0,:ntop],nw,ovlp,verbose=False)
print('fmax={}\tfmin={}'.format(fmax,fmin))


fff,csd=scipy.signal.csd(y0,y1,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
st=fff/0.3

f,a=getFig('csd');figs.append(f);axs.append(a);nfig+=1
a.loglog(st,np.abs(csd))

f,a=getFig('angle');figs.append(f);axs.append(a);nfig+=1
angle=np.arctan2(csd.imag,csd.real)
a.plot(st,angle/(pi))

f,a=getFig(prefix+'cos(angle)'+suffix);figs.append(f);axs.append(a);nfig+=1
a.plot(st,np.cos(angle),lw=2,label='cos')
hdl,lbl,lgd=getLabels(ax=a,ncol=2,fontsize=15)
hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
if save:
    savePlotFile(path=spath,ax=a,sameX=True)