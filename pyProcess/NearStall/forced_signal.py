# -*- coding: utf-8 -*-
import os
import math
import numpy as np
from scipy.signal import spectrogram 
from scipy.signal import welch as psdw
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
save=False
bar=True
cutoff=2
filt=False
detrend=False
rmAvg=True
dtrndString='_original'
media='082F-63FE'
#%% Options
sPattern='solT*.q'
vnames=['r','u','v','w','p']; # Variable names
vname='w'
#%% Simulation Options
A=15 #WLE Amplitude, if SLE A=0
AoA=17 #Angle of Attack
iRange=np.asarray(range(30,91,10))
kRange=np.asarray(range(6,6+24*7+1,24))
nwave=8 #Number of LE wavelengths
#%% Paths set-up
if A>0:
    wavy=True
    sfolder='{}WLE'.format(nwave)
    subpath='ss001/'
    block=1 #Block to look at
else:
    wavy=False
    sfolder='{}SLE'.format(nwave)
    subpath='heaving/ss003/'
    block=4 #Block to look at


simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/{}/post/{}/{}".format(user,media,simfolder,subpath)
spath='/media/{}/{}/phd/thesis/data/nearStall/{}{}History/'.format(user,media,sfolder,vname)
sfpath='/media/{}/{}/phd/thesis/figures/nearStall/{}{}History/'.format(user,media,sfolder,vname)

dataPath=spath+'data.dat'
tPath=spath+'t.dat'

if not os.path.exists(spath):
    os.makedirs(spath)
if not os.path.exists(sfpath):
    os.makedirs(sfpath)
print('Reading data from:\n {}'.format(path))
#%%
gfile='grid.xyz' # Grid file name
files=p3d.getFileNames(path=path,pattern=sPattern);nt = len(files)
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
for i in range(nx):
    for k in range(nz):
        Data[i,k,:],tt,nsam,fsam=rsample(data[i,k,:n],t[:n],verbose=False,rmAvg=rmAvg,tnew=tnew)
        if filt:        
            #Data_filt[i,k,:]=fourierFilter(Data[i,k,:],'low',fsam,cutoff*0.3,15,False)
            Data_filt[i,k,:]=smooth(Data[i,k,:],int(64/(cutoff*0.3)))
        if detrend:
            Data[i,k,:]=rmvLS(tt,Data[i,k,:],2)
            dtrndString='_detrend'

#%%
plt.close('all')
for k in [0,1,2,3,4,5,6,7]:
    plt.plot(tnew,Data[0,k,:],lw=2,label="peak {}".format(k))
handle,labels,legend = getLabels()
#%% Signal Correlation
i=0
tt0 = tnew
ks=[1,2,3,4]
mode='full'
plt.close('all')
lag=[]
DD=Data
wLength=10*64
nShift=4*64
f,a=getFig('Histories12');figs.append(f),axs.append(a);nfig+=1
a.plot(tt0,DD[0,ks[0],:],lw=2,color='blue')
a.plot(tt0,DD[0,ks[1],:],lw=2,color='red')
f,a=getFig('Histories34');figs.append(f),axs.append(a);nfig+=1
a.plot(tt0,DD[0,ks[2],:],lw=2,color='blue')
a.plot(tt0,DD[0,ks[3],:],lw=2,color='red')