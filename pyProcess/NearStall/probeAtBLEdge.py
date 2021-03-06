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
#%%
importlib.reload(p3d)

#%% Options
save=False
bar=True
A=15 #WLE Amplitude, if SLE A=0
AoA=10 #Angle of Attack
block=0 #Block to look at
kk=0
nwave=8
vnames=['r','u','v','w','p']; # Variable names
sPattern='solT*.*.q'
ssl=True
#%% Paths set-up
if A>0:
    nw=32;ovlp=0.5
    sfolder='{}WLE'.format(nwave)
    if ssl:
        subpath='heaving/ss005/'
        sfolder+='SSL'
    else:
        subpath='heaving/ss006/'
        sfolder+='Central'
    wavy=True
else:
    wavy=False   
    subpath='ss002/'
    sfolder='{}SLE'.format(nwave)

simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/sonyHDD/post/{}/{}".format(user,simfolder,subpath)
spath='/home/rperezt/Documents/thesis/data/nearStall/{}LSBProbes/'.format(sfolder)
sfpath='/home/rperezt/Documents/thesis/figures/nearStall/{}LSBProbes/'.format(sfolder)
blPath='/home/rperezt/Documents/thesis/data/nearStall/LSBHistory{}/BL.dat'.format(sfolder)
if not os.path.exists(spath) and save:
    os.makedirs(spath)
if not os.path.exists(sfpath) and save:
    os.makedirs(sfpath)
print('Reading data from:\n {}'.format(path))
#%%
gfile='grid.xyz' # Grid file name
files=p3d.getFileNames(path=path,pattern=sPattern)
sfile=files[0]#'solTavgCf+tw+Cp.q'# Solution file name

fl=p3d.flow(path,gfile,sfile) # Check-out flow object
fl.rdHdr() # Read header from grid file and set-up flow object properties
fl.rdGrid() # Read Grid
fl.rdSol(vnames=vnames) # Read solution file
flInfo=fl.rdFlowInfo() # Read solution file info, i.e. Mach, AoA, Re and time
mach=flInfo[0]
up=fl.blk[block]
#%% Read BL integral quantities

dataBL=np.loadtxt(blPath,skiprows=1,unpack=True)
nprob=dataBL.shape[1]
nSTA=int(dataBL.shape[0]/2)
xprob=dataBL[range(0,2*nSTA,2),:]
yprob=dataBL[range(1,2*nSTA,2),:]

xiprob=np.zeros((nSTA,nprob),dtype='int')
etprob=np.zeros((nSTA,nprob),dtype='int')
xBK,yBK=up.var['x'].getValues()[:,:,0],up.var['y'].getValues()[:,:,0]
nxi,net=xBK.shape
for n in range(nSTA):
    for i in range(nprob):
        r=((xBK-xprob[n,i])**2+(yBK-yprob[n,i])**2)**0.5
        indices=np.unravel_index([np.argmin(r)],[nxi,net])
        xiprob[n,i]=indices[0][0]
        etprob[n,i]=indices[1][0]

#%% Extract data
dtSTA=0.5
ivar=['p']
nvar=len(ivar)
nt=int(nSTA*dtSTA*64)
t=np.zeros(nt)
data=np.zeros((nt,nprob,nvar))
for N in range(nSTA):
    N0=int(N*dtSTA*64);N1=int((N+1)*dtSTA*64)
    for n in range(N0,N1):
        t[n]=fl.rdFlowInfo(sfile=files[n])[3]
        fl.rdSol(sfile=files[n])
        for i in range(nvar):
            var=ivar[i]
            data[n,:,i]=fl.blk[block].var[var].getValues()[xiprob[N,:],etprob[N,:],kk]
            
#%%
nBL=2
nplot=int(nBL*64*dtSTA*0.5)
fl.rdSol(sfile=files[nplot])
f,a,im=fl.blk[block].contourf(varname=var,vmin=0.5,vmax=0.6,nlvl=21,cmap=plt.cm.hot,bar=True);nfig+=1
figs.append(f);axs.append(a)
a.plot(xprob[nBL,:],yprob[nBL,:],lw=2,color='k',linestyle='--',marker='x')
a.set_xlim(-0.5,-0.3)
a.set_ylim(0.0,0.2)
axShow(a)
#%% Plot and save histories
probRange=range(2,7)#[0,1,2]

for i in range(nvar):
    var=ivar[i]
    f,a=getFig(sfolder+'_Hist_{}'.format(var));nfig+=1
    figs.append(f);axs.append(a)
    for j in probRange:
        axs[nfig].plot(t[:]-230,data[:,j,i],lw=2,label='p{:03d}'.format(j))
    fit(axs[nfig])
    hdl,lbl,lgd=getLabels(ax=axs[nfig])
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)  
    if save:
        savePlotFile(ax=axs[nfig],sameX=True,path=spath)

#%% PSD
plt.close('all')
scale=False; step=False;
#nw=64;ovlp=0.125;sclg='density'
nw=256;ovlp=0.5;sclg='density'
sgn,tn,nsam,fsam=rsample(data[:,0,0],t)
nseg,novlp,ntt,fmax,fmin=defWin(tn,sgn,nw,ovlp,verbose=False)
fdata=np.zeros((int(nseg/2+1),nprob,nvar))
Data=data.copy()
for i in range(nvar):  
    for j in probRange:
        sgn,tn,nsam,fsam=rsample(data[:,j,i],t,verbose=True,rmAvg=True)
        #sgn=rmvLS(tn,sgn,3)
        Data[:,j,i]=sgn
        ff,fdata[:,j,i]=psdw(sgn,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)

#%Plot and save frequencies histories
st=ff/mach

for i in range(nvar):
    var=ivar[i]
    f,a=getFig(sfolder+'_Freq_{}'.format(var));nfig+=1
    figs.append(f);axs.append(a)
    for j in probRange:
        psgn=fdata[:,j,i]
        if (step):
            spsgn=psgn*(10**j)
            if (scale):
                var=np.var(data[:,j,i])
                print(j,var,i)
                spsgn/=var
        else:
            spsgn=psgn
        axs[nfig].loglog(st,spsgn,lw=2,label='p{:03d}'.format(j))
    fit(axs[nfig])
    hdl,lbl,lgd=getLabels(ax=axs[nfig])
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
    if save:
        savePlotFile(ax=axs[nfig],sameX=True,path=spath)


#%% Spectogram
# Define windows
save=False
nw=256;ovlp=0.9;sclg='density'
nseg,novlp,ntt,fmax,fmin=defWin(tn,Data[:,0,0],nw,ovlp,verbose=False)
print('fmax={}\tfmin={}'.format(fmax/0.3,fmin/0.3))
t0=t[0]
plt.close('all')
vmin=[1e-8,1e-9,1e-9,1e-9,1e-9];vmax=[1e-5,1e-4,1e-4,1e-4,1e-4]
for ii,i in enumerate(probRange):
    f,a=getFig('{}Spectogram{:d}'.format(var,i));figs.append(f),axs.append(a);nfig+=1
    sgnl=Data[:,i,0]
    # Perform STFT
    ff2,tt2,psd=spectrogram(sgnl,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
    # Plot Spectogram
    TT,FF=np.meshgrid((tt2+t0-230),ff2/0.3)
    im=a.pcolor(TT,FF,psd,norm=colors.LogNorm(vmin=vmin[ii], vmax=vmax[ii]),cmap=kbw,vmin=vmin[ii],vmax=vmax[ii])
    if bar and not save:
        cb=f.colorbar(im,orientation='vertical')
        cb.set_label('PSD',fontsize=fs)
    a.set_ylabel('f',fontsize=fs)
    a.set_xlabel('t',fontsize=fs)
    a.set_ylim(5,100)
    a.set_xlim(46,117)
    a.set_yscale('log')
    if save:
        figName='{}Spectogram{:d}'.format(var,i)
        saveFigOnly(path=sfpath,fig=f,ax=a,name=figName,ext='.pdf')