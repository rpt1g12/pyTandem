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
var='pp';
save=False
A=15 #WLE Amplitude, if SLE A=0
AoA=10 #Angle of Attack
block=0 #Block to look at
kk=0
vnames=['r','u','v','w','p']; # Variable names
sPattern='solT*.*.q'
#%% Paths set-up
if A>0:
    wavy=True 
    useTheta=True
    if AoA==6:
        xLSB=-0.4350
        irange=range(90,121,5)
        nw=18;ovlp=0.5
        nwave=1 #Number of LE wavelengths
        vmin=[1e-16,1e-16,1e-16,1e-12,1e-8,1e-8];
        vmax=[1e-9,1e-9,1e-9,1e-5,1e-5,1e-5]
    elif AoA==10:
        xLSB=-0.4465
        irange=range(70,101,5)
        nw=32;ovlp=0.5
        nwave=1 #Number of LE wavelengths
        vmin=[1e-12,1e-12,1e-12,1e-8,1e-7,1e-7];
        vmax=[1e-7,1e-7,1e-7,1e-5,1e-5,1e-5]
    sfolder='{}WLE'.format(nwave)
    subpath='ss004/'
    blPath='/home/rpt1g12/Documents/thesis/data/lowAoA/data/lowAoA/blAnalysis/CpLines1WLE{:02d}.dat'.format(AoA)
else:
    wavy=False   
    useTheta=True
    subpath='ss002/'
    if AoA==6:
        xLSB=-0.3594
        irange=range(100,149,8)
        nwave=1
        nw=16;ovlp=0.5

    elif AoA==10:
        nwave=1
        xLSB=-0.4180
        irange=range(90,121,5)
        nw=16;ovlp=0.5
    sfolder='{}SLE'.format(nwave)

simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/dellHDD/post/{}/unsteady/{}".format(user,simfolder,subpath)
spath='/home/rpt1g12/Documents/thesis/data/lowAoA/{}{:02d}LSBProbes/'.format(sfolder,AoA)
sfpath='/home/rpt1g12/Documents/thesis/data/lowAoA/{}{:02d}LSBProbes/'.format(sfolder,AoA)
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
#%%
flavg=fl.getAvg(files,vnames)
flavg.wrSol(sfile='solTA.qa',vnames=vnames)
up=flavg.blk[block]
#%%

varmax=flavg.blk[0].var[var].getValues().max()
varmin=flavg.blk[0].var[var].getValues().min()

f,a,im=flavg.contourf(varname=var,vmin=varmin,vmax=varmax,k=0,nlvl=21,cmap=plt.cm.hot_r);nfig+=1
figs.append(f);axs.append(a)
axs[nfig].set_aspect('equal')

#%% Read BL integral quantities
xprob,yprob=[],[]
xiprob,etprob=[],[]
if useTheta:
    hName='{}{:02d}'.format(sfolder,AoA)
    datasetBL='/home/'+user+'/Documents/thesis/data/lowAoA/blAnalysis/{}.dat'.format(hName);
    xBL,thetaBL,UeBL=np.loadtxt(datasetBL,skiprows=1,unpack=True,usecols=[0,3,6])
#    theta=rsample(thetaBL,xBL,force=True,tnew=[xLSB])[0]
#    Ue=rsample(UeBL,xBL,force=True,tnew=[xLSB])[0]
    datasetBL='/home/'+user+'/Documents/thesis/data/lowAoA/blAnalysis/CpLines{}.dat'.format(hName);
    xBL,yBL=np.loadtxt(datasetBL,skiprows=1,unpack=True,usecols=[2,3])
    ii=0;
    while not math.isnan(xBL[ii]):
        xprob.append(xBL[ii])
        yprob.append(yBL[ii])
        ii+=1
    nprob=len(xprob)
    xBK,yBK=up.var['x'].getValues()[:,:,0],up.var['y'].getValues()[:,:,0]
    nxi,net=xBK.shape
    for i in range(nprob):
        r=((xBK-xprob[i])**2+(yBK-yprob[i])**2)**0.5
        indices=np.unravel_index([np.argmin(r)],[nxi,net])
        xiprob.append(indices[0][0])
        etprob.append(indices[1][0])
    theta=1;Ue=0.3
else:
    #Probing at max RMS
    for i in irange:
        print('i={:}'.format(i))
        x,y=up.var['x'].getValues()[i,:80,kk],up.var['y'].getValues()[i,:80,kk]
        val=up.var[var].getValues()[i,:80,kk]
        ip=np.argmax(val)
        xiprob.append(i);etprob.append(ip)
        xprob.append(x[ip]);yprob.append(y[ip])
        nprob=len(xprob)
    
xprob=np.asarray(xprob)
yprob=np.asarray(yprob)
theta=1;Ue=0.3

#%% Plot and Save probe locations
probRange=np.linspace(0,nprob-1,6,dtype='int')
xNaca,yNaca=up.var['x'].getValues()[:,0,0],up.var['y'].getValues()[:,0,0]
axs[nfig].plot(xprob[probRange],yprob[probRange],lw=2,color='blue',label='probes',marker='x',ms=10)
axs[nfig].plot(xNaca[:],yNaca[:],lw=1,color='black',label='wall')
axs[nfig].set_aspect('equal')
hdl,lbl,lgd=getLabels(ax=axs[nfig])
hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
if save:
    hName='{}{:02d}'.format(sfolder,AoA)
    savePlotFile(ax=axs[nfig],path=spath,name=hName+'.dat')
axShow(axs[nfig])
#%%
ivar=['p']
nvar=len(ivar)
nprob=len(xprob)
nt=len(files)
#nt=64*60+1
t=np.zeros(nt)
data=np.zeros((nt,nprob,nvar))
for n in range(nt):
    t[n]=fl.rdFlowInfo(sfile=files[n])[3]
    fl.rdSol(sfile=files[n])
    for i in range(nvar):
        var=ivar[i]
        data[n,:,i]=fl.blk[block].var[var].getValues()[xiprob,etprob,kk]
#%% Plot and save histories
hName='{}{:02d}'.format(sfolder,AoA)

for i in range(nvar):
    var=ivar[i]
    f,a=getFig(hName+'_Hist_{}'.format(var));nfig+=1
    figs.append(f);axs.append(a)
    for j in probRange:
        axs[nfig].plot(t[:],data[:,j,i],lw=2,label='p{:03d}'.format(j))
    fit(axs[nfig])
    hdl,lbl,lgd=getLabels(ax=axs[nfig])
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)  
    if save:
        savePlotFile(ax=axs[nfig],sameX=True,path=spath)

#%%

scale=False; step=False; sclg='density'

sgn,tn,nsam,fsam=rsample(data[:,0,0],t)
nseg,novlp,ntt,fmax,fmin=defWin(tn,sgn,nw,ovlp,verbose=False)
nfreq=len(sgn)
fdata=np.zeros((int(nseg/2+1),nprob,nvar))
Data=data.copy()
for i in range(nvar):  
    for j in range(nprob):
        sgn,tn,nsam,fsam=rsample(data[:,j,i],t,verbose=True,rmAvg=True)
        Data[:,j,i]=sgn
        ff,fdata[:,j,i]=psdw(sgn,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)

#%% Plot and save frequencies histories
st=ff*theta/Ue

for i in range(nvar):
    var=ivar[i]
    f,a=getFig(hName+'_Freq_{}'.format(var));nfig+=1
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

#%%
t0=t[0]
plt.close('all')
nw=64;ovlp=0.9;sclg='density'
nseg,novlp,ntt,fmax,fmin=defWin(tn,Data[:,0,0],nw,ovlp,verbose=True)
print('fmax={}\tfmin={}'.format(fmax/0.3,fmin/0.3))
for ii,i in enumerate(probRange):
    f,a=getFig('{}Spectogram{:d}'.format(var,i));figs.append(f),axs.append(a);nfig+=1
    sgnl=Data[:,i,0]
    # Perform STFT
    ff2,tt2,psd=spectrogram(sgnl,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
    # Plot Spectogram
    TT,FF=np.meshgrid((tt2),ff2/0.3)
    im=a.pcolor(TT,FF,psd,norm=colors.LogNorm(vmin=vmin[ii], vmax=vmax[ii]),cmap=kbw,vmin=vmin[ii],vmax=vmax[ii])
    if bar and not save:
        cb=f.colorbar(im,orientation='vertical')
        cb.set_label('PSD',fontsize=fs)
    a.set_ylabel('f',fontsize=fs)
    a.set_xlabel('t',fontsize=fs)
    a.set_ylim(7,100)
    a.set_xlim((tt2[0]),(tt2[-1]))
    a.set_yscale('log')
    if save:
        figName='{}Spectogram{:d}'.format(var,i)
        saveFigOnly(path=sfpath,fig=f,ax=a,name=figName,ext='.pdf')