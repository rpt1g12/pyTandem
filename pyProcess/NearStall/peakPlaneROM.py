# -*- coding: utf-8 -*-
import os
import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
from lib.myPlot3dOperator import *
import lib.Plot3DClasses as p3d
from scipy import stats
from scipy.interpolate import griddata
from scipy.signal import argrelextrema as locmaxmin
pi=np.pi
import getpass
user=getpass.getuser()
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')
figs=[];axs=[];hdls=[];lbls=[];lgds=[];nfig=-1;fs=18
import importlib
#%%
importlib.reload(p3d)
#%% Options
read=False
save=False
sPattern='solT*.q'
vnames=['r','u','v','w','p']; # Variable names
plane=2
autoMinMax=False
bar=True
#%% Simulation Options
vname='w'

A=15 #WLE Amplitude, if SLE A=0
AoA=10 #Angle of Attack
block=0 #Block to look at
kk=6 #Spanwise slice to look at
nwave=8 #Number of LE wavelengths
#xiRange=np.asarray(range(50,109))
#etRange=np.asarray(range(1,30))
xiRange=np.asarray(range(40,80))
etRange=np.asarray(range(0,30))
zeRange=np.asarray(range(kk,kk+1))
#%% Paths set-up
if A>0:
    wavy=True
    sfolder='{}WLEk{:d}'.format(nwave,kk)
    subpath='heaving/ss003/'
else:
    wavy=False
    sfolder='{}SLEk{:d}'.format(nwave,kk)
    subpath='heaving/ss003/'

simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/sonyHDD/post/{}/{}".format(user,simfolder,subpath)
spath='/home/rperezt/Documents/thesis/figures/nearStall/ROM{}/'.format(sfolder)

dpath='/home/rperezt/Documents/thesis/data/nearStall/ROM{}/'.format(sfolder)
if not os.path.exists(spath) and save:
    os.makedirs(spath)
if not os.path.exists(dpath) and save:
    os.makedirs(dpath)
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
nt=len(files)
nxi,net,nze=fl.blk[block].size
#%%
flROM=fl.getSubsets(xRanges=[xiRange],yRanges=[etRange],zRanges=[zeRange],ssName='ROM',sspath=path+'ROM/')
temp=flROM.blk[block].var['x'].getValues()
xiD=len(xiRange)
etD=len(etRange)
xiD0=110
etD0=49
D0=np.zeros((xiD0*etD0,nt))
D=np.zeros((xiD*etD,nt))
t=np.zeros(nt)
#%%
if read:
    for n,file in enumerate(files):
        print('{:3.2f}%'.format(100*n/nt))
        fl.rdSol(vnames=vnames,sfile=file)
        t[n]=fl.rdFlowInfo(file)[3]
        for i,ii in enumerate(xiRange):
            for j,jj in enumerate(etRange):
                D0[i+j*xiD,n]=fl.blk[block].var[vname].getValues()[ii,jj,kk]
    
    D0.tofile(dpath+'Dkk{:d}.dat'.format(kk))
    t.tofile(dpath+'tkk{:d}.dat'.format(kk))
else:
    print('Reading previously stored data...')
    t=np.fromfile(dpath+'tkk{:d}.dat'.format(kk))
    nt=len(t)
    D0=np.reshape(np.fromfile(dpath+'Dkk{:d}.dat'.format(kk)),(xiD0*etD0,nt))
    for n,file in enumerate(files):
        temp=np.reshape(D0[:,n],(etD0,xiD0))
        for i,ii in enumerate(xiRange):
                for j,jj in enumerate(etRange):
                    D[i+j*xiD,n]=temp[jj-1,ii]
    
  
#%% Remove average
Dmean=np.zeros(xiD*etD)
Dvar=Dmean.copy()
Dnorm=D.copy()
for i in range(xiD*etD):
    Dmean[i]=np.mean(D[i,:])
    Dvar[i]=(np.var(D[i,:]))**0.5
    D[i,:]-=Dmean[i]
    Dnorm[i,:]=D[i,:]/Dvar[i]

#%% Define frequency bins
fbinLim=np.linspace(-0.5,99.5,101)
fbin=np.zeros(101)
flen=len(fbin)
#%%

#for niter,t0 in enumerate(range(0,601,50)):
for niter,t0 in enumerate([400]):
    nmodes=20
    plt.close('all')
    print('iteration: {}'.format(niter))
    tstep=1
    # you get nice modes with 250<t0<800 and t1=2369 (t'=80)
    #t1=2369
    t1=2000#np.where(t-230>78)[0][0]
    print('using from t0={:d} to t1={:d}'.format(t0,t1))
    tRange=range(t0,t1,tstep)
    U,s,V,eV,eVec,Phi,Psi,PhiPOD=ROM(D[:,tRange],mode=2,r=400,npro=nmodes)
    #%Identify the unstable modes
    ntmode=U.shape[1]
    unst_mode_index = []
    eVNorm=np.zeros(ntmode)
    for i in range(ntmode):
        temp = np.linalg.norm(eV[i])
        eVNorm[i]=temp
    
    #% Calculate Frequency and Growth rate
    tn=t[tRange][:-1]-230
    dt=tstep/64.0
    freqs = np.zeros((ntmode))
    grate = np.zeros((ntmode))
    for i in range(ntmode):
        temp = np.log(eV[i])/(dt*mach)
        freqs[i] = temp.imag/(2*pi)
        grate[i] = temp.real
        if grate[i]>=0 and freqs[i] >= 0:
            unst_mode_index.append(i)
    
    inst_max=grate.argmax()
    
    print ("\nUnstable modes indeces\n", unst_mode_index)
    totUnst=len(unst_mode_index)
    totPosFreq=len(np.asarray(np.where(freqs>=0))[0])
    print ("\nNumber of unstable modes: {:d} out of {:d} total".format(totUnst,totPosFreq))  
    print("{:2.2f}%".format(100*totUnst/totPosFreq))
    
    print('Most unstable mode is:'+str(grate.argmax()))
    print('Growth rate='+str(grate[inst_max]))
    print('Freq='+str(freqs[inst_max]))
    
    #% Get least damped modes
    leastDamp=[]
    temp=list(grate)
    minG=grate.min()
    while len(leastDamp)<nmodes:
        ii=np.argmax(temp)
        if freqs[ii]>=0:
            leastDamp.append(ii)
        temp[ii]=minG
    print('\nLeast damped modes are:\nmode\tG\tf')
    for i in leastDamp:
        print('{:03d}\t{:+02.4f}\t{:02.4f}'.format(i,grate[i],freqs[i]))
    #%Fill frequency bins
    for i in unst_mode_index:
        freq=freqs[i]
        if freq>=0 and freq<=100.5:
            ibin=int(np.round(freq))
            fbin[ibin]+=1
    for i in range(flen):
        fbin[i]=min(niter+1,fbin[i])
#%% Plot unstable frequency PDF
f,a=getFig('fBins')
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
#axs[nfig].plot(fbin/len(range(250,601,25)),lw=2,label='fPDF')
#if save:
#    savePlotFile(ax=axs[nfig],sameX=True,path=dpath)
axs[nfig].bar(fbinLim,fbin/(niter+1),width=1,color='blue')
axs[nfig].set_xlim(-0.5,50)
#axs[nfig].set_ylim(0,1)
axs[nfig].set_xlabel(r'$f$',fontsize=28)
axs[nfig].set_ylabel(r'$Probability$',fontsize=28)
axs[nfig].set_xticks(range(102))
plt.show()
#%% Normalise energy
sn=s/sum(s)
# Plot POD energy
f,a=getFig('Sigma')
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
axs[nfig].loglog(sn,lw=2) 
sSum=s.copy() 
for i in range(len(s)):
    sSum[i]=sum(sn[:i+1])
f,a=getFig('SigmaSum')
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
axs[nfig].loglog(sSum,lw=2)
axs[nfig].axhline(y=0.9,lw=2,color='k')

# Plot DMD coefficients for least-damped modes
if len(leastDamp)<len(unst_mode_index) or len(unst_mode_index)==0:
    plotRange=leastDamp
else:
    plotRange=unst_mode_index

#plotRange=np.where(freqs>0)[0]

f,a=getFig('DMDtime');figs.append(f),axs.append(a);nfig+=1
for i in plotRange:
    axs[nfig].plot(tn,Psi.real[i,:],lw=2,label='m{:02d}R'.format(i))
hdl,lbl,lgd=getLabels(ax=a,ncol=3,fontsize=15)
hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
fit(a)

 #%% Plot POD time coefficients
plt.close('all')
nrange=len(plotRange)
f,a=getFig('PODtime');figs.append(f),axs.append(a);nfig+=1
for i in range(5):
    axs[nfig].plot(tn,s[i]*V[:,i],lw=2,label='m{:02d}'.format(i))
hdl,lbl,lgd=getLabels(ax=a,ncol=3,fontsize=15)
hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
fit(a)

#%% Plot POD
s90=min(nmodes,np.where(sSum>0.9)[0][0])
for nn,nmode in enumerate(range(s90)):
    temp=np.reshape(U[:,nn],(xiD,etD,1),order='F')
    flROM.blk[block].setData(vname='POD{:d}'.format(nn),val=temp/np.max(temp))
    vmax=1
    vmin=-vmax
    f,a,im=flROM.blk[block].contourf(varname='POD{:d}'.format(nn),vmin=vmin,vmax=vmax,nlvl=21,cmap=plt.cm.jet,bar=True);nfig+=1
    figs.append(f);axs.append(a)
#%% Plot DMD
#plt.close('all')
#plotRange=range(20)
#plotRange=[164,51] 
#plotRange=np.where(PhiPOD>0.95)[0]   
for nn,nmode in enumerate(plotRange):

    temp=np.reshape(Phi[:,nmode],(xiD,etD,1),order='F')
    flROM.blk[block].setData(vname='DMD{:d}'.format(nmode),val=temp/np.max(temp))
    f,a,im=flROM.blk[block].contourf(varname='DMD{:d}'.format(nmode),vmin=-1,vmax=1,nlvl=21,cmap=plt.cm.jet,bar=True);nfig+=1
    figs.append(f);axs.append(a)
    #flROM.blk[block].contour(varname='DMD{:d}'.format(nmode),vmin=-1,vmax=1,nlvl=11,ax=a)
    ftag='f={:3.2f}\ng={:3.2f}'.format(freqs[nmode],grate[nmode])  
    a.text(0.1,0.75,ftag,ha='center',va='bottom',transform=a.transAxes,fontsize=28)    

# DMD spectrum
f,a=getFig('DMDspectrum');figs.append(f),axs.append(a);nfig+=1
a.scatter(freqs,grate,s=PhiPOD*100,c=eVNorm,cmap=plt.cm.jet)
a.axhline(y=0,lw=2,linestyle='--')
fabs=np.max(np.abs(freqs))
a.set_xlim(-fabs,fabs)  
  
#%%
print ("\nUnstable modes indeces\n", unst_mode_index)
print ("\nNumber of unstable modes: {:d} out of {:d} total".format(len(unst_mode_index),ntmode))  

print('Most unstable mode is:'+str(grate.argmax()))
print('Growth rate='+str(grate[inst_max]))
print('Freq='+str(freqs[inst_max]))

print('\nLeast damped modes are:\nmode\tG\tf')
for i in leastDamp:
    print('{:03d}\t{:+02.4f}\t{:02.4f}'.format(i,grate[i],freqs[i]))