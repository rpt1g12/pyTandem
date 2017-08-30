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
read=True
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
xiRange=np.asarray(range(50,109))
etRange=np.asarray(range(1,30))
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
path="/media/{}/dellHDD/post/{}/{}".format(user,simfolder,subpath)
spath='/home/rpt1g12/Documents/thesis/figures/nearStall/ROM{}/'.format(sfolder)

dpath='/home/rpt1g12/Documents/thesis/data/nearStall/ROM{}/'.format(sfolder)
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
flROM=fl.getSubsets(xRanges=[xiRange],yRanges=[etRange],zRanges=[range(1)],ssName='ROM',sspath=path+'ROM/')
temp=flROM.blk[block].var['x'].getValues()
xiD=len(xiRange)
etD=len(etRange)
D=np.zeros((xiD*etD,nt))
t=np.zeros(nt)
#%%
if read:
    for n,file in enumerate(files):
        print('{:3.2f}%'.format(n/nt))
        fl.rdSol(vnames=vnames,sfile=file)
        t[n]=fl.rdFlowInfo(file)[3]
        for i,ii in enumerate(xiRange):
            for j,jj in enumerate(etRange):
                D[i+j*xiD,n]=fl.blk[block].var[vname].getValues()[ii,jj,kk]
    
    D.tofile(dpath+'Dkk{:d}.dat'.format(kk))
    t.tofile(dpath+'tkk{:d}.dat'.format(kk))
else:
    print('Reading previously stored data...')
    t=np.fromfile(dpath+'tkk{:d}.dat'.format(kk))
    nt=len(t)
    D=np.reshape(np.fromfile(dpath+'Dkk{:d}.dat'.format(kk)),(xiD*etD,nt))
D0=D.copy()   
#%% Remove average
Dmean=np.zeros(xiD*etD)
Dvar=Dmean.copy()
Dnorm=D.copy()
for i in range(xiD*etD):
    Dmean[i]=np.mean(D[i,:])
    Dvar[i]=np.var(D[i,:])
    D[i,:]-=Dmean[i]
    Dnorm[i,:]=D[i,:]/Dvar[i]
#%%
nmodes=20
plt.close('all')
tstep=1
tRange=range(250,nt,tstep)
U,s,V,eV,eVec,Phi,Psi,PhiPOD=ROM(D[:,tRange],r=-1)
#%Identify the unstable modes
ntmode=U.shape[1]
unst_mode_index = []
eVNorm=np.zeros(ntmode)
for i in range(ntmode):
    temp = np.linalg.norm(eV[i])
    eVNorm[i]=temp
    if temp > 1:
        unst_mode_index.append(i)


print ("\nUnstable modes indeces\n", unst_mode_index)
print ("\nNumber of unstable modes: {:d} out of {:d} total".format(len(unst_mode_index),ntmode))

#% Calculate Frequency and Growth rate
tn=t[tRange][:-1]-230
dt=tstep/64.0
freqs = np.zeros((ntmode))
grate = np.zeros((ntmode))
for i in range(ntmode):
    temp = np.log(eV[i])/(dt*mach)
    freqs[i] = temp.imag/(2*pi)
    grate[i] = temp.real  

inst_max=grate.argmax()

print('Most unstable mode is:'+str(grate.argmax()))
print('Growth rate='+str(grate[inst_max]))
print('Freq='+str(freqs[inst_max]))

#% Get least damped modes
leastDamp=[]
temp=list(grate)
minG=grate.min()
for i in range(nmodes):
    ii=np.argmax(temp)
    leastDamp.append(ii)
    temp[ii]=minG
print('\nLeast damped modes are:\n{}'.format(leastDamp)+
      '\nwith growth rate:\n{}\n'.format(grate[leastDamp])+
      'and frequencies:{}\n'.format(freqs[leastDamp]))

#%% Normalise energy
sn=s/sum(s)
# Plot POD energy
f,a=getFig('Sigma')
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
axs[nfig].plot(sn,lw=2)      
#%% Plot DMD coefficients for least-damped modes
if len(leastDamp)<len(unst_mode_index) or len(unst_mode_index)==0:
    plotRange=leastDamp
else:
    plotRange=unst_mode_index
f,a=getFig('DMDtime');figs.append(f),axs.append(a);nfig+=1
for i in plotRange:
    axs[nfig].plot(tn,Psi.real[i,:],lw=2,label='m{:02d}R'.format(i))
hdl,lbl,lgd=getLabels(ax=a,ncol=3,fontsize=15)
hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
fit(a)
#%%
#plt.close('all')
for nn,nmode in enumerate(plotRange):
##    #Plot POD
#    temp=np.reshape(U[:,nn],(xiD,etD,1),order='F')
#    flROM.blk[block].setData(vname='POD{:d}'.format(nn),val=temp/np.max(temp))
#    vmax=1
#    vmin=-vmax
#    f,a,im=flROM.blk[block].contourf(varname='POD{:d}'.format(nn),vmin=vmin,vmax=vmax,nlvl=21,cmap=plt.cm.jet,bar=True);nfig+=1
#    figs.append(f);axs.append(a)
 
   #Plot DMD
    temp=np.reshape(Phi[:,nmode],(xiD,etD,1),order='F')
    flROM.blk[block].setData(vname='DMD{:d}'.format(nmode),val=temp/np.max(temp))
    f,a,im=flROM.blk[block].contourf(varname='DMD{:d}'.format(nmode),vmin=-1,vmax=1,nlvl=21,cmap=plt.cm.jet,bar=True);nfig+=1
    figs.append(f);axs.append(a)
    #flROM.blk[block].contour(varname='DMD{:d}'.format(nmode),vmin=-1,vmax=1,nlvl=11,ax=a)
    ftag='f={:3.2f}\ng={:3.2f}'.format(freqs[nmode],grate[nmode])  
    a.text(0.1,0.75,ftag,ha='center',va='bottom',transform=a.transAxes)    

# DMD spectrum
f,a=getFig('DMDspectrum');figs.append(f),axs.append(a);nfig+=1
a.scatter(freqs,grate,c=PhiPOD*100,s=eVNorm*100/eVNorm.max(),cmap=plt.cm.jet)
a.axhline(y=0,lw=2,linestyle='--')
fabs=np.max(np.abs(freqs))
a.set_xlim(-fabs,fabs)  
  
#%%
print ("\nUnstable modes indeces\n", unst_mode_index)
print ("\nNumber of unstable modes: {:d} out of {:d} total".format(len(unst_mode_index),ntmode))  

print('Most unstable mode is:'+str(grate.argmax()))
print('Growth rate='+str(grate[inst_max]))
print('Freq='+str(freqs[inst_max]))

print('\nLeast damped modes are:\n{}'.format(leastDamp)+
      '\nwith growth rate:\n{}\n'.format(grate[leastDamp])+
      'and frequencies:{}\n'.format(freqs[leastDamp]))