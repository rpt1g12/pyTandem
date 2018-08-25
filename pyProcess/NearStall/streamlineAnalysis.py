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
#%%
scale=1
save=True
trange=range(000,1301)
xrange=range(0,96)
play=True;playRange=range(64)
t0=235
ts=230
dt=1.0/32
npeak=4
maxE=3
r=maxE
dataset='/home/rperezt/Desktop/paraviewExports/P{:1d}.dat'.format(npeak)
spath='/home/rperezt/Desktop/paraviewExports/P{:1d}/'.format(npeak)
xyz=np.loadtxt(dataset,skiprows=2,unpack=True)
if scale:
    save=False
    
xyz=np.loadtxt(dataset,skiprows=2,unpack=True)
fin=open(dataset)
fin.readline()
line=fin.readline()
nsteps,nt=line.rstrip('\n').rsplit('\t')[:2]
nsteps=int(nsteps);nt=int(nt)
fin.close()
#%% reshape array
xyz=np.reshape(xyz,(3,nt,nsteps))
z0=xyz[2,0,0]
print('\n z0={}\n'.format(z0))
xyz[2,:,:]=xyz[2,:,:]-z0
#%%
idx=[2]
nt=len(trange)
lx=len(xyz[0,0,xrange])
t=np.array([t0-ts+i*(dt) for i in range(nt-1)])
A=np.zeros((len(idx)*lx,nt))
for nn,n in enumerate(trange):
    temp=[]
    for i in range(len(idx)):
        j=idx[i]
        temp+=list(xyz[j,n,xrange])    
    A[:,nn]=temp


#%% Remove average
Amean=np.zeros((len(idx)*lx))
for i in range(lx):
    Amean[i,]=np.mean(A[i,:])
    A[i,:]-=Amean[i]
    #A[i,:]/=np.var(A[i,:])

Amean[:]=0
#%%
U,s,V=POD(A[:,:-1],r=-1)
# Normalise energy
sn=s/sum(s)
i=1;tot=0
while tot<0.99:
    tot=sum(sn[:i])
    i+=1
rank=i
eV,eVec,Phi,Psi,PhiPOD=DMD(A,U,s,V,mode=1,r=rank,rPro=rank)
nmodes=rank
#%%Identify the unstable modes
ntmode=Phi.shape[1]
unst_mode_index = []
eVNorm=np.zeros(ntmode)
for i in range(ntmode):
    temp = np.linalg.norm(eV[i])
    eVNorm[i]=temp
    if temp > 1:
        unst_mode_index.append(i)

print ("\nNumber of unstable modes: {:d}".format(len(unst_mode_index)))
print ("Unstable modes indeces\n", unst_mode_index)
#%% Calculate Frequency and Growth rate
freqs = np.zeros((ntmode))
grate = np.zeros((ntmode))
for i in range(ntmode):
    temp = np.log(eV[i])/dt
    freqs[i] = temp.imag/(2*pi)
    grate[i] = temp.real  

inst_max=grate.argmax()

print('Most unstable mode is:'+str(grate.argmax()))
print('Growth rate='+str(grate[inst_max]))
print('Freq='+str(freqs[inst_max]))

#%% Get least damped modes
leastDamp=[]
temp=list(grate)
minG=grate.min()
for i in range(nmodes):
    ii=np.argmax(temp)
    leastDamp.append(ii)
    temp[ii]=minG
print('\nLeast damped modes are:\n{}\nwith growth rate:\n{}\nand frequencies:{}\n'.format(leastDamp,grate[leastDamp],freqs[leastDamp]))
#%%

f,a=getFig('PODtime');figs.append(f),axs.append(a);nfig+=1
for i in range(nmodes):
    axs[nfig].plot(t,V[:,i]*scale*s[i],lw=2,label='m{:02d}'.format(i))
#%%
#nmodes=10
f,a=getFig('PODmode');figs.append(f),axs.append(a);nfig+=1
for i in range(nmodes):
    axs[nfig].plot(U[:,i]*scale*s[i],lw=2,label='m{:02d}'.format(i))
#%%
if r<0:
    r=len(s)
f,a=getFig('DMDmode');figs.append(f),axs.append(a);nfig+=1

plotRange=leastDamp
for i in plotRange:
    axs[nfig].plot(Phi.real[:,i],lw=2,label='m{:02d}R'.format(i))
#%%
#nmodes=10
f,a=getFig('DMDtime');figs.append(f),axs.append(a);nfig+=1
for i in plotRange:
    axs[nfig].plot(t,Psi.real[i,:],lw=2,label='m{:02d}R'.format(i))
#%% DMD spectrum
f,a=getFig('DMDspectrum');figs.append(f),axs.append(a);nfig+=1
a.scatter(freqs,grate,s=PhiPOD*100,c=eVNorm,cmap=plt.cm.jet)
a.axhline(y=0,lw=2,linestyle='--')
fmin=min(freqs);fmax=max(freqs);fabs=max(fmax,np.abs(fmin))
a.set_xlim(-fabs,fabs)    
#%%
nfigs=len(figs)-1
for i in range(nfigs):
    fit(axs[i])
    hdl,lbl,lgd=getLabels(ax=axs[i],ncol=2,fontsize=15)
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
    #Saving
    if save:
        savePlotFile(path=spath,ax=axs[i])
#%%
rModes=range(nmodes)
Sr=np.diag(s)[:nmodes,:nmodes]
Ur=U[:,:nmodes]
Vr=V.conj().T[:nmodes,:]
rPOD=np.dot(Ur,np.dot(Sr,Vr))  
#rDMD=np.dot(Phi[:,leastDamp],Psi[leastDamp,:])    
rDMD=np.dot(Phi,Psi).real   
#%%
if play:
    f,a=getFig('streamline');figs.append(f),axs.append(a);nfig+=1
    xmax=np.max(xyz[0,playRange,:])
    #a.set_xlim(-0.5,xmax)
    zmax=np.max(np.abs(rPOD[:,playRange]))
    #a.set_ylim(-zmax,zmax)
    for n in playRange:
        a.plot(A[:,n]+Amean[:],lw=2,color='red',label='original')
        a.plot(rDMD[:,n]+Amean[:],lw=2,color='green',linestyle=':',label='DMD')
        a.plot(rPOD[:,n]+Amean[:],lw=2,color='blue',linestyle='--',label='POD')
        axShow(a)
        plt.pause(0.1)
        cll(a)
