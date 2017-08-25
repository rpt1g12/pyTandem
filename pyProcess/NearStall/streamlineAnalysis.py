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
play=False;playRange=range(400,600)
t0=290
ts=230
dt=1.0/64
npeak=5
maxE=8
dataset='/home/rpt1g12/Desktop/paraviewExports/P{:1d}/coordinatesExport.dat'.format(npeak)
spath='/home/rpt1g12/Desktop/paraviewExports/P{:1d}/'.format(npeak)
xyz=np.loadtxt(dataset,skiprows=1,unpack=True)
if scale:
    save=False
#%%

nt=int(xyz.shape[1]/100)
xyz=np.reshape(xyz,(3,nt,100))
z0=xyz[2,0,0]
t=np.array([t0-ts+i*(dt) for i in range(nt-1)])
print('\n z0={}\n'.format(z0))
xyz[2,:,:]=xyz[2,:,:]-z0

#%%Remove average
#for j in [2]:
#    for i in range(100):
#       xyz[j,:,i]=xyz[j,:,i]-np.mean(xyz[j,:,i])
#%%
idx=[2]
A=np.zeros((nt-1,len(idx)*97))
A2=np.zeros((nt-1,len(idx)*97))
for i in range(len(idx)):
        j=idx[i]
        A[0,i*100:(i+1)*100]=xyz[j,0,3:]
for n in range(1,nt-1):
    for i in range(len(idx)):
        j=idx[i]
        A2[n-1,i*100:(i+1)*100]=A[n,i*100:(i+1)*100]=xyz[j,n,3:]
for i in range(len(idx)):
        j=idx[i]
        A[nt-2,i*100:(i+1)*100]=xyz[j,nt-1,3:]
A1=np.transpose(A)
A2=np.transpose(A2)
#%% POD decomposition A=U.S.V*
U,s,Vt=np.linalg.svd(A1,0,1)
#%% DMD decomposition
S_1=np.diag(1/s)
V=np.transpose(Vt)
M1=np.dot(A2,V)
del A2,A
M2=np.dot(M1,S_1)
del M1,S_1
A=np.dot(U,M2)
#Eigen-decomposition of A
eV,eVec=np.linalg.eig(A)
DMD=np.dot(M2,eVec.real)
#%%Identify the unstable modes
unst_mode_index = []
eVNorm=np.zeros_like(s)
for i in range(len(s)):
    temp = np.linalg.norm(eV[i])
    eVNorm[i]=temp
    if temp > 1:
        unst_mode_index.append(i)

print ("\nNumber of unstable modes: {:d}".format(len(unst_mode_index)))
print ("Unstable modes indeces\n", unst_mode_index)
#%% Calculate Frequency and Growth rate
ntmode=len(s)
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
#%%
s0=s/sum(s)
if maxE<1:
    i=0
    while (sum(s0[:i])<maxE):
        i+=1
    nmodes=i
else:
    nmodes=maxE
    maxE=sum(s0[:nmodes])
print('\nFirst {:2d} modes contain {:2.2f}% of the energy\n'.format(nmodes,maxE*100))
f,a=getFig('S');figs.append(f),axs.append(a);nfig+=1
axs[nfig].semilogy(s0,lw=2,label='s')
#%%
#nmodes=10
f,a=getFig('time');figs.append(f),axs.append(a);nfig+=1
for i in range(nmodes):
    axs[nfig].plot(t,V[:,i]*scale*s[i],lw=2,label='m{:02d}'.format(i))
#%%
#nmodes=10
f,a=getFig('PODmode');figs.append(f),axs.append(a);nfig+=1
for i in range(nmodes):
    axs[nfig].plot(U[:,i]*scale*s[i],lw=2,label='m{:02d}'.format(i))
#%%
#nmodes=10
f,a=getFig('DMDmode');figs.append(f),axs.append(a);nfig+=1
for i in range(nmodes):
    axs[nfig].plot(DMD[:,i],lw=2,label='m{:02d}'.format(i))
#%%
nfigs=len(figs)
for i in range(nfigs):
    fit(axs[i])
    hdl,lbl,lgd=getLabels(ax=axs[i],ncol=2,fontsize=15)
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
    #Saving
    if save:
        savePlotFile(path=spath,ax=axs[i])
#%%
rModes=range(nmodes)
rA=np.zeros((97,nt-1))
for n in range(nt-1):
    for m in rModes:
        rA[:,n]+=U[:,m]*s[m]*Vt[m,n]        
#%%
if play:
    f,a=getFig('streamline');figs.append(f),axs.append(a);nfig+=1
    xmax=np.max(xyz[0,playRange,:])
    a.set_xlim(-0.5,xmax)
    zmax=np.max(np.abs(rA[:,playRange]))
    a.set_ylim(-zmax,zmax)
    for n in playRange:
        a.plot(xyz[0,n,3:],rA[:,n],lw=2,color='blue',label='rec')
        a.plot(xyz[0,n,3:],xyz[2,n,3:],lw=2,color='red',label='original')
        axShow(a)
        plt.pause(0.1)
        cll(a)
