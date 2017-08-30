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
nst=64
nt=64*12
sPattern='solT*.q'
vnames=['r','u','v','w','p']; # Variable names
#%% Options
save=True
comp=True; #Compressible??
A=15 #WLE Amplitude, if SLE A=0
AoA=10 #Angle of Attack
block=0 #Block to look at
kk=0 #Spanwise slice to look at
nwave=8 #Number of LE wavelengths
ssl=True
#%% Rake set-up
ibounds=[40,81] #Manually selected bounds for rakes to be placed inside LSB
nn=1024; #Points in rake
l=0.12; #Rake's length
step=5 #Iteration step
auto=True #Automatic LSB index detection
#%% Paths set-up
if A>0:
    wavy=True
    sfolder='{}WLE'.format(nwave)
    if ssl:
        subpath='heaving/ss005/STA/'
        sfolder+='SSL'
    else:
        subpath='heaving/ss006/STA/'
        sfolder+='Central'
else:
    wavy=False
    sfolder='{}SLE'.format(nwave)
    subpath='heaving/ss001/'
simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/dellHDD/post/{}/{}".format(user,simfolder,subpath)
spath='/home/rpt1g12/Documents/thesis/data/nearStall/LSBHistory{}/'.format(sfolder)
if not os.path.exists(spath) and save:
    os.makedirs(spath)
print('Reading data from:\n {}'.format(path))
#%% Prepare STA flow object

gfile='grid.xyz' # Grid file name
files=p3d.getFileNames(path=path,pattern=sPattern);nt=len(files)


fl=p3d.flow(path,gfile,files[0]) # Check-out flow object
fl.rdHdr() # Read header from grid file and set-up flow object properties
fl.rdGrid() # Read Grid
fl.rdSol(vnames=vnames) # Read solution file
flInfo=fl.rdFlowInfo() # Read solution file info, i.e. Mach, AoA, Re and time

#%% Obtain surface normals
fl.blk[block].getMetrics()
nx,ny,xo,yo=(fl.blk[block].mets[3].getValues()[:,0,kk],fl.blk[block].mets[4].getValues()[:,0,kk],
            fl.blk[block].var['x'].getValues()[:,0,kk],fl.blk[block].var['y'].getValues()[:,0,kk])

#%% Get indices inside the LSB

lsb_idx=np.asarray(range(ibounds[0],ibounds[1])) # Ignore automatic LSB index detection

print('from i={:03d} to {:03d}'.format(lsb_idx[0],lsb_idx[-1]))
#%% Set-up rakes
dl=l/(nn-1)
rl=np.array([i*dl for i in range(nn)])
u=np.zeros_like(rl)
dudy=np.zeros(nn)
nrakes=len(lsb_idx[::step])
npoint=len(rl)
u99=np.zeros((nrakes,npoint))
#%% Loop inside LSB indices and create rakes
ii=0
rakes=[]
for i in lsb_idx[::step]:
    print('i={:}'.format(i))
    rake=(p3d.rake(xo[i],yo[i],nx[i],ny[i],nn,l)) #Checkout a rake object at (xo,yo) with direction (nx,ny) and length l
    rakes.append(rake)
    ii+=1
nrakes=len(rakes)
#%% Prepare Arrays
xbl=np.zeros((nt,nrakes))
ybl=xbl.copy()
x0=xbl.copy()
Ue=xbl.copy()
theta=xbl.copy()
dstr=xbl.copy()
delta=xbl.copy()
H=xbl.copy()
#%%
for n in range(nt):
    fl.rdSol(vnames=vnames,sfile=files[n]) # Read solution file  
    Wz=fl.blk[block].derive('v','x',True)-fl.blk[block].derive('u','y',True)
    fl.blk[block].setData('Wz',Wz)      
    # Loop inside LSB indices
    for ii,rake in enumerate(rakes):
        print('rake={:}'.format(ii))
        x,y=rake.x,rake.y  #Extract xy coordinates of rake
        tx,ty=rake.tx,rake.ty #Extract tangent direction
        wz=fl.blk[block].interpolate2dk('Wz',x,y,kk,mode='data')[:,0,0] #Obtain interpolated values of omega_z at rake's points
        u=-myIntegral(wz,rl);u[np.isnan(u)]=max(u) #Define a pseudo velocity based on omega_z
        r=fl.blk[block].interpolate2dk('r',x,y,kk,mode='data')[:,0,0] #Obtain density at rake's points
        dudy=fcbFD(u,dl) #Compute pseudo velocity derivative in the wall-normal direction
        j99=np.where(u<0.99*u[-1])[0][-1]+1;jmax=-1 #Find out index at which delta99 happens
        delta[n,ii]=rl[j99] #Save delta, i.e. boundary layer thickness
        um=u[j99] #Save boundary layer edge velocity Ue
        Ue[n,ii]=um
        x0[n,ii]=x[0]
        u99[ii,:]=u/um #Normalise velocity by Ue
        xbl[n,ii]=x[j99]
        ybl[n,ii]=y[j99]
        if comp:        
            rm=r[j99]
            r/=rm
            tmp=1-(r[:j99])*u[:j99]
            tmp2=rl[:j99]
            dstr[n,ii]=np.trapz(tmp,tmp2)
            tmp=((r[:j99])*u[:j99])*(1-u[:j99])
            theta[n,ii]=np.trapz(tmp,tmp2)
        else:
            tmp=1-u[:j99]
            tmp2=rl[:j99]
            dstr[n,ii]=np.trapz(tmp,tmp2)
            tmp=u[:j99]*(1-u[:j99])
            theta[n,ii]=np.trapz(tmp,tmp2)

    H[n,:]=dstr[n,:]/theta[n,:]
    n+=1

#%% Plot delta and delta_star
f,a=getFig('BLquantities')
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
for n in range(nt):
    axs[nfig].plot(x0[n,:],delta[n,:],lw=2,label='d{:02d}'.format(n))
    axs[nfig].plot(x0[n,:],dstr[n,:],lw=2,label='dstr{:02d}'.format(n))
    axs[nfig].plot(x0[n,:],theta[n,:],lw=2,label='theta{:02d}'.format(n))
    axs[nfig].plot(x0[n,:],H[n,:],lw=2,label='H{:02d}'.format(n))

i=nfig
hdl,lbl,lgd=getLabels(ax=axs[i])
hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
fit(axs[i])
if save:
    savePlotFile(path=spath,ax=axs[i],sameX=False)
#%% edge
f,a=getFig('BL')
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays

for n in range(nt):
    axs[nfig].plot(xbl[n,:],ybl[n,:],lw=2,label='BL{:02d}'.format(n))
axs[nfig].set_xlim(-0.5,-0.3)
axs[nfig].set_aspect('equal')

i=nfig
hdl,lbl,lgd=getLabels(ax=axs[i])
hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
fit(axs[i])
if save:
    savePlotFile(path=spath,ax=axs[i],sameX=False)

#%
axs[nfig].plot(xo,yo,lw=2,color='k',label='wall')