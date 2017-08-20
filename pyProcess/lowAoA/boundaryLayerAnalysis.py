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
sPattern='solTA*'
vnames=['r','u','v','w','p']; # Variable names
#%% Options
save=True
comp=True; #Compressible??
A=0 #WLE Amplitude, if SLE A=0
AoA=10 #Angle of Attack
block=4 #Block to look at
kk=24 #Spanwise slice to look at
nwave=1 #Number of LE wavelengths
#%% Rake set-up
ibounds=[10,271] #Manually selected bounds for rakes to be placed inside LSB
nn=1024; #Points in rake
l=0.13; #Rake's length
step=1 #Iteration step
auto=True #Automatic LSB index detection
#%% Paths set-up
if A>0:
    wavy=True
    sfolder='{}WLE'.format(nwave)
    subpath='average/ss006/'
else:
    wavy=False
    sfolder='{}SLE'.format(nwave)
    subpath='average/ss004/'
simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/dellHDD/post/{}/{}".format(user,simfolder,subpath)
spath='/home/rpt1g12/Documents/thesis/data/lowAoA/blAnalysis/'
if not os.path.exists(spath) and save:
    os.makedirs(spath)
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
cosa,sina=np.cos(np.deg2rad(flInfo[1])),np.sin(np.deg2rad(flInfo[1]))
up=fl.blk[block]
del fl
#%% Compute spanwise vorticity
up.getMetrics()
Wz=up.derive('v','x',True)-up.derive('u','y',True)
U=np.sqrt(up.var['u'].getValues()**2+up.var['v'].getValues()**2)
up.setData('Wz',Wz)
up.setData('U',U)
Cp=((up.var['p'].getValues()-1/1.4)/(0.5*mach**2))
up.setData('Cp',Cp)

#%% Obtain surface normals
nx,ny,xo,yo=(up.mets[3].getValues()[:,0,kk],up.mets[4].getValues()[:,0,kk],
            up.var['x'].getValues()[:,0,kk],up.var['y'].getValues()[:,0,kk])

wz_wall=Wz[:,0,kk]
#%% Get indices inside the LSB
rev_idx=np.where(wz_wall>0)[0]
if auto: 
    if not wavy:
        rev_idx_diff=rev_idx[1:]-rev_idx[0:-1]
        bigdiff=np.where(rev_idx_diff>10)[0]
        if len(bigdiff)>0:
            rev_idx[-1]=rev_idx[bigdiff[0]-1]
    lsb_idx=np.asarray(range(rev_idx[0]-10,rev_idx[-1]+11))
else:
    lsb_idx=np.asarray(range(ibounds[0],ibounds[1])) # Ignore automatic LSB index detection

print('from i={:03d} to {:03d}'.format(lsb_idx[0],lsb_idx[-1]))
#%% Plot contour of omega_z
if AoA==0:
    vmin=-1;vmax=1
elif AoA==6:
    vmin=-2;vmax=1
elif AoA==10:
    vmin=-3.5;vmax=1
f,a,im=up.contourf(varname='Cp',vmin=vmin,vmax=vmax,k=kk,nlvl=21,cmap=plt.cm.hot_r,bar=False)
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
up.contour(varname='Cp',vmin=vmin,vmax=vmax,k=kk,nlvl=21,ax=axs[nfig])

a.set_aspect('equal')
axs[nfig].set_xlim(-0.5,0.5)
axs[nfig].set_ylim(0.0,0.25)
if save:
    saveFigOnly(path=spath,fig=figs[nfig],ax=axs[nfig],name='Cp{}{:02d}'.format(sfolder,AoA),ext='.pdf')

cll(axs[nfig])
#%% Set-up rakes
dstr=[];theta=[];jdelta=[];xbl=[];ybl=[];delta=[];Ue=[];x0=[] #Emply lists to store: deltaStar, theta, ...
dl=l/(nn-1)
rl=np.array([i*dl for i in range(nn)])
u=np.zeros_like(rl)
dudy=np.zeros(nn)
nrakes=len(lsb_idx[::step])
npoint=len(rl)
u99=np.zeros((nrakes,npoint))
#%% Loop inside LSB indices
ii=0
for i in lsb_idx[::step]:
    print('i={:}'.format(i))
    rake=(p3d.rake(xo[i],yo[i],nx[i],ny[i],nn,l)) #Checkout a rake object at (xo,yo) with direction (nx,ny) and length l
    x,y=rake.x,rake.y  #Extract xy coordinates of rake
    tx,ty=rake.tx,rake.ty #Extract tangent direction
    wz=up.interpolate2dk('Wz',x,y,kk,mode='data')[:,0,0] #Obtain interpolated values of omega_z at rake's points
    u=-myIntegral(wz,rl);u[np.isnan(u)]=max(u) #Define a pseudo velocity based on omega_z
    r=up.interpolate2dk('r',x,y,kk,mode='data')[:,0,0] #Obtain density at rake's points
    dudy=fcbFD(u,dl) #Compute pseudo velocity derivative in the wall-normal direction
    j99=np.where(u<0.99*u[-1])[0][-1]+1;jmax=-1 #Find out index at which delta99 happens
    jdelta.append(j99)
    delta.append(rl[j99]) #Save delta, i.e. boundary layer thickness
    um=u[j99] #Save boundary layer edge velocity Ue
    Ue.append(um)
    x0.append(x[0])
    u99[ii,:]=u/um #Normalise velocity by Ue
    xbl.append(x[j99])
    ybl.append(y[j99])
    if comp:        
        rm=r[j99]
        r/=rm
        tmp=1-(r[:j99])*u[:j99]
        tmp2=rl[:j99]
        dstr.append(np.trapz(tmp,tmp2))
        tmp=((r[:j99])*u[:j99])*(1-u[:j99])
        theta.append(np.trapz(tmp,tmp2))  
    else:
        tmp=1-u[:j99]
        tmp2=rl[:j99]
        dstr.append(np.trapz(tmp,tmp2))
        tmp=u[:j99]*(1-u[:j99])
        theta.append(np.trapz(tmp,tmp2))
    ii+=1

    
#%%Lists to numpy arrays
xbl=np.asarray(xbl) 
ybl=np.asarray(ybl) 
x0=np.asarray(x0)    
Ue=np.asarray(Ue)    
theta=np.asarray(theta)
dstr=np.asarray(dstr)
delta=np.asarray(delta)
H=dstr/theta

#%%Interpolate at delta

Udelta=up.interpolate2dk('U',xbl,ybl,kk,mode='data')[:,0,0]
CpDelta=up.interpolate2dk('Cp',xbl,ybl,kk,mode='data')[:,0,0]

#%%  Plot boundary layer edge
figs[0].canvas.set_window_title('CpLines{}{:02d}'.format(sfolder,AoA))
axs[nfig].plot(xo,yo,lw=2,color='black',label='wall')
axs[0].plot(xbl,ybl,lw=2,color='black',linestyle='--',label='BL')
a.set_aspect('equal')
axs[nfig].set_xlim(-0.5,0.5)
axs[nfig].set_ylim(0.0,0.25)
#%% Plot delta and delta_star
f,a=getFig('deltas')
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
axs[nfig].plot(x0,delta,lw=2,label='delta')
axs[nfig].plot(x0,dstr,lw=2,label='delta_star')

#%% Plot momentum thickness theta
f,a=getFig('theta')
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
axs[nfig].plot(x0,theta,lw=2,label='theta')

#%% Plot shape factor H
f,a=getFig('ShapeFactor') 
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
axs[nfig].plot(x0,H,lw=2,label='H')

#%% Plot edge velocity
f,a=getFig('Ue') 
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
axs[nfig].plot(x0,Ue,lw=2,label='ue')
axs[nfig].plot(x0,Udelta,lw=2,label='Udelta')
#%% Plot edge pressure
f,a=getFig('Cp') 
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
axs[nfig].plot(x0,CpDelta,lw=2,label='Cp')
#%% Plot edge velocity
f,a=getFig('{}{:02d}'.format(sfolder,AoA))
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
axs[nfig].plot(x0,delta,lw=2,label='delta')
axs[nfig].plot(x0,dstr,lw=2,label='delta_star')
axs[nfig].plot(x0,theta,lw=2,label='theta')
axs[nfig].plot(x0,H,lw=2,label='H')
axs[nfig].plot(x0,Ue,lw=2,label='ue')
axs[nfig].plot(x0,Udelta,lw=2,label='uDelta')
axs[nfig].plot(x0,CpDelta,lw=2,label='CpDelta')
#%%
i=0
hdl,lbl,lgd=getLabels(ax=axs[i])
hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
fit(axs[i])
if save:
    savePlotFile(path=spath,ax=axs[i],sameX=False)
for i in range(nfig,nfig+1):
    hdl,lbl,lgd=getLabels(ax=axs[i])
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
    fit(axs[i])
    if save:
        savePlotFile(path=spath,ax=axs[i],sameX=True)
        