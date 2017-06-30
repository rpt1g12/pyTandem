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
A=15 #WLE Amplitude, if SLE A=0
AoA=0 #Angle of Attack
block=4 #Block to look at
krange=range(1,49) #Spanwise slice to look at
nwave=1 #Number of LE wavelengths
#%% Rake set-up
ibounds=[10,271] #Manually selected bounds for rakes to be placed inside LSB
nn=1024; #Points in rake
l=0.13; #Rake's length
step=1 #Iteration step
auto=True #Automatic LSB index detection
#%% xLSB set-up
if A>0:
    sfolder='{}WLE'.format(nwave)
    subpath='average/ss006/'
    if AoA==0:
        x0,x1=-0.2580, 0.1268
    elif AoA==6:
       x0,x1=-0.4350, -0.2388
    elif AoA==10:
       x0,x1=-0.4612, -0.33
else:
    sfolder='{}SLE'.format(nwave)
    subpath='average/ss004/'
    if AoA==0:
       x0,x1=-0.1406, 0.3516 
    elif AoA==6:
       x0,x1=-0.3594, -0.0039
    elif AoA==10:
       x0,x1=-0.4180, -0.1328       
xLSB=(x0+x1)/2       
#%% Paths set-up
simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/dellHDD/post/{}/{}".format(user,simfolder,subpath)
spath='/home/rpt1g12/Documents/thesis/data/lowAoA/deltaHalfLSB/'
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
up.setData('Wz',Wz)

#%%
zLSB=[];deltaLSB=[];zOUT=[];deltaOUT=[]
for kk in krange:
    
    #%% Obtain surface normals
    nx,ny,xo,yo=(up.mets[3].getValues()[:,0,kk],up.mets[4].getValues()[:,0,kk],
                up.var['x'].getValues()[:,0,kk],up.var['y'].getValues()[:,0,kk])
    
    wz_wall=Wz[:,0,kk]
    #%% Get indices inside the LSB
    rev_idx=np.where(wz_wall>0)[0]
    xLSB_idx=np.where(xo>xLSB)[0][0]
    lsb_idx=np.asarray(range(xLSB_idx-2,xLSB_idx+3))
    
    if len(rev_idx)>0 and rev_idx[0]<xLSB_idx:
        inout=True
        zLSB.append(up.var['z'].getValues()[0,0,kk])
        if zLSB[-1]>0.0275:
            zLSB[-1]=zLSB[-1]-0.11
    else:
        inout=False
        zOUT.append(up.var['z'].getValues()[0,0,kk])
        
   
    print('from i={:03d} to {:03d}'.format(lsb_idx[0],lsb_idx[-1]))
    #%% Set-up rakes
    jdelta=[];delta=[];x0=[] #Emply lists to store: deltaStar, theta, ...
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
        print('i={:}, k={:}'.format(i,kk))
        rake=(p3d.rake(xo[i],yo[i],nx[i],ny[i],nn,l)) #Checkout a rake object at (xo,yo) with direction (nx,ny) and length l
        x,y=rake.x,rake.y  #Extract xy coordinates of rake
        tx,ty=rake.tx,rake.ty #Extract tangent direction
        wz=up.interpolate2dk('Wz',x,y,kk,mode='values')[:,0,0] #Obtain interpolated values of omega_z at rake's points
        u=-myIntegral(wz,rl);u[np.isnan(u)]=max(u) #Define a pseudo velocity based on omega_z
        dudy=fcbFD(u,dl) #Compute pseudo velocity derivative in the wall-normal direction
        j99=np.where(u<0.99*u[-1])[0][-1]+1;jmax=-1 #Find out index at which delta99 happens
        jdelta.append(j99)
        delta.append(rl[j99]) #Save delta, i.e. boundary layer thickness
        x0.append(x[0])
        ii+=1
    delta=np.asarray(delta)
    deltaK=rsample(delta,x0,force=True,tnew=[xLSB])[0]
    if inout:
        deltaLSB.append(deltaK)
    else:
        deltaOUT.append(deltaK)
        
#%% Join lists and sort
dummy=list(deltaOUT)
deltaOUT=[x for (y,x) in sorted(zip(zOUT,dummy))]
zOUT=sorted(zOUT)
dummy=list(deltaLSB)
deltaLSB=[x for (y,x) in sorted(zip(zLSB,dummy))]
zLSB=sorted(zLSB)
zALL=list(zOUT)+list(zLSB)
dummy=list(deltaOUT)+list(deltaLSB)
deltaALL=[x for (y,x) in sorted(zip(zALL,dummy))]
zALL=sorted(zALL)
#%%Lists to numpy arrays
zLSB=np.asarray(zLSB)
deltaLSB=np.asarray(deltaLSB)
zOUT=np.asarray(zOUT)
deltaOUT=np.asarray(deltaOUT)
#%% Parabolic fit
p_coeff=np.polyfit(zLSB,deltaLSB,2)    


fp=1/(4*p_coeff[0])
hp=-2*p_coeff[1]*fp
kp=p_coeff[2]-hp**2/(4*fp)
Vp=(-hp,kp);Fp=(-hp,kp+fp)
Pp=(-hp+2*fp,kp+fp)    
print('\nParabola coefficients:\n p(x)=(x-h)^2/(4f)+k\n (f,h,k)=({:2.4f},{:2.4f},{:2.4f})\n'.format(fp,hp,kp))

width=max([hp-min(zLSB),max(zLSB)-hp])
zLSBm=hp-1.2*width;zLSBM=hp+1.2*width
zfit=np.linspace(zLSBm,zLSBM,51)
p2=np.poly1d(p_coeff)
deltaFit=p2(zfit)
#%% Plot delta and delta_star
f,a=getFig('{}{:02d}'.format(sfolder,AoA))
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
axs[nfig].plot(zALL,deltaALL,lw=2,label='delta')
#axs[nfig].scatter(zLSB,deltaLSB,lw=2,label='LSB',c='blue',s=150)
#axs[nfig].scatter(zOUT,deltaOUT,lw=2,label='OUT',c='red',s=150)
axs[nfig].plot(zfit,deltaFit,lw=2,linestyle='--',color='black',label='Fit')
#%%
for i in range(nfig+1):
    hdl,lbl,lgd=getLabels(ax=axs[i])
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
    fit(axs[i])
    if save:
        savePlotFile(path=spath,ax=axs[i])
        