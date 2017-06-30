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
styles=['-','--',':','-.']
#%%
importlib.reload(p3d)

#%% Options
save=False
auto=False #Automatic LSB index detection
ibounds=[11,101] #Manually selected bounds for rakes to be placed inside LSB
comp=True; #Compressible??
A=15 #WLE Amplitude, if SLE A=0
AoA=20 #Angle of Attack
block=4 #Block to look at
kk=94 #Spanwise slice to look at
nwave=8 #Number of LE wavelengths
nn=1024; #Points in rake
l=0.13; #Rake's length
step=5 #Iteration step
#%% Paths set-up
if A>0:
    sfolder='{}WLE'.format(nwave)
else:
    sfolder='{}SLE'.format(nwave)
subpath='average/vsmallDomain/'
simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/dellHDD/post/{}/{}".format(user,simfolder,subpath)
spath='/home/rpt1g12/Documents/thesis/data/StalledWLE/{}_LSBProfiles/'.format(sfolder)
if not os.path.exists(spath) and save:
    os.makedirs(spath)
print('Reading data from:\n {}'.format(path))
#%%
vnames=['r','u','v','w','p']; # Variable names
gfile='grid.xyz' # Grid file name
sfile='solTA_shifted.qa'# Solution file name

fl=p3d.flow(path,gfile,sfile) # Check-out flow object
fl.rdHdr() # Read header from grid file and set-up flow object properties
fl.rdGrid() # Read Grid
fl.rdSol(vnames=vnames) # Read solution file
flInfo=fl.rdFlowInfo() # Read solution file info, i.e. Mach, AoA, Re and time
mach=flInfo[0]
cosa,sina=np.cos(np.deg2rad(flInfo[1])),np.sin(np.deg2rad(flInfo[1]))
up=fl.blk[block]
del fl
#%% Obtain surface normals
up.getMetrics()

nx,ny,xo,yo=(up.mets[3].getValues()[:,0,kk],up.mets[4].getValues()[:,0,kk],
            up.var['x'].getValues()[:,0,kk],up.var['y'].getValues()[:,0,kk])
#%% Compute spanwise vorticity
Wz=up.derive('v','x',True)-up.derive('u','y',True)
up.setData('Wz',Wz)
wz_wall=Wz[:,0,kk]
#%% Get indices inside the LSB
rev_idx=np.where(wz_wall>0)[0]
if auto:
    lsb_idx=[];ii=0
    for i in rev_idx:
        lsb_idx.append(i)
        test=(lsb_idx[ii]-lsb_idx[ii-1])>1
        if test:
            del lsb_idx[-1]
            break
        ii+=1
else:
    lsb_idx=np.asarray(range(ibounds[0],ibounds[1])) # Ignore automatic LSB index detection
#%% Plot contour of omega_z
#f,a,im=up.contourf(varname='p',vmin=0.65,vmax=0.7,k=kk,nlvl=50,cmap=plt.cm.hot,bar=False)
f,a,im=up.contourf(varname='Wz',vmin=-10,vmax=10,k=kk,nlvl=50,cmap=plt.cm.bwr,bar=False)
a.set_aspect('equal')
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
axs[0].set_xlim(-0.6,-0.1)
axs[0].set_ylim(0.0,0.5)
saveFigOnly(path=spath,fig=figs[0],ax=axs[0],name='omega_z')
#%% Check out figures and axes
f,a=getFig('v_Profiles') #Check out a figure for velocity profiles
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
f,a=getFig('cp_Profiles') #Check out a figure for cp profiles
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
#f,a=getFig('Ue') #Check out a figure Ue along x
#figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
#%% Set-up rakes
rakes=[]; #Create empty array to store rake objects
dstr=[];jdelta=[];xbl=[];ybl=[];delta=[];Ue=[] #Emply lists to store: deltaStar, theta, ...
dl=l/(nn-1)
rl=np.array([i*dl for i in range(nn)])
u=np.zeros_like(rl)
dudy=np.zeros(nn)
nrakes=len(lsb_idx[::step])
npoint=len(rl)
u99=np.zeros((nrakes,npoint))
cp99=u99.copy()
#%% Loop inside LSB indices
ii=0
for i in lsb_idx[::step]:
    print('i={:}'.format(i))
    rakes.append(p3d.rake(xo[i],yo[i],nx[i],ny[i],nn,l)) #Checkout a rake object at (xo,yo) with direction (nx,ny) and length l
    x,y=rakes[-1].x,rakes[-1].y  #Extract xy coordinates of rake
    tx,ty=rakes[-1].tx,rakes[-1].ty #Extract tangent direction
    wz=up.interpolate2dk('Wz',x,y,kk,mode='values')[:,0,0] #Obtain interpolated values of omega_z at rake's points
    u=-myIntegral(wz,rl);u[np.isnan(u)]=max(u) #Define a pseudo velocity based on omega_z
    cp=up.interpolate2dk('p',x,y,kk,mode='values')[:,0,0] #Obtain pressure at rake's points
    cp99[ii,:]=2*(cp-1/1.4)/mach**2
    rakes[-1].var=u #Set rake's values to pseudo velocity
    dudy=fcbFD(rakes[-1].var,dl) #Compute pseudo velocity derivative in the wall-normal direction
    j99=np.where(u<0.99*u[-1])[0][-1]+1;jmax=-1 #Find out index at which delta99 happens
    jdelta.append(j99)
    delta.append(rl[j99]) #Save delta, i.e. boundary layer thickness
    um=u[j99] #Save boundary layer edge velocity Ue
    Ue.append(um)
    rakes[-1].var/=um #Normalise velocity by Ue
    u99[ii,:]=rakes[-1].var[:]
    style=styles[int(ii/7)]
    axs[1].plot(u99[ii,:],rl[:]/delta[ii],lw=2,linestyle=style,label='r{:02d}'.format(ii))
    axs[0].plot(x,y,lw=2,linestyle=style,label='r{:02d}'.format(ii)) 
    xbl.append(rakes[-1].x[j99])
    ybl.append(rakes[-1].y[j99])
    ii+=1
#%%  Plot boundary layer edge
axs[0].plot(xbl,ybl,lw=2,color='black',linestyle='--',label='BL')
axs[1].set_ylim(0,1)
fit(axs[0],spg=(0.1,0.2))
##%% Plot boundary layer edge velocity
#axs[3].plot(xbl,Ue,lw=2,color='blue',label='Ue')
#fit(axs[3])
#%% Cp average
cp_avg=0
for i in range(nrakes):
    cp_avg+=np.mean(cp99[i,:jdelta[i]])
cp_avg/=nrakes
for i in range(nrakes):
    style=styles[int(i/7)]
    axs[2].plot((cp99[i,:]-cp_avg)/np.abs(cp_avg),rl[:]/delta[i],lw=2,linestyle=style,label='r{:02d}'.format(i))
axs[2].set_ylim(0,1)
axs[2].set_xlim(-0.5,0.5)
#%%
if save:
    for i in range(nfig+1):
        hdl,lbl,lgd=getLabels(ax=axs[i])
        hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
        savePlotFile(path=spath,ax=axs[i])