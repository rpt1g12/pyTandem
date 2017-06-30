# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
from lib.myPlot3dOperator import *
import lib.Plot3DClasses as p3d
from scipy import stats
from scipy.interpolate import griddata
pi=np.pi
import getpass
user=getpass.getuser()
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all');
figs=[];axs=[];hdls=[];lbls=[];lgds=[];nfig=-1;fs=18
import importlib
#%%
importlib.reload(p3d)
#%% Options
save=False
side='both'
var='cf' # Variable to plot
A=15 #WLE Amplitude, if SLE A=0
AoA=10 #Angle of Attack
nwave=1 #Number of WLE wavelenghts
topBlk=1 # Suction side block
botBlk=0 # Pressure side block
#%%
if var=='cp':
    slices=[]
    if A>0:
        sim='WLE'
        if AoA==0:
           slices+=[-0.2580, 0.1268, -0.3734, 0.0460, 0.1075]
        elif AoA==6:
           slices+=[-0.4350, -0.2388, -0.4619, -0.2849, -0.2580]
        elif AoA==10:
           slices+=[-0.4465, -0.3311, -0.4735, -0.3734, -0.3503]
    else:
        sim='SLE'
        if AoA==0:
           slices+=[-0.1406, 0.3516, -0.3555, 0.2500, 0.3047]
        elif AoA==6:
           slices+=[-0.3594, -0.0039, -0.4570, -0.0781, -0.0391]
        elif AoA==10:
           slices+=[-0.4180, -0.1328, -0.4727, -0.2109, -0.1719]
    #slices.sort()
    nslice=len(slices) # Number of slices
elif var=='cf':
    nslice=17
    if A>0:
        sim='WLE'
        if AoA==0:
           x0,x1=-0.2580, 0.1268
        elif AoA==6:
           x0,x1=-0.4350, -0.2388
        elif AoA==10:
           x0,x1=-0.4612, -0.33
    else:
        sim='SLE'
        if AoA==0:
           x0,x1=-0.1406, 0.3516 
        elif AoA==6:
           x0,x1=-0.3594, -0.0039
        elif AoA==10:
           x0,x1=-0.4180, -0.1328
    
    xm=(x0+x1)/2
    slices=np.linspace(x0,xm,nslice)
#%% Paths set-up
if sim=='SLE':
    subpath='average/ss003/ss001/'
else:
    subpath='average/ss005/ss001/'
simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/dellHDD/post/{}/{}".format(user,simfolder,subpath)
spath='/home/rpt1g12/Documents/thesis/data/lowAoA/{:1d}WLE{:d}{}spanProfile/'.format(nwave,AoA,var)
print('Reading data from:\n {}'.format(path))
#%%
vnames=['cf','twx','twy','twz','cp']; # Variable names
gfile='grid.xyz' # Grid file name
files=p3d.getFileNames(path=path,pattern='solTavgCf+tw+Cp*')
sfile=files[0]#'solTavgCf+tw+Cp.q'# Solution file name

fl=p3d.flow(path,"grid.xyz",sfile) # Check-out flow object
fl.rdHdr() # Read header from grid file and set-up flow object properties
fl.rdGrid() # Read Grid
fl.rdSol(vnames=vnames) # Read solution file
flInfo=fl.rdFlowInfo() # Read solution file info, i.e. Mach, AoA, Re and time
#%% Extract slices
nze=fl.blk[botBlk].size[2]
nxi=fl.blk[botBlk].size[0]
z=fl.blk[botBlk].var['z'].getValues()[:nslice,0,:]
x=fl.blk[topBlk].var['x'].getValues()[:,0,:]
v_top=fl.blk[topBlk].var[var].getValues()[:,0,:]
v_bot=fl.blk[botBlk].var[var].getValues()[:,0,:]
iv_top=np.zeros((nslice,nze))
iv_bot=np.zeros((nslice,nze))
for i in range(nslice):
    for k in range(nze):
        iv_top[i,k]=rsample(v_top[:,k],x[:,k],tnew=slices[i])
        iv_bot[i,k]=rsample(v_bot[:,k],x[:,k],tnew=slices[i])
#%%
fname='{}{}'.format(simfolder,var)
f,a=getFig(fname);nfig+=1 # Get figure and axis objects
figs.append(f);axs.append(a) # Append them to the figures and axes arrays

for i in range(0,nslice,4):
    axs[nfig].plot(z[i,:],iv_top[i,:],lw=2,label='i{:03d}'.format(i))

#%%
if var=='cf':
    za=np.zeros((nslice));zb=za.copy()
    for i in range(nslice):
        v=iv_top[i,:];v_rev=np.flipud(v)
        ia=np.where(v[12:]<0)[0][0]+12
        z0=float(z[i,ia]);z1=float(z[i,ia-1]);dz=z1-z0        
        v0=v[ia];v1=v[ia-1];dv=v1-v0;dvdz=dv/dz
        za[i]=z0+(-v0/dvdz);
        ib=np.where(v_rev<0)[0][0]
        if ib==0:
            ib=np.where(v>0)[0][0]
            z0=float(z[i,ib]+z[i,-1]-z[i,0]);z1=float(z[i,ib-1]+z[i,-1]-z[i,0]);dz=z1-z0        
            v0=v[ib];v1=v[ib-1];dv=v1-v0;dvdz=dv/dz
        else:
            z0=float(z[i,-(ib+1)]);z1=float(z[i,-(ib)]);dz=z1-z0        
            v0=v[-(ib+1)];v1=v[-(ib)];dv=v1-v0;dvdz=dv/dz
        
        zb[i]=z0+(-v0/dvdz);
    za=za-z[0,36]
    zb=zb-z[0,36]
    #%%
    zab=np.concatenate((np.flipud(zb),za))
    xab=np.concatenate((np.flipud(slices),slices))
    p_coeff=np.polyfit(zab,xab,2)    
    zabm=min(zab)*1.2;zabM=max(zab)*1.2
    zfit=np.linspace(zabm,zabM,51)
    p2=np.poly1d(p_coeff)
    xfit=p2(zfit)
    #%%
    fname='{}Zab'.format(simfolder)
    f,a=getFig(fname);nfig+=1 # Get figure and axis objects
    figs.append(f);axs.append(a) # Append them to the figures and axes arrays
    
    
    fp=1/(4*p_coeff[0])
    hp=-2*p_coeff[1]*fp
    kp=p_coeff[2]-hp**2/(4*fp)
    Vp=(-hp,kp);Fp=(-hp,kp+fp)
    Pp=(-hp+2*fp,kp+fp)    
    print('\nParabola coefficients:\n p(x)=(x-h)^2/(4f)+k\n (f,h,k)=({:2.4f},{:2.4f},{:2.4f})\n'.format(fp,hp,kp))
    
    axs[nfig].plot(za[:],slices,lw=2,label='za',color='blue',marker='o')
    axs[nfig].plot(zb[:],slices,lw=2,label='zb',color='red',marker='o')
    

#%%
for i in range(nfig+1):    
    hdl,lbl,lgd=getLabels(axs[i],fontsize=fs,ncol=int(nslice/2))
    hdls.append(hdl)
    lbls.append(lbl)
    lgds.append(lgd)
    fit(axs[i])


    #%% Saving
    if save:
        if i==0:
            sameX=True
        else:
            sameX=False
        savePlotFile(path=spath,ax=axs[i],sameX=sameX)
#%% If ploting Cf, add x-axis line
if var=='cf':
    axs[0].axhline(y=0,color='black',lw=2)
    axs[1].plot(zfit,xfit,lw=2,label='fit',color='black',linestyle='--')
    axs[1].scatter(Vp[0],Vp[1],c='red',s=100,label='V')
    axs[1].scatter(Fp[0],Fp[1],c='blue',s=100,label='F')
    axs[1].scatter(Pp[0],Pp[1],c='green',s=100,label='P')