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
var='cp' # Variable to plot
A=0 #WLE Amplitude, if SLE A=0
AoA=6 #Angle of Attack
nwave=1 #Number of WLE wavelenghts
topBlk=1 # Suction side block
botBlk=0 # Pressure side block
if A>0:
    sim='WLE'
    slices=[0,12,24,36]
else:
    sim='SLE'
    slices=np.asarray(range(0,48,1)) # Array of slice indices

nslice=len(slices) # Number of slices
#%% Paths set-up
if sim=='SLE':
    subpath='average/ss003/ss001/'
else:
    subpath='average/ss005/ss001/'
simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/dellHDD/post/{}/{}".format(user,simfolder,subpath)
spath='/home/rpt1g12/Documents/thesis/data/lowAoA/{:1d}WLE{:d}_{}Dist/'.format(nwave,AoA,var)
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
nxi=[fl.blk[botBlk].size[0],fl.blk[topBlk].size[0]]
x_bot=np.zeros((nxi[0],nslice))
x_top=np.zeros((nxi[1],nslice))
v_bot=x_bot.copy()
v_top=x_top.copy()
for kk in range(nslice):
    k=slices[kk]
    x_top[:,kk]=fl.blk[topBlk].var['x'].getValues()[:,0,k] # Extract top x coordinates
    x_bot[:,kk]=fl.blk[botBlk].var['x'].getValues()[:,-1,k] # Extract bottom x coordinates
    v_top[:,kk]=fl.blk[topBlk].var[var].getValues()[:,0,k] # Extract top variable
    v_bot[:,kk]=fl.blk[botBlk].var[var].getValues()[:,-1,k] # Extract bottom variable
#%% Plotting
fname='{}_{}'.format(simfolder,var)
f,a=getFig(fname);nfig+=1 # Get figure and axis objects
figs.append(f);axs.append(a) # Append them to the figures and axes arrays

v_top_avg=np.zeros_like(v_top[:,0])
v_bot_avg=np.zeros_like(v_top[:,0])
for kk in range(nslice):
    k=slices[kk]
    if sim=='SLE':
        v_top_avg+=v_top[:,kk]
        v_bot_avg+=v_bot[:,kk]
        if kk==(nslice-1):
            v_top_avg/=nslice
            v_bot_avg/=nslice
            if side=='top' or side=='both':
                axs[nfig].plot(x_top[:,kk],v_top_avg[:],lw=2,label='top')
            if side=='bot' or side=='both':
                axs[nfig].plot(x_bot[:,kk],v_bot_avg[:],lw=2,label='bot')
    else:
        if side=='top' or side=='both':
            axs[nfig].plot(x_top[:,kk],v_top[:,kk],lw=2,label='k{:03}_top'.format(k))
        if side=='bot' or side=='both':
            axs[nfig].plot(x_bot[:,kk],v_bot[:,kk],lw=2,label='k{:03}_bot'.format(k))

#%%
for i in range(nfig+1):    
    hdl,lbl,lgd=getLabels(axs[i],fontsize=fs,ncol=int(nslice/2))
    hdls.append(hdl)
    lbls.append(lbl)
    lgds.append(lgd)
    fit(axs[i])

#%% Saving
if save:
    savePlotFile(path=spath,ax=axs[nfig],name=fname+'.dat')
#%% If ploting Cf, add x-axis line
if var=='cf':
    axs[nfig].axhline(y=0,color='black',lw=2)
#%%
if sim=='SLE':
    v=v_top_avg
    x=x_top[:,-1]
    kk=3
else:
    kk=3
    v=v_top[:,kk]
    x=x_top[:,kk]    

vRsmp,xRsmp,nRsmp,kx=rsample(v,x,257,force=True); dx=1/kx
#%%
if var=='cp':
    dvdx=fcbFD(vRsmp,dx)
    dvMax=max(dvdx)
    dvdx/=dvMax
    dvdx2=fcbFD2(vRsmp,dx)
    dvMax2=max(dvdx2[int((nRsmp-1)/8):])
    dvdx2/=dvMax2
    #%%
    if kk==3:
        fname='{}_d{}dx'.format(simfolder,var)
        f,a=getFig(fname);nfig+=1 # Get figure and axis objects
        figs.append(f);axs.append(a) # Append them to the figures and axes arrays
        
        axs[nfig].plot(xRsmp,dvdx,lw=2,color='blue',label='d{}dx'.format(var))
        
        
        dvMaxIdx=np.argmax(dvdx)
        xm=xRsmp[dvMaxIdx]
#        axs[nfig].axvline(x=xm,lw=2,linestyle='--',color='red')
        
        i=dvMaxIdx-1;
        test0=dvdx[i+1];test=dvdx[i]
        while test>0.0:
            test0=dvdx[i];i-=1
            test=dvdx[i]    
        x0=xRsmp[i]
        
        i=dvMaxIdx+1;
        test0=dvdx[i-1];test=dvdx[i]
        while test>0.1:
            test0=dvdx[i];i+=1
            test=dvdx[i]    
        x1A=xRsmp[i]
        
        i=dvMaxIdx+1;
        test0=dvdx[i-1];test=dvdx[i]
        while test0>test:
            test0=dvdx[i];i+=1
            test=dvdx[i]    
        x1B=xRsmp[i]
        
        x1=xm#min(x1A,x1B)    
        
#        axs[nfig].axvline(x=x0,lw=2,linestyle='--',color='black')
#        axs[nfig].axvline(x=x1,lw=2,linestyle='--',color='black')
#        axs[nfig-1].axvline(x=x0,lw=2,linestyle='--',color='black')
#        axs[nfig-1].axvline(x=x1,lw=2,linestyle='--',color='black')
    
        print('Max grad: {:02.4f}'.format(dvMax))
        print('Pressure revery length: {:02.4f}'.format(x1-x0))
        print('from {:02.4f} to {:02.4f}'.format(x0,x1))
 #%%   
    pmin=min(vRsmp)
    imin=int(np.argmin(vRsmp))
    xmin=xRsmp[imin]
    print('min Cp: {:02.4f}'.format(pmin))
    print('min Cp loc: {:02.4f}'.format(xmin))

    
    axs[nfig].set_ylim(-1,1)
    axShow(axs[nfig])
#%%
else:
    i=10
    test=vRsmp[i]
    while test>0:
        i+=1
        test=vRsmp[i]
    i0=i
    x0=xRsmp[i0]

    indices=[]
    fctr=1
    for i in range(i0,nRsmp):
        if fctr*vRsmp[i]>0 and i<(nRsmp-1):
            indices.append(i)
            fctr*=(-1)
        

    i1=indices[-1]
    x1=xRsmp[i1]
        
    print('LSB length: {:02.4f}'.format(x1-x0))
    print('from {:02.4f} to {:02.4f}'.format(x0,x1))
    
    axs[nfig].axvline(x=x0,lw=2,linestyle='--',color='black')
    axs[nfig].axvline(x=x1,lw=2,linestyle='--',color='black')
    


