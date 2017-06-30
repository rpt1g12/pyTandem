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
save=True
prandtl=False
side='both'
var='cp' # Variable to plot
M_inf=0.3 #Mach number
A=0 #WLE Amplitude, if SLE A=0
AoA=6 #Angle of Attack
nwave=1 #Number of WLE wavelenghts
topBlk=1 # Suction side block
botBlk=0 # Pressure side block
if A>0:
    sim='WLE'
    slices=np.asarray(range(0,49,1)) # Array of slice indices
else:
    sim='SLE'
    slices=np.asarray(range(0,49,1)) # Array of slice indices

nslice=len(slices) # Number of slices

#%%Scale coefficients with prandtl rule
if prandtl:
    beta=np.sqrt(1-M_inf**2)
else:
    beta=1

pinf=1/1.4
pdyn=0.5*M_inf**2
#%% Paths set-up
if sim=='SLE':
    subpath='average/ss003/ss001/'
else:
    subpath='average/ss005/ss001/'
simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/dellHDD/post/{}/{}".format(user,simfolder,subpath)
spath='/home/rpt1g12/Documents/thesis/data/lowAoA/{:1d}WLE{:d}_{}maxCp/'.format(nwave,AoA,var)
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
#%%Get maximum pressure
maxCp=np.zeros(len(slices))
maxLoc=maxCp.copy()
for k in slices:
    maxCp[k]=max(v_bot[:,k])
    imax=np.argmax(v_bot[:,k])    
    maxLoc[k]=x_bot[imax,k]

#%% Plotting
fname='{}maxCp'.format(simfolder)
f,a=getFig(fname);nfig+=1 # Get figure and axis objects
figs.append(f);axs.append(a) # Append them to the figures and axes arrays    
Lz=0.11*nwave
z=np.flipud(np.linspace(-0.5*Lz,0.5*Lz,nslice))

#axs[nfig].plot(z,beta*(maxCp*pdyn+pinf),lw=2,color='blue',label='maxp')
axs[nfig].plot(z,maxCp,lw=2,color='blue',label='maxCp')

fname='{}maxLoc'.format(simfolder)
f,a=getFig(fname);nfig+=1 # Get figure and axis objects
figs.append(f);axs.append(a) # Append them to the figures and axes arrays    

axs[nfig].plot(z,maxLoc,lw=2,color='blue',label='maxLoc')

#%%
for i in range(nfig+1):    
    hdl,lbl,lgd=getLabels(axs[i],fontsize=fs,ncol=1)
    hdls.append(hdl)
    lbls.append(lbl)
    lgds.append(lgd)
    fit(axs[i])
    #%% Saving
    if save:
        savePlotFile(path=spath,ax=axs[i])