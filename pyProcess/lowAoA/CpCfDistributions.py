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
topBlk=1 # Suction side block
botBlk=0 # Pressure side block
slices=np.asarray(range(0,48,12)) # Array of slice indices
nslice=len(slices) # Number of slices
#%% Paths set-up
subpath='average/ss005/ss001/'
simfolder='1A{:02d}W11AoA{:02d}'.format(A,AoA)
path="/media/{}/dell\'s_hdd/post/{}/{}".format(user,simfolder,subpath)
spath='/home/rpt1g12/Documents/thesis/data/lowAoA/{}Distribution/'.format(var)
print('Reading data from:\n {}'.format(path))
#%%
vnames=['cf','twx','twy','twz','cp']; # Variable names
gfile='grid.xyz' # Grid file name
sfile='solTavgCf+tw+Cp.q'# Solution file name

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
    x_bot[:,kk]=fl.blk[botBlk].var['x'].getValues()[:,0,k] # Extract bottom x coordinates
    v_top[:,kk]=fl.blk[topBlk].var[var].getValues()[:,0,k] # Extract top variable
    v_bot[:,kk]=fl.blk[botBlk].var[var].getValues()[:,0,k] # Extract bottom variable
#%% Plotting
fname='{}_{}'.format(simfolder,var)
f,a=getFig(fname);nfig+=1 # Get figure and axis objects
figs.append(f);axs.append(a) # Append them to the figures and axes arrays

for kk in range(nslice):
    k=slices[kk]
    if side=='top' or side=='both':
        axs[nfig].plot(x_top[:,kk],v_top[:,kk],lw=2,label='k{:03}_top'.format(k))
    if side=='bot' or side=='both':
        axs[nfig].plot(x_bot[:,kk],v_bot[:,kk],lw=2,label='k{:03}_bot'.format(k))

for i in range(nfig+1):    
    hdl,lbl,lgd=getLabels(axs[i],fontsize=fs,ncol=int(nslice/2))
    hdls.append(hdl)
    lbls.append(lbl)
    lgds.append(lgd)
    fit(axs[i])

#%% Saving
if save:
    savePlotFile(path=spath,ax=axs[nfig],name=fname)
#%% If ploting Cf, add x-axis line
if var=='cf':
    axs[nfig].axhline(y=0,color='black',lw=2)

