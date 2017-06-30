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
save=False
A=15 #WLE Amplitude, if SLE A=0
AoA=6 #Angle of Attack
nwave=1 #Number of WLE wavelenghts

if A>0:
    sim='WLE'
else:
    sim='SLE'

#%%
importlib.reload(p3d)
#%% Options
sPattern='solTA*'
vnames=['r','u','v','w','p']; # Variable names
#%% Sub-set settings
ssName='aerofoil'
bkXYZ=[1,2,1]
blocks=[1,4]
xRanges=[None]
yRanges=[[-1],[0]]
zRanges=[None]
    
#%% Paths set-up
if sim=='SLE':
    subpath='average/ss003/'
else:
    subpath='average/fullDomain/ss001/'
simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/dellHDD/post/{}/{}".format(user,simfolder,subpath)
sspath=path+ssName+'/'
print('Reading data from:\n {}'.format(path))
#%%
gfile='grid.xyz' # Grid file name
files=p3d.getFileNames(path=path,pattern=sPattern)
sfile=files[0]#'solTavgCf+tw+Cp.q'# Solution file name

fl=p3d.flow(path,"grid.xyz",sfile) # Check-out flow object
fl.rdHdr() # Read header from grid file and set-up flow object properties
fl.rdGrid() # Read Grid
fl.rdSol(vnames=vnames) # Read solution file
flInfo=fl.rdFlowInfo() # Read solution file info, i.e. Mach, AoA, Re and time
#%% Compute gradients
nBlocks=fl.nbk

for nb in range(nBlocks):
    fl.blk[nb].derive('p','x')
    fl.blk[nb].derive('p','y')
    fl.blk[nb].derive('p','z')
#%% Save function file
fl.wrFun(vnames=['dpdx','dpdy','dpdz'],ffile='gradP.f')