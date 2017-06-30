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
A=0 #WLE Amplitude, if SLE A=0
AoA=10 #Angle of Attack
nwave=8 #Number of WLE wavelenghts

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
    subpath='average/ss005/'
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
#%%SubSet set-up
ss=p3d.flow(sspath,"grid.xyz",files[0])
ssBlocks=[]
nxi=[];neta=[];nzeta=[]
for kb in range(bkXYZ[2]):    
    for jb in range(bkXYZ[1]):
        for ib in range(bkXYZ[0]):
            nb=kb*bkXYZ[1]*bkXYZ[0]+jb*bkXYZ[0]+ib
            ssBlocks.append(fl.blk[blocks[nb]].getSubset(xlim=xRanges[ib],ylim=yRanges[jb],zlim=zRanges[kb],link=True))
            if kb==0:            
                neta.append(ssBlocks[-1].size[1])
                if jb==0:            
                    nxi.append(ssBlocks[-1].size[0])
    nzeta.append(ssBlocks[-1].size[2])
#%%
ss.setHdr(bkXYZ,nxi,neta,nzeta)

for nb in range(len(blocks)):
    ss.blk[nb]=ssBlocks[nb]
#%%
ss.wrGrid()
ss.wrSol(vnames=vnames,sfile=sfile,flInfo=flInfo)