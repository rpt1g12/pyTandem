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
sPattern='*.qa'
vnames=['Q','wx','wy','wz','D']; # Variable names
#%% Options
save=False
A=15 #WLE Amplitude, if SLE A=0
AoA=10 #Angle of Attack
block=4 #Block to look at
krange=range(1,49) #Spanwise slice to look at
nwave=1 #Number of LE wavelengths
#%% xLSB set-up
if A>0:
    sfolder='{}WLE'.format(nwave)
    subpath='average/ss006/'
else:
    sfolder='{}SLE'.format(nwave)
    subpath='average/ss004/'

#%% Paths set-up
simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/dellHDD/post/{}/{}".format(user,simfolder,subpath)
spath='/home/rpt1g12/Documents/thesis/figures/lowAoA/OmegaX{:d}/'.format(AoA)
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
#%%
if AoA==0:
    vmin=-5;vmax=5
elif AoA==6:
    vmin=-15;vmax=15
elif AoA==10:
    vmin=-15;vmax=15
f,a,im=fl.blk[block].contourf(varname='wx',vmin=vmin,vmax=vmax,k=kk,nlvl=21,cmap=plt.cm.bwr,bar=False)
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
#up.contour(varname='Cp',vmin=vmin,vmax=vmax,k=kk,nlvl=21,ax=axs[nfig])

a.set_aspect('equal')
axs[nfig].set_xlim(-0.5,0.5)
axs[nfig].set_ylim(0.0,0.25)
if save:
    saveFigOnly(path=spath,fig=figs[nfig],ax=axs[nfig],name='Cp{}{:02d}'.format(sfolder,AoA),ext='.pdf')