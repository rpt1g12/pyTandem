# -*- coding: utf-8 -*-
import os
import math
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
var='p';
save=True
A=15 #WLE Amplitude, if SLE A=0
AoA=10 #Angle of Attack
block=0 #Block to look at
nwave=1
vnames=['r','u','v','w','p']; # Variable names
sPattern='solT*.*.q'
#%% Paths set-up
if A>0:
    wavy=True 
    sfolder='{}WLE'.format(nwave)
    subpath='ss004/'
else:
    wavy=False   
    subpath='ss002/'
    sfolder='{}SLE'.format(nwave)

simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/dellHDD/post/{}/unsteady/{}".format(user,simfolder,subpath)
spath='/home/rpt1g12/Documents/thesis/figures/lowAoA/{}{:02d}LSBPics/'.format(sfolder,AoA)
if not os.path.exists(spath) and save:
    os.makedirs(spath)
print('Reading data from:\n {}'.format(path))
#%%
gfile='grid.xyz' # Grid file name
files=p3d.getFileNames(path=path,pattern=sPattern)
sfile=files[62]#'solTavgCf+tw+Cp.q'# Solution file name

fl=p3d.flow(path,gfile,sfile) # Check-out flow object
fl.rdHdr() # Read header from grid file and set-up flow object properties
fl.rdGrid() # Read Grid
fl.rdSol(vnames=vnames) # Read solution file
flInfo=fl.rdFlowInfo() # Read solution file info, i.e. Mach, AoA, Re and time
mach=flInfo[0]
#%%

vmax=-1.5
vmin=-3.5

#cp=fl.blk[0].var['p'].getValues()
#cp=(cp-1/1.4)/(0.5*mach**2)
#fl.blk[0].setData(vname='Cp',val=cp)
#f,a,im=fl.blk[0].contourf(varname='Cp',vmin=vmin,vmax=vmax,k=0,nlvl=21,cmap=plt.cm.hot,bar=False);nfig+=1
#figs.append(f);axs.append(a)
#axs[nfig].set_xlim(-0.485,-0.25)
#axs[nfig].set_ylim(0,0.15)
#axs[nfig].set_aspect('equal')
#saveFigOnly(path=spath,fig=figs[nfig],ax=axs[nfig],name='p62',ext='.pdf')

for ii in range(62,93):
    fl.rdSol(vnames=vnames,sfile=files[ii])
    cp=fl.blk[0].var['p'].getValues()
    cp=(cp-1/1.4)/(0.5*mach**2)
    fl.blk[0].setData(vname='Cp',val=cp)
    f,a,im=fl.blk[0].contourf(varname='Cp',vmin=vmin,vmax=vmax,k=0,nlvl=21,cmap=plt.cm.hot,bar=False);nfig+=1
    figs.append(f);axs.append(a)
    axs[nfig].set_xlim(-0.485,-0.25)
    axs[nfig].set_ylim(0,0.15)
    axs[nfig].set_aspect('equal')
    saveFigOnly(path=spath,fig=figs[nfig],ax=axs[nfig],name='p{:2d}'.format(ii),ext='.pdf')
    
