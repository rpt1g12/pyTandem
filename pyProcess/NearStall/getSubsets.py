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
importlib.reload(p3d)

#%%
save=False

#%% Sub-set settings
ssName='upSurface'
bkXYZ=[1,1,1]
blocks=[4]
xRanges=[range(110)]
yRanges=[[1]]
zRanges=[None]
vnames=['r','u','v','w','p']; # Variable names
sPattern='solT*.*.q'
var='p'    
#%% Paths set-up
path='/media/rpt1g12/My Passport/results/quickPitchUp/ss001/T45T165/'
sspath="/media/{}/dellHDD/post/4A15W11AoA10/heaving/{}/".format(user,ssName)
#sspath=path+ssName+'/'
if not os.path.exists(sspath) and save:
    os.makedirs(sspath)
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
ss=fl.getSubsets(fromBlocks=blocks,bkXYZ=bkXYZ,xRanges=xRanges,yRanges=yRanges,zRanges=zRanges,ssName=ssName,sspath=sspath,ssfile=files[0])

ss.wrGrid()

#%% Loop trough time 1731
for sfile in files[1731:]:
    fl.rdSol(vnames=vnames,sfile=sfile)
    flInfo=fl.rdFlowInfo()
    fName='{}T{:03.4f}.f'.format(var,flInfo[3])
    ss.blk[:]=fl.getSubsets(fromBlocks=blocks,bkXYZ=bkXYZ,xRanges=xRanges,yRanges=yRanges,zRanges=zRanges,
                            ssName=ssName,sspath=sspath,ssfile=sfile,ssfl=ss)
    ss.wrSol(vnames=vnames,sfile=sfile,flInfo=flInfo)
    #ss.wrFun(vnames=[var],ffile=fName)