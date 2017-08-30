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
avgL=int(64*0.5)
save=True
sPattern='solT*.q'
vnames=['r','u','v','w','p']; # Variable names
#%% Simulation Options
A=15 #WLE Amplitude, if SLE A=0
AoA=10 #Angle of Attack
nwave=8 #Number of LE wavelengths
ssl=True
#%% Paths set-up
if A>0:
    wavy=True
    sfolder='{}WLE'.format(nwave)
    if ssl:
        subpath='heaving/ss005/'
        sfolder+='SSL'
    else:
        subpath='heaving/ss006/'
        sfolder+='Central'
else:
    wavy=False
    sfolder='{}SLE'.format(nwave)
    subpath='heaving/ss003/'

simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/dellHDD/post/{}/{}".format(user,simfolder,subpath)
spath=path+'STA/'
if not os.path.exists(spath) and save:
    os.makedirs(spath)
print('Reading data from:\n {}'.format(path))
#%% Flow object set up
gfile='grid.xyz' # Grid file name
files=p3d.getFileNames(path=path,pattern=sPattern)
sfile=files[0]#'solTavgCf+tw+Cp.q'# Solution file name

fl=p3d.flow(path,gfile,sfile) # Check-out flow object
fl.rdHdr() # Read header from grid file and set-up flow object properties
fl.rdGrid() # Read Grid
fl.rdSol(vnames=vnames) # Read solution file
flInfo=fl.rdFlowInfo() # Read solution file info, i.e. Mach, AoA, Re and time
mach=flInfo[0]
#%%
havgL=int(avgL/2)
nt=len(files)
for n in range(int(nt/avgL)):
    sfile=files[havgL+n*avgL]
    flInfo=fl.rdFlowInfo(sfile=sfile)
    flavg=fl.getAvg(files=files[n*avgL:(n+1)*avgL+1],vnames=vnames)
    flavg.wrSol(path=spath,vnames=vnames,sfile=sfile,flInfo=flInfo)