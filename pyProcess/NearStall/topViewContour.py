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
tSteps=range(0,32*12,32)
#tSteps=[0]
save=False
sPattern='solT*.q'
vnames=['r','u','v','w','p']; # Variable names
vname='twz'
plane=1
autoMinMax=False
bar=True
#%% Simulation Options
A=15 #WLE Amplitude, if SLE A=0
AoA=10 #Angle of Attack
kk=0 #Spanwise slice to look at
nwave=8 #Number of LE wavelengths
Lz_2=nwave*0.11/2
if plane==1:
    yrange=(-Lz_2,Lz_2)
    xrange=(-0.5-A/1000.0,0.5)
elif plane==2:
    yrange=(0,0.2)
    xrange=(-0.5-A/1000.0,0.5)
#%% Paths set-up
if A>0:
    #tSteps=np.asarray(range(0,64*3*11+1,64*3))  
    tSteps=range(11)
    wavy=True
    sfolder='{}WLESTA'.format(nwave)
    subpath='heaving/ss001/STA/'
    block=1 #Block to look at
    if vname=='twx':
        cmap=plt.cm.bwr
        vmin,vmax=(-1e-4,1e-4)
        nlvl=4
    elif vname=='twz':
        cmap=plt.cm.coolwarm
        #vmin,vmax=(-1e-3,1e-3)
        vmin,vmax=(-5e-4,5e-4)
        xrange=(-0.5-A/1000.0,-0.25)
        nlvl=6
    elif vname=='Cp':
        cmap=plt.cm.hot
        vmin,vmax=(-4.0,0.0)
        nlvl=10
    else:
        cmap=plt.cm.jet
        autoMinMax=True
        nlvl=21
else:
    tSteps=range(0,32*12,32)
    wavy=False
    sfolder='{}SLE'.format(nwave)
    subpath='heaving/ss003/'
    block=4 #Block to look at
    if vname=='twx':
        cmap=plt.cm.bwr
        vmin,vmax=(-5e-4,5e-4)
        nlvl=6
    elif vname=='twz':
        cmap=plt.cm.bwr
        vmin,vmax=(-1e-4,1e-3)
        nlvl=6
    elif vname=='Cp':
        cmap=plt.cm.hot
        vmin,vmax=(-3.5,0)
        nlvl=11
    else:
        cmap=plt.cm.jet
        autoMinMax=True
        nlvl=21
if plane==1:
    view='Top'
elif plane==2:
    view='Side'
simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/dellHDD/post/{}/{}".format(user,simfolder,subpath)
spath='/home/rpt1g12/Documents/thesis/figures/nearStall/{}{}{}/'.format(sfolder,vname,view)
if not os.path.exists(spath) and save:
    os.makedirs(spath)
print('Reading data from:\n {}'.format(path))
#%%
gfile='grid.xyz' # Grid file name
files=p3d.getFileNames(path=path,pattern=sPattern)
sfile=files[tSteps[0]]#'solTavgCf+tw+Cp.q'# Solution file name

fl=p3d.flow(path,gfile,sfile) # Check-out flow object
fl.rdHdr() # Read header from grid file and set-up flow object properties
fl.rdGrid() # Read Grid
fl.rdSol(vnames=vnames) # Read solution file
flInfo=fl.rdFlowInfo() # Read solution file info, i.e. Mach, AoA, Re and time
mach=flInfo[0]
#%%
xwall=fl.blk[block].var['x'].getValues()[0,kk,:]
ywall=fl.blk[block].var['z'].getValues()[0,kk,:]

for ii in tSteps:
    fl.rdSol(vnames=vnames,sfile=files[ii])
    flInfo=fl.rdFlowInfo()
    if vname=='Cp':
        cp=fl.blk[block].var['p'].getValues()
        cp=(cp-1/1.4)/(0.5*mach**2)
        fl.blk[block].setData(vname=vname,val=cp)
    if vname not in vnames+['Cp']:
        fl.blk[block].getWSS(Re=flInfo[2])
    var=fl.blk[block].var[vname].getValues()
    if autoMinMax:
        vmin=var.min()#-0.5*(np.abs(var.min())+var.max())
        vmax=var.max()
    if save:
        bar=False
    f,a,im=fl.blk[block].contourf(varname=vname,vmin=vmin,vmax=vmax,plane=plane,k=kk,nlvl=nlvl,cmap=cmap,bar=bar);nfig+=1
    figs.append(f);axs.append(a)
    fl.blk[block].contour(varname=vname,vmin=vmin,vmax=vmax,plane=plane,k=kk,nlvl=nlvl,ax=axs[nfig])
    axs[nfig].plot(xwall,ywall,lw=2,color='k')
    if not save:
        anotation=r'$t={:3.4f}$'.format(flInfo[3])
        axs[nfig].text(0.0,-0.1,anotation,ha='center',va='bottom',transform=axs[nfig].transAxes)
    axs[nfig].set_xlim(xrange)
    axs[nfig].set_ylim(yrange)
    axs[nfig].set_aspect('equal')
    if save:
        saveFigOnly(path=spath,fig=figs[nfig],ax=axs[nfig],name='{}{:04d}'.format(vname,ii),ext='.pdf')       
#%%
print('max={},min={}'.format(var.max(),var.min()))