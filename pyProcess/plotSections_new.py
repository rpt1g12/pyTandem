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
from lib.matplotlib2tikz import save as tikz_save
#plt.close('all')

import importlib
#%%
importlib.reload(p3d)
user='rpt1g12'
path='/home/'+user+'/Desktop/post/'
opt=4
vnames=['Cf','twx','twy','twz','Cp']    
varname=vnames[0];block=4;
if opt==2:
    sim='A15W11AoA20/'
    name='2WLE'
    section=[37,87]
    sName=['T1','T2']
elif opt==3:
    sim='3A15W11AoA20/'
    name='3WLE'
    section=[37,87,137]
    sName=['T1','T2','T3']
elif opt==4:
    sim='A4A15W11AoA20/f64/'
    name='4WLE'
    section=[37,87,137,187]
    sName=['T1','T2','T3','T4']
elif opt==8:
    sim='8A15W11AoA20/vsmallDomain/'
    name='8WLE'
    section=[19,44,69,119]
    sName=['T4','T5','T6','T8']
#%%
fl=p3d.flow(path+sim,"grid.xyz","solTCf+tw+Cp.qa")     
fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames)
bk=fl.blk[block]
surf=bk.getSubset(ylim=[0])
#%%    
fg,ax=getFig('{}_{}'.format(name,varname))
#%%
ii=0
for k in section:
    ax.plot(surf.var['x'].getValues()[:,0,k],surf.var[varname].getValues()[:,0,k],lw=2,label=sName[ii])
    ii+=1

ax.grid(True)
handle,labels,legend=getLabels(ax,1,(1,1))
fit(ax)
#%%
save=True
name=ax.figure.canvas.get_window_title()
if (save==True):
    plt.sca(ax)
    path='pgfPlots/'
    savePlotFile(path=path+name+'.dat',vary=labels)
