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
opt=80
sname='solTA.qa'
vnames=['Cf','twx','twy','twz','Cp']    
varname=vnames[4];block=4;
pinf=0#1/1.4
kinf=1#0.5*0.3**2
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
    sim='8A15W11AoA20/fullDomain/'
    name='8WLE'
    section=[19,44,69,119]
    sName=['T4','T5','T6','T8']
    fctr=1
elif opt==80:
    sim='8A00W11AoA20/fullDomain/'
    name='8SLE'
    section=[19,44,69,119]
    sName=['T4','T5','T6','T8']
    fctr=np.sqrt(1-0.3**2)
#%%
fl=p3d.flow(path+sim,"grid.xyz",sname)     
fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames)
bk=fl.blk[block]
bk0=fl.blk[1]
surf=bk.getSubset(ylim=[0])
surf0=bk0.getSubset(ylim=[120])

#%%    
#fg,ax=getFig('{}_{}'.format(name,varname))
#%%
ii=0;k0=1;k1=-1

v=((k0*surf0.var[varname].avgDir()-pinf)/kinf+(k1*surf.var[varname].avgDir()-pinf/kinf))*fctr
#ax.plot(surf.var['x'].getValues()[:,0,k],surf.var[varname].getValues()[:,0,k],lw=2,label=sName[ii])
ax.plot(surf.var['x'].avgDir(),v,lw=2,label=name)


ax.grid(True)
handle,labels,legend=getLabels(ax,1,(1,1))
#if varname=='Cp':
#    ax.invert_yaxis()
axShow(ax)
#%%
save=False
name=ax.figure.canvas.get_window_title()
if (save==True):
    plt.sca(ax)
    path='pgfPlots/'
    savePlotFile(path=path+name+'.dat',vary=labels)
