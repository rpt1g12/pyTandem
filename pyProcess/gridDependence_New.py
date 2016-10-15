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
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')

import importlib
#%%
importlib.reload(p3d)
user='rpt1g12'
path='/media/'+user+'/My Passport/results/6Blocks/post/'
grid=['G1','G2','G3','G4']
ll=0
#%%
vnames=['Cf','twx','twy','twz','Cp']
    
varname=vnames[4];block=4
for ll in range(len(grid)):
    fl=p3d.flow(path+grid[ll]+'/',"grid.xyz","solTCf+tw+Cp.qa")
    

    
    fl.rdHdr()
    fl.rdGrid()
    fl.rdSol(vnames=vnames)
    #%%
    bk=fl.blk[block]
    #%%
    surf=bk.getSubset(ylim=[0])
    
    vavg=surf.var[varname].avgDir()[:,0]
    xavg=surf.var['x'].avgDir()[:,0]
    if (ll==0):
        f0,a0=getFig('grid{}'.format(varname))
    ax=plt.gca()
    ax.plot(xavg,vavg,lw=2,label=grid[ll])
handle,labels=ax.get_legend_handles_labels()
legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=1)
fit(ax)
#%%
save=True
name=ax.figure.canvas.get_window_title()
if (save==True):
    plt.sca(ax)
    path='pgfPlots/'
    savePlotFile(path=path+name+'.dat',vary=labels)