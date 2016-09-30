# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
import lib.Plot3DClasses as p3d
from scipy import stats
pi=np.pi
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')

import importlib
importlib.reload(p3d)

#%%
#fl=p3d.flow("/home/rpt1g12/Desktop/post/8A00W11AoA20/surface/","grid.xyz","solTavgCf+tw+Cp0641.q")
fl=p3d.flow("/home/rperezt/Desktop/post/8A15W11AoA20/surface/","grid.xyz","solTavgCf+tw+Cp0480.q")
#fl=p3d.flow("/home/rpt1g12/Desktop/post/A4A15W11AoA20/f64/","grid.xyz","solTCf+tw+Cp.qa")

fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=['Cf','twx','twy','twz','Cp'])

#%%
zlim=range(19,200,25)
varname='Cp';bk=fl.blk[0]
surf=bk.getSubset(ylim=[0])
x=bk.var['x'].getValues()

z=bk.var['z'].getValues()

Cp=bk.var[varname].getValues()

Cpavg=surf.var[varname].avgDir()[:,0]
xavg=surf.var['x'].avgDir()[:,0]


fig,ax=getFig(varname)

ii=0
for i in zlim:
    ii+=1
    ax.plot(x[:,0,i],Cp[:,0,i],label='T{:d}'.format(ii))

ax.plot(xavg,Cpavg,label='avg')

handle,labels=ax.get_legend_handles_labels()
legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=3)
#plt.contourf(x[:,0,:],z[:,0,:],Cp[:,0,:])
#plt.axis('equal')

#%%
save=True
name='8WLE_'+varname
if (save==True):
    plt.sca(ax)
    path='pgfPlots/'
    savePlotFile(path=path+name+'.dat',vary=labels)

