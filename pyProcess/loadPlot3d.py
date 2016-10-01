# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
from lib.myPlot3dOperator import *
import lib.Plot3DClasses as p3d
from scipy import stats
pi=np.pi
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')

import importlib
importlib.reload(p3d)

#%%
user='rpt1g12'
path='/home/'+user+'/Desktop/post/8A15W11AoA20/vsmallDomain/'
fl=p3d.flow(path,"grid.xyz","solTA.qa")

#vnames=['Cf','twx','twy','twz','Cp']
vnames=['r','u','v','w','p']
varname=vnames[4];block=4

fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames)


#%%
zlim=range(37,200,50)
bk=fl.blk[block]
surf=bk.getSubset(ylim=[0])
x=bk.var['x'].getValues()

z=bk.var['z'].getValues()

vr=bk.var[varname].getValues()
size=bk.var[varname].getSize()
size.append(5)
vr2=np.zeros(size)
for nb in range(fl.nbk):
    for n in vnames:
        vr2=fl.blk[nb].var[n].getValues()
        lze=fl.blk[nb].size[2]
        for k in range(lze):
            kk=k-76
            fl.blk[nb].var[n].val[:,:,k]=vr2[:,:,kk].copy()

#%%
#fl.wrSol("solTA_shifted.qa")
#%%
Cpavg=surf.var[varname].avgDir()[:,0]
xavg=surf.var['x'].avgDir()[:,0]


fig,ax=getFig(varname)

ii=0
for i in zlim:
    ii+=1
    ax.plot(x[:,0,i],vr[:,0,i],label='T{:d}'.format(ii))

ax.plot(xavg,Cpavg,label='avg')

handle,labels=ax.get_legend_handles_labels()
legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=3)
figc,axc=getFig(varname+' contour')
axc.contourf(x[:,0,:],z[:,0,:],fl.blk[block].var[varname].val[:,0,:])
axc.invert_yaxis()
axc.axis('equal')
#%%
xx=bk.var['x'].clone()
yy=bk.var['y'].clone()
zz=bk.var['z'].clone()
dydx=yy.derVar(0)

figd,axd=getFig('Derivative')
axd.contourf(xx.val[:,:,0],yy.val[:,:,0],dydx[:,:,0])
#%%
save=False
name='4WLE_'+varname
if (save==True):
    plt.sca(ax)
    path='pgfPlots/'
    savePlotFile(path=path+name+'.dat',vary=labels)

