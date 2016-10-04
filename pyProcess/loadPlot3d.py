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
#%%
#surf=bk.getSubset(ylim=[0])
#x=bk.var['x'].getValues()
#
#z=bk.var['z'].getValues()
#
#vr=bk.var[varname].getValues()
#sfl=fl.shiftK(125)

#%%
#sfl.wrSol("solTA_shiftedTest.qa")
#%%
#Cpavg=surf.var[varname].avgDir()[:,0]
#xavg=surf.var['x'].avgDir()[:,0]
#
#
#fig,ax=getFig(varname)
#
#ii=0
#for i in zlim:
#    ii+=1
#    ax.plot(x[:,0,i],vr[:,0,i],label='T{:d}'.format(ii))
#
#ax.plot(xavg,Cpavg,label='avg')
#
#handle,labels=ax.get_legend_handles_labels()
#legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=3)
#figc,axc=getFig(varname+' contour')
#axc.contourf(x[:,0,:],z[:,0,:],sfl.blk[block].var[varname].val[:,0,:])
#axc.invert_yaxis()
#axc.axis('equal')
#%%
#bk.getMetrics()
#
#vname1='w';vname2='y'
#dname='d{}d{}'.format(vname1,vname2)
#bk.derive(vname1,vname2)
#
#figd,axd=getFig(dname)
#axd.contourf(bk.var['x'].val[:,0,:],bk.var['z'].val[:,0,:],bk.var[dname].val[:,0,:])
#%%
x,y,z=np.mgrid[-0.5:0.5:100j,0.1:0.1:1j,-0.4:0.4:50j]
print('Interpolating...')
iblock=bk.interpolate('p',x,y,z,method='nearest')
plt.contourf(iblock.var['x'].val[:,0,:],iblock.var['z'].val[:,0,:],iblock.var['p'].val[:,0,:])
#%%
save=False
name='4WLE_'+varname
if (save==True):
    plt.sca(ax)
    path='pgfPlots/'
    savePlotFile(path=path+name+'.dat',vary=labels)

