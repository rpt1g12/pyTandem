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
def getCirc(x0,y0,r,theta0,theta1,thetas):
    theta=np.linspace(theta0,theta1,thetas)
    fctr=np.pi/180;theta*=fctr
    ntheta=theta.size
    x=np.zeros((theta.size,1,1))
    y=np.zeros((theta.size,1,1))
    for n in range(ntheta):
        x[n,0,0]=(r*np.cos(theta[n]))+x0
        y[n,0,0]=(r*np.sin(theta[n]))+y0
    return x,y
    
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
#xlim=(-0.5,0.5);ylim=(0,1);kplane=125
#x,y,z=np.mgrid[xlim[0]:xlim[1]:100j,ylim[0]:ylim[1]:100j,0:0:1j]
#lvl=np.linspace(0.5,1,51)
#print('Interpolating...')
#iblock=bk.interpolate2dk('p',x,y,kplane,method='cubic')
#ifig,iax=getFig('Interpolated')
#iax.contourf(iblock.var['x'].val[:,:,0],iblock.var['y'].val[:,:,0],iblock.var['p'].val[:,:,0],levels=lvl)
#sfig,sax=getFig('Original')
#sax.contourf(bk.var['x'].val[:,:,kplane],bk.var['y'].val[:,:,kplane],bk.var['p'].val[:,:,kplane],levels=lvl)
#sax.set_xlim(xlim[0],xlim[1])
#sax.set_ylim(ylim[0],ylim[1])
#%%
mfl=fl.mergeBlocks(3,2)
#%%
save=False
name='4WLE_'+varname
if (save==True):
    plt.sca(ax)
    path='pgfPlots/'
    savePlotFile(path=path+name+'.dat',vary=labels)

