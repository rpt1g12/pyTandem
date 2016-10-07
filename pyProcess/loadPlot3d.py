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
    thetad=theta.copy()
    fctr=np.pi/180;theta*=fctr
    ntheta=theta.size
    x=np.zeros((theta.size,1,1))
    y=np.zeros((theta.size,1,1))
    for n in range(ntheta):
        x[n,0,0]=(r*np.cos(theta[n]))+x0
        y[n,0,0]=(r*np.sin(theta[n]))+y0
    return x,y,thetad
    
#%%
user='rpt1g12'
path='/home/'+user+'/Desktop/post/8A15W11AoA20/fullDomain/'
fl=p3d.flow(path,"grid.xyz","solTA.qa")

#vnames=['Cf','twx','twy','twz','Cp']
vnames=['r','u','v','w','p']
zlim=range(37,200,50)
varname=vnames[4];block=4

fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames)


#%%
#bk=fl.blk[block]
#%%
#surf=bk.getSubset(ylim=[0])
#%%
sfl=fl.shiftK(125)
#%%
sfl.wrSol("solTA_shiftedTest.qa")
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
#    ax.plot(bk.var['x'].getValues()[:,0,i],bk.var[varname].getValues()[:,0,i],label='T{:d}'.format(ii))
#
#ax.plot(xavg,Cpavg,label='avg')
#
#handle,labels=ax.get_legend_handles_labels()
#legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=3)
#figc,axc=getFig(varname+' contour')
#axc.contourf(bk.var['x'].getValues()[:,0,:],bk.var['z'].getValues()[:,0,:],sfl.blk[block].var[varname].val[:,0,:])
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
#normz=(fl.blk[4].var['z'].getValues()[0,0,:])/(0.11)
#mfl=fl.mergeBlocks(3,2)
#del fl
##%%
#mfig,mrax=getFig('Merged')
#cfig,cax=getFig('Circle')
#mrax.contourf(mfl.blk[0].var['x'].val[:,:,0],mfl.blk[0].var['y'].val[:,:,0],mfl.blk[0].var['u'].val[:,:,0])
#r=1;nr=360;x0=0.0;y0=0.0
#xc,yc,theta=getCirc(x0,y0,r,0,359,nr)
##%%
#ncirc=201#8*4+1
#circ=[]
##%%
#for kk in range(ncirc):
#    #kk=round(nk*6.25)
#    circle=mfl.blk[0].interpolate2dk('u',xc,yc,kk)
#    v=mfl.blk[0].interpolate2dk('v',xc,yc,kk,mode='values')
#    circle.setData(vname='v',val=v)
#    circle.setData(vname='theta',val=theta)
#
#    dl=(2*pi*r/nr)
#    rad=np.deg2rad(circle.var['theta'].getValues())
#    tx=dl*np.sin(rad)
#    ty=-dl*np.cos(rad)
#    circle.setData(vname='tx',val=tx)
#    circle.setData(vname='ty',val=ty)
#    
#    cr=np.sum(circle.var['tx'].getValues(True)*circle.var['u'].getValues(True)+
#               circle.var['ty'].getValues(True)*circle.var['v'].getValues(True))
#    circ.append(cr)
#    print('circulation k={:d}={:8.5f}'.format(kk,cr))
#
##%%
#circulation=circ.copy()
#for k in range(201):
#    kk=k+(125-201)
#    circulation[k]=circ[kk]  
##%%
#
#crfig,crax=getFig('8SLE_circulation')
#crax.plot(normz,circulation,label='circulation')
#handle,labels=crax.get_legend_handles_labels()
#legend=crax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=1)
#fit(crax)
#%%
#
#mrax.scatter(circle.var['x'].getValues(True),circle.var['y'].getValues(True))
##mrax.quiver(circle.var['x'].getValues(),circle.var['y'].getValues(),circle.var['tx'].getValues(),circle.var['ty'].getValues())
#cax.plot(circle.var['theta'].getValues(True),circle.var['u'].getValues(True),lw=2,color='blue',label='u')
#cax.plot(circle.var['theta'].getValues(True),circle.var['v'].getValues(True),lw=2,color='red',label='v')
#handle,labels=cax.get_legend_handles_labels()
#legend=cax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=1)
##mrax.axis('equal')
#mrax.set_xlim(-1.5,1.5)
#mrax.set_ylim(-1.5,1.5)
#mrax.set_xlabel('x',fontsize=20)
#mrax.set_ylabel('y',fontsize=20)
#cax.set_xlabel(r'$\theta$',fontsize=20)
#cax.set_ylabel(r'$u$ & $v$',fontsize=20)
#fit(cax)
#%%
save=False
name=crfig.canvas.get_window_title()
ax=crax
if (save==True):
    plt.sca(ax)
    path='pgfPlots/'
    savePlotFile(path=path+name+'.dat',vary=labels)

