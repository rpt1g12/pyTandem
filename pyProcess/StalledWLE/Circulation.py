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
import getpass
user=getpass.getuser()
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
path='/home/'+user+'/Desktop/post/8A00W11AoA20/fullDomain/'
fl=p3d.flow(path,"grid.xyz","solTA.qa")


vnames=['r','u','v','w','p']

fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames)

#%%
normz=(fl.blk[4].var['z'].getValues()[0,0,:])/(0.11)
mfl=fl.mergeBlocks(3,2)
del fl
#%%
mfig,mrax=getFig('Merged')
mrax.contourf(mfl.blk[0].var['x'].val[:,:,0],mfl.blk[0].var['y'].val[:,:,0],mfl.blk[0].var['u'].val[:,:,0])
r=1;nr=360;x0=0.0;y0=0.0
xc,yc,theta=getCirc(x0,y0,r,0,359,nr)
#%%
ncirc=201#8*4+1
circ=[]
#%%
for kk in range(ncirc):
    circle=mfl.blk[0].interpolate2dk('u',xc,yc,kk)
    v=mfl.blk[0].interpolate2dk('v',xc,yc,kk,mode='values')
    circle.setData(vname='v',val=v)
    circle.setData(vname='theta',val=theta)

    dl=(2*pi*r/nr)
    rad=np.deg2rad(circle.var['theta'].getValues())
    tx=dl*np.sin(rad)
    ty=-dl*np.cos(rad)
    circle.setData(vname='tx',val=tx)
    circle.setData(vname='ty',val=ty)
    
    cr=np.sum(circle.var['tx'].getValues()*circle.var['u'].getValues(True)+
               circle.var['ty'].getValues()*circle.var['v'].getValues(True))
    circ.append(cr)
    print('circulation k={:d}={:8.5f}'.format(kk,cr))

#%%
circulation=circ.copy()
for k in range(201):
    kk=k+(124)
    if (kk>200):
        kk-=200
    circulation[k]=circ[kk]  
#%%

crfig,crax=getFig('8SLE_circulation')
crax.plot(normz,circulation,label='circulation')
handle,labels=crax.get_legend_handles_labels()
legend=crax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=1)
fit(crax)

#%%
save=False
name=crfig.canvas.get_window_title()
ax=crax
if (save==True):
    plt.sca(ax)
    path='pgfPlots/'
    savePlotFile(path=path+name+'.dat',vary=labels)


