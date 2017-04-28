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
from scipy.interpolate import interp2d
from scipy.integrate import odeint
pi=np.pi
import getpass
user=getpass.getuser()
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')
figs=[];axs=[];hdls=[];lbls=[];lgds=[];nfig=-1;fs=18
import importlib
#%%
importlib.reload(p3d)

#%%
path="/media/rpt1g12/dellHDD/post/8A15W11AoA20/average/vsmallDomain/"

vnames=['r','u','v','w','p'];

fl=p3d.flow(path,"grid.xyz",'solTA.qa')
fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames);flInfo=fl.rdFlowInfo()
#%%
#fl.blk[4].getMetrics()
#J=fl.blk[4].J.getValues()[:,0,:]
#%%

ss=fl.blk[4].getSubset(ylim=[1])
f,a,im=ss.contourf(varname='u',k=0,plane=1,vmin=-0.05,vmax=0.05)
ss.drawMeshPlane(direction=1,skp=(5,1),showBlock=True,ax=a)
#%%
x=ss.var['x'].getValues()[:,0,:]
z=ss.var['z'].getValues()[:,0,:]
u=ss.var['u'].getValues()[:,0,:]
w=ss.var['w'].getValues()[:,0,:]
u=u/np.mean(np.abs(u))
w=w/np.mean(np.abs(w))
#u,w=normalise(u,w)
nxi,nze=u.shape
xi,ze=np.linspace(0,nxi-1,nxi),np.linspace(0,nze-1,nze)
x_xize=interp2d(ze,xi,x)
z_xize=interp2d(ze,xi,z)
u_xize=interp2d(ze,xi,u)
w_xize=interp2d(ze,xi,w)

uw = lambda p,t: [u_xize(p[1],p[0])[0],w_xize(p[1],p[0])[0]]

dt = 1
t0 = 0
nt=100
t1 = nt*dt
t = np.linspace(t0,t1,nt)
xz=np.zeros((2,nt))
for kk in range(12,25,2):
    for ii in range(20,81,4):
        print(ii,kk)
        streamline=odeint(uw,(ii,kk),t)
        for i in range(nt):
            xz[0,i]=x_xize(streamline[i,1],streamline[i,0])
            xz[1,i]=z_xize(streamline[i,1],streamline[i,0])      
        a.plot(xz[0,:],xz[1,:],lw=2,color='green')
axShow(a)

#%%

def normalise(u,v):
    m=np.sqrt(u**2+v**2)
    u=u/m
    v=v/m
    return u,v

a.figure