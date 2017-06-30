# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
from lib.myPlot3dOperator import *
import lib.Plot3DClasses as p3d
from scipy import stats
from scipy.interpolate import interp2d
from scipy.integrate import odeint
pi=np.pi
import getpass
user=getpass.getuser()

figs=[];axs=[];hdls=[];lbls=[];lgds=[];nfig=-1;fs=18
import importlib
#%%
importlib.reload(p3d)

#%%
path="/media/rpt1g12/dellHDD/post/4A15W11AoA20/average/"

vnames=['r','u','v','w','p'];

fl=p3d.flow(path,"grid.xyz",'solTA.qa')
fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames);flInfo=fl.rdFlowInfo()
#%%
#fl.blk[4].getMetrics()
#J=fl.blk[4].J.getValues()[:,0,:]
#%%
plt.close('all')
ss=fl.blk[4].getSubset(ylim=[1])
f,a,im=ss.contourf(varname='u',k=0,plane=1,vmin=-0.05,vmax=0.05)
a.invert_yaxis()
a.set_aspect('equal')
ss.drawMeshPlane(direction=1,skp=(5,1),showBlock=True,ax=a)
#%%
x=ss.var['x'].getValues()[:,0,:]
z=ss.var['z'].getValues()[:,0,:]
u=ss.var['u'].getValues()[:,0,:]
w=ss.var['w'].getValues()[:,0,:]
#u=u/np.mean(np.abs(u))
#w=w/np.mean(np.abs(w))
#un,wn=normalise(u,w)
nxi,nze=u.shape
xi,ze=np.linspace(0,nxi-1,nxi),np.linspace(0,nze-1,nze)
x_xize=interp2d(ze,xi,x)
z_xize=interp2d(ze,xi,z)
u_xize=interp2d(ze,xi,u)
w_xize=interp2d(ze,xi,w)

uw = lambda p,t: [u_xize(p[1],p[0])[0],-w_xize(p[1],p[0])[0]]
#%%
nlines=0
critPt=[]
#%%
critPt.append(np.array([10,87]))
#%%
iDelta=13
kDelta=-4
iRange=[critPt[-1][0]+iDelta]
kRange=[critPt[-1][1]+kDelta]
dt = -10; t0 = 0; nt=2000
t1 = nt*dt
t = np.linspace(t0,t1,nt)
xz=np.zeros((2,nt))
for kk in kRange:
    for ii in iRange:
        print('\nStart indices: ({:3.2f},{:3.2f})'.format(ii,kk))
        streamline=odeint(uw,(ii,kk),t)
        for i in range(nt):
            xz[0,i]=x_xize(streamline[i,1],streamline[i,0])
            xz[1,i]=z_xize(streamline[i,1],streamline[i,0])      
        a.plot(xz[0,:],xz[1,:],lw=2,color='green')
        line=a.lines[-1];nlines+=1
axShow(a)        

print('End indices: ({:3.2f},{:3.2f})'.format(streamline[-1,0],streamline[-1,1]))
print('Critical Points:\n')
for n in range(len(critPt)):
    print('n = {:2} : ({:+04.2f},{:+04.2f})'.format(n,critPt[n][0],critPt[n][1]))
print('\n')

#%% Clean last line
a.lines.pop(-1);axShow(a)
#%% Clean all lines
for i in range(nlines):
    a.lines.pop(-1)
axShow(a);nlines=0

#%%
critPt.append(np.array([iRange[0],kRange[0]]));
a.scatter(xz[0,-1],xz[1,-1]);axShow(a)
#%% Append Last Point
critPt.pop(-1)
critPt.append(streamline[-1,:]);
a.scatter(xz[0,-1],xz[1,-1]);axShow(a)
print('Critical Points:\n')
for n in range(len(critPt)):
    print('n = {:2} : ({:+04.2f},{:+04.2f})'.format(n,critPt[n][0],critPt[n][1]))
#%% Append First point
critPt.pop(-1)
critPt.append(np.array([iRange[0],kRange[0]]));
a.scatter(xz[0,0],xz[1,0]);axShow(a)
print('Critical Points:\n')
for n in range(len(critPt)):
    print('n = {:2} : ({:+04.2f},{:+04.2f})'.format(n,critPt[n][0],critPt[n][1]))
#%%
print('Critical Points:\n')
for n in range(len(critPt)):
    print('n = {:2} : ({:+04.2f},{:+04.2f})'.format(n,critPt[n][0],critPt[n][1]))
#%% Save result
#np.savetxt(path+'/critPoints.dat',np.asarray(critPt),delimiter='\t',header='xcrit\tzcrit')