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
cosa,sina=np.cos(np.deg2rad(5)),np.sin(np.deg2rad(5))
import importlib
#%%
importlib.reload(p3d)
#%%
visible=False;save=True;rotate=True

path='/home/'+user+'/Desktop/post/naca0012G3/k25/'


files=p3d.getFileNames(path=path)

vnames=['Q','wx','wy','wz','Cp'];Re=125000;M=0.4

fl=p3d.flow(path,"grid.xyz",'solTQ+W+Cp0000.q')
fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames)
#%%
if rotate:
    for bk in fl.blk[:]:
        x0=bk.var['x'].getValues()-0.5
        y0=bk.var['y'].getValues()
        x,y=x0*cosa+y0*sina+1.0,y0*cosa-x0*sina
        bk.setData('x',x)    
        bk.setData('y',y)
#%%

    
f,a,im=fl.contourf('wz',vmin=-60,vmax=60,nlvl=20,cmap=plt.cm.bwr,bar=False)
a.plot(fl.blk[4].var['x'].getValues()[:,0,0],fl.blk[4].var['y'].getValues()[:,0,0],color='black',lw=2)
a.plot(fl.blk[1].var['x'].getValues()[:,-1,0],fl.blk[1].var['y'].getValues()[:,-1,0],color='black',lw=2)
a.set_aspect('equal')
#a.autoscale(tight=True)
a.set_xlim(0.0,1.25)
a.set_ylim(-0.1,0.3)
a.get_xaxis().set_visible(visible)
a.get_yaxis().set_visible(visible)
f.set_frameon(False)
f.tight_layout(pad=0)
if save:
    f.savefig('/home/'+user+'/anaconda3/pyTandem/pgfPlots/wz_G3.pdf', bbox_inches='tight', pad_inches=0,dpi=300)

#%%
#def naca(x,c=1.0,t0=12):
#    """NACA 00XX series equation"""
#    t=abs(t0)*0.01
#    k=0.991148635
#    if (t0<0.0):
#       fctr=-1.0
#    else:
#       fctr=1.0
#    xc=x/c
#    y= fctr*(t*c*k/0.2) \
#    * (0.298222773*np.sqrt(xc)-0.127125232*xc-0.357907806*xc**2 \
#    +0.291984971*xc**3-0.105174606*xc**4)
#    return y
#def x2xhat(x,y):
#    a=np.deg2rad(5)
#    y=naca(x+0.5)
#    xh=(x-0.5)*np.cos(a)+y*np.sin(a)+1
#    yh=y*np.cos(a)-(x-0.5)*np.sin(a)
#    return xh,yh
#def xhat2x(xh,yh):
#    a=np.deg2rad(5)
#    A=np.array([[np.cos(a),np.sin(a)],[-np.sin(a),np.cos(a)]])
#    Ainv=np.linalg.inv(A)
#    [x,y]=np.matmul(Ainv,[xh-1,yh]
#    x+=0.5
#    return x,y