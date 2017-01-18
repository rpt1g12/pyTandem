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
#cosa,sina=np.cos(np.deg2rad(5)),np.sin(np.deg2rad(5))
import importlib
#%%
importlib.reload(p3d)
#%%
path='/home/'+user+'/Desktop/post/naca0012G2/k25/'


files=p3d.getFileNames(path=path)

vnames=['Q','wx','wy','wz','Cp'];Re=125000;M=0.4

fl=p3d.flow(path,"grid.xyz",'solTQ+W+Cp0000.q')
fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames)
#%%
visible=False
x0,y0=-0.2589,0.0701;H=1e-8
fl.setTouch(x0=x0,y0=y0,R=0.020,H=H)
f,a,im=fl.contourf(varname='fu',vmin=-H,vmax=H,nlvl=256,bar=False)
fl.drawMeshPlane(ax=a)
a.set_aspect('equal')
#a.autoscale(tight=True)
a.set_xlim(x0-0.1,x0+0.1)
a.set_ylim(y0-0.1,y0+0.1)
#a.get_xaxis().set_visible(visible)
#a.get_yaxis().set_visible(visible)
#f.set_frameon(False)
f.tight_layout(pad=0)
fl.wrSol(sfile='force.q',vnames=['FA','fu','fv','fU','Q'])
