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
figs=[];axs=[];hdls=[];lbls=[];lgds=[];nfig=-1;fs=18
import importlib
#%%
importlib.reload(p3d)

#%%
path="/media/rpt1g12/dell's_hdd/post/8A15W11AoA20/average/vsmallDomain/"

vnames=['r','u','v','w','p'];

fl=p3d.flow(path,"grid.xyz",'solTA.qa')
fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames);flInfo=fl.rdFlowInfo()
#%%
ss=fl.blk[4].getSubset(ylim=[1])

x=ss.var['x'].getValues()[:,0,:]
z=ss.var['z'].getValues()[:,0,:]
u=ss.var['u'].getValues()[:,0,:]
w=ss.var['w'].getValues()[:,0,:]
nxi,nze=u.shape
xi,ze=np.linspace(0,nxi-1,nxi),np.linspace(0,nze-1,nze)

###%%
#x = np.linspace(0,nxi-1,nxi)
#z = np.linspace(0,nze-1,nze)
#u = -1-x**2+z[:,np.newaxis]
#w = 1+x-z[:,np.newaxis]**2
#speed = np.sqrt(u*u + w*w)
#%%
f,a=getFig('Streamlines')
#a.streamplot(xi,ze,u.transpose(),w.transpose())

strm=Streamlines(xi,ze,u.transpose(),w.transpose(),maxLen=1000)
strm.plot(ax=a);axShow(a)