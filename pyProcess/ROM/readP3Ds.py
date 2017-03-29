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
path='/home/'+user+'/Desktop/post/1A15W11AoA10/ss002/'
path_ss='/home/'+user+'/Desktop/post/rom/fl_ss/'

files=p3d.getFileNames(path=path)

vnames=['r','u','v','w','p']#['r','dpdz','drdz','w','p']

fl=p3d.flow(path,"grid.xyz",files[0])
fl.rdHdr()
fl.rdGrid()
mach,aoa,Re,time=fl.rdFlowInfo(sfile=files[0])
print(mach,aoa,Re,time)

#%%
rx,ry,rz=range(90),range(1),range(fl.blk[0].size[2])
lxib,letb,lzeb=[len(rx)],[len(ry)],[len(rz)]
fl_ss=p3d.flow(path_ss,"grid.xyz",files[0])
fl_ss.setHdr([1,1,1],lxib,letb,lzeb)

for n in range(len(files)):
    fl.rdSol(sfile=files[n],vnames=vnames)
    ss=fl.blk[1].getSubset(xlim=rx,ylim=ry,zlim=rz)
    ss.getMetrics()
    ss.derive('p','ze')
    ss.derive('r','ze')    
    fl_ss.blk[0]=ss
    if n==0:
        fl_ss.wrGrid()
    fl_ss.wrSol(files[n],vnames=['r','dpdze','drdze','w','p'])