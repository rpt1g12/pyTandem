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
save=False
path='/home/rperezt/Desktop/post/1A15W11AoA15.8/ss001/'

files=p3d.getFileNames(path=path)

vnames=['r','u','v','w','p'];

fl=p3d.flow(path,"grid.xyz",files[0])
fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames)
flInfo=fl.rdFlowInfo()

flavg=fl.getAvg(files,vnames)

f1,a1,im1=fl.contourf(varname='u',vmin=-0.3,vmax=0.6,k=0,nlvl=256,cmap=plt.cm.hot)
a1.set_aspect('equal')

f,a,im=flavg.contourf(varname='u',vmin=-0.3,vmax=0.6,k=0,nlvl=256,cmap=plt.cm.hot)
a.set_aspect('equal')