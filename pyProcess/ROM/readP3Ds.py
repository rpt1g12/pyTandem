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
path='/home/'+user+'/Desktop/post/rom/'


files=p3d.getFileNames(path=path)

vnames=['Cf','twx','twy','twz','Cp'];Re=125000;M=0.4

fl=p3d.flow(path,"grid.xyz",'solTCf+tw+Cp0000.q')
fl.rdHdr()
fl.rdGrid()
mach,aoa,Re,time=fl.rdFlowInfo(sfile='solTCf+tw+Cp1919.q')
print(mach,aoa,Re,time)
fl.rdSol(vnames=vnames)