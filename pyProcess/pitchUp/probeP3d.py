# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import fcbFD
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

path='/home/'+user+'/Desktop/post/4A15W11AoA10/AoA10_20_20_10/'


files=p3d.getFileNames(path=path,pattern='solTCf+tw+Cp????.q')
nt=len(files)
time=np.zeros(nt)
probe=np.zeros((nt,4))

vnames=['Cf','twx','twy','twz','Cp'];Re=125000;M=0.4

fl=p3d.flow(path,"grid.xyz",files[0])
fl.rdHdr()
fl.rdGrid()
for n in range(nt):
    time[n]=fl.rdFlowInfo()[3]
    fl.rdSol(vnames=vnames,sfile=files[n])
    probe[n,:]=fl.blk[0].var['twz'].getValues()[80,0,12:180:48]
    
#%%
f,a=getFig()
tn=450
for i in range(4):
    a.plot(time[:tn],probe[:tn,i]*1000+i)
fit(a)