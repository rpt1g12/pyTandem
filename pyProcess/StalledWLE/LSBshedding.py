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
#folder='/media/{:}/082F-63FE/phd/thesis/data/StalledWLE/fig25_3WLE/'
folder='/media/{:}/082F-63FE/phd/thesis/data/StalledWLE/fig25/'.format(user)
#folder='/media/{:}/082F-63FE/phd/thesis/data/StalledWLE/fig24/'
file='4WLE_T4'
path=folder+file+'.dat'
columns=list(range(1,14,2))
st=np.loadtxt(path,skiprows=1,unpack=True,usecols=[0])
psd=np.loadtxt(path,skiprows=1,unpack=True,usecols=columns)
nprobe,n=psd.shape
#%%
#rescale frequency
M=0.3;sin20=np.sin(np.deg2rad(20))
#4WLE_T2
#Ue=0.4685246100212967
#th=0.000695462637932316
#4WLE_T4
#Ue=0.4020089494807683
#th=0.0007055318336547253
Ue=0.3;th=1
##3WLE_T1
#Ue=0.4852434254676617
#th=0.0007028910971796682
##3WLE_T2
#Ue=0.4856653212420562
#th=0.0007012175418189963
##3WLE_T3
#Ue=0.41372536072052213
#th=0.0006990588421679554
#2WLE_T1
#Ue=0.4833127562549126
#th=0.0007049064446625458
##2WLE_T2
#Ue=0.4294783186914878
#th=0.0006931799226859244
##2SLE_T1
#Ue=0.4194793168266415
#th=0.0006522503825071501
##2SLE_T2
#Ue=0.419489160747591
#th=0.0006522068651498401
#%%
st*=(M/sin20)*th/Ue
f,a=getFig('BLscaled_'+file)
varx=['St{:1d}'.format(i) for i in range(nprobe)]
vary=['i{:1d}'.format(i) for i in range(nprobe)]
for i in range(nprobe):
    a.loglog(st,psd[i,:])
#%%
if save==True:
    savePlotFile(path=folder,ax=a,varx=varx,vary=vary)