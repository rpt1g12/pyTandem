# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from lib.myPlots import *
pi=np.pi
sina=np.sin(np.deg2rad(20))
cosa=np.cos(np.deg2rad(20))
#plt.close('all')
#%%
path='/home/rperezt/anaconda3/pyTandem/signalData/TE_Signal/'
#%%
sim='8A00W11AoA20/'
t0,p0,u0=np.loadtxt(path+sim+'psignal.dat',skiprows=1,unpack=True)
t0-=t0[0]
#%%
f0,a0=getFig(sim)
a0.plot(t0,u0,lw=2,color='blue',label=sim)
#%%
sim='8A15W11AoA20/'
t1,p1,u1=np.loadtxt(path+sim+'psignal.dat',skiprows=1,unpack=True)
t1-=t1[0]
#%%
a0.plot(t1,u1,lw=2,color='red',label=sim)
axShow(a0)