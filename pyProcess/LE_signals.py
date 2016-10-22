# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from lib.myPlots import *
import getpass
user=getpass.getuser()
pi=np.pi
sina=np.sin(np.deg2rad(20))
cosa=np.cos(np.deg2rad(20))
plt.close('all')
#%%
path='/home/'+user+'/anaconda3/pyTandem/signalData/TE_Signal/'
#%%
sim='8A00W11AoA20/';var='u'
t0,p0,u0=np.loadtxt(path+sim+'psignal.dat',skiprows=1,unpack=True)
t0-=t0[0]
#%%
f0,a0=getFig(var)
if var=='u':
    v=u0
else:
    v=p0
a0.plot(t0,v,lw=2,color='blue',label=sim)
#%%
sim='8A15W11AoA20/'
t1,p1,u1=np.loadtxt(path+sim+'psignal.dat',skiprows=1,unpack=True)
t1-=t1[0]
#%%
if var=='u':
    v=u1
    name='U_TE_signal'
else:
    v=p1
    name='p_TE_signal'
    
a0.plot(t1,v,lw=2,color='red',label=sim)

a0.set_xlabel(r'$t^*$')
a0.set_ylabel(r'$|U|$')
getLabels(ax=a0)
axShow(a0)

#%%
save=False

if save:
    savePlotFile(ax=a0,name=name)