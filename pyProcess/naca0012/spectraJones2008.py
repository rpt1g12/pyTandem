# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from lib.myPlots import *
from scipy.signal import welch as psdw
from lib.stats import *
pi=np.pi

save=False
#%%
fJ,SkJ=np.loadtxt('spectralData/fig7b_3.dat',skiprows=1,unpack=True)
u,v,w,t0=np.loadtxt('spectralData/signal.dat',skiprows=1,unpack=True)
k0=0.5*(u**2+v**2+w**2)
k1,t,nsam,fsam=rsample(k0,t0*0.4,verbose=True)
u,v,w,t1=np.loadtxt('spectralData/signal3.dat',skiprows=1,unpack=True)
k0=0.5*(u**2+v**2+w**2)
k2,t2,nsam,fsam=rsample(k0,t1*0.4,verbose=True)


#%%
scale=False;sclg='density'
nw=6;ovlp=0.5;
    
nseg,novlp,ntt,fmax,fmin=defWin(t1,k1,nw,ovlp,verbose=False)
f,Sk=psdw(k1,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
nseg,novlp,ntt,fmax,fmin=defWin(t2,k2,nw,ovlp,verbose=False)
f2,Sk2=psdw(k2,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)

#%%
plt.close('all')
fig,a=getFig('spectraComparison')

a.semilogy(f,Sk,label='LES',lw=2)
a.semilogy(f2,Sk2,label='LES2',lw=2)
a.semilogy(fJ,SkJ,label='Jones2008',lw=2)
a.set_xlim(0,150)
if save:
    savePlotFile(ax=a)
