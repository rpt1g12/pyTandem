# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import welch as psdw
from scipy.signal import hilbert
import pywt #Wavelets
import matplotlib.pyplot as plt;plt.close('all')
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
figs=[];axs=[];hdls=[];lbls=[];lgds=[];nfig=-1;fs=18
plt.close('all')
import importlib

import scipy

#%%
mode='low'
noiseLevel=0.8
cutOff=4
width=11
n=101;fsam=(n-1)
x=np.linspace(0,2*pi,n)

y=3*np.sin(2*x)+np.sin(2*x)*np.cos(16*x)+np.cos(3.5*x)+noiseLevel*np.random.rand(n)


#%%
Y=scipy.fftpack.fft(y)

ff=np.linspace(-fsam/2,fsam/2,n)

Y_shift=scipy.fftpack.fftshift(Y)

#%%
hWidth=int(np.ceil(width/2))
fRlim=hWidth
if width%2==0:
    fLlim=hWidth
else:
    fLlim=hWidth-1
window=np.hanning(width)
windowR = np.hanning(width)[0:fRlim]
windowL = np.hanning(width)[fLlim:]
fWindow=np.ones(n)
if mode=='high':
    fWindow[np.abs(ff)<cutOff]=0
    nL=np.where(fWindow==0)[0][0]    
    nR=np.where(fWindow==0)[0][-1]+1
    fWindow[nL:nL+hWidth]=windowL
    fWindow[nR-hWidth:nR]=windowR
elif mode=='low':
    fWindow[np.abs(ff)>cutOff]=0
    nL=np.where(fWindow==1)[0][0]+1    
    nR=np.where(fWindow==1)[0][-1]
    fWindow[nR:nR+hWidth]=windowL
    fWindow[nL-hWidth:nL]=windowR
        
Y_filt=Y_shift*fWindow


#%%

Y_filt_shift=scipy.fftpack.ifftshift(Y_filt)
y_filt=scipy.fftpack.ifft(Y_filt_shift)

#%%

fctr=np.max([np.real(Y_shift).max(),np.real(Y_filt).max(),np.imag(Y_filt).max(),np.imag(Y_shift).max()])
f,a=getFig('Frequency Domain');nfig+=1
figs.append(f);axs.append(a)
axs[nfig].plot(ff,np.real(Y_shift),color='blue',lw=2,label='Re{y}');
axs[nfig].plot(ff,np.real(Y_filt),color='cyan',lw=2,linestyle='--',label='Re{y_filt}');
axs[nfig].plot(ff,np.imag(Y_shift),color='red',lw=2,label='Im{y}');
axs[nfig].plot(ff,np.imag(Y_filt),color='orange',lw=2,linestyle='--',label='Im{y_filt}');
axs[nfig].plot(ff,fWindow*fctr,color='green',lw=2,linestyle='-.',label='window');


f,a=getFig('Window');nfig+=1
figs.append(f);axs.append(a)
axs[nfig].plot(window,lw=2,color='black',label='full')
axs[nfig].plot(range(fRlim),windowR,lw=2,color='cyan',linestyle='--',label='right')
axs[nfig].plot(range(fLlim,width),windowL,lw=2,color='magenta',linestyle='--',label='left')

f,a=getFig('Time Domain');nfig+=1
figs.append(f);axs.append(a)

y_filt2=fourierFilter(y,mode=mode,fsam=fsam,cutOff=cutOff,width=2*width)

axs[nfig].plot(x,y,color='blue',lw=2,label='y')
axs[nfig].plot(x,np.real(y_filt),color='red',lw=2,linestyle='--',label='y_filt');
axs[nfig].plot(x,np.real(y_filt2),color='cyan',lw=2,linestyle='-.',label='y_filt2');


#%%
nfigs=len(figs)
for i in range(nfigs):
    fit(axs[i])
    hdl,lbl,lgd=getLabels(ax=axs[i],ncol=3,fontsize=15)
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)