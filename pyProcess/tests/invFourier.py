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
n=1001;fs=(n-1)
x=np.linspace(0,2*pi,n)

y=np.sin(10*x)

f,a=getFig('Original Signal');nfig+=1
figs.append(f);axs.append(a)

axs[nfig].plot(x,y);fit(axs[nfig])

#%%
Y=scipy.fftpack.fft(y)

f,a=getFig('Signal FFT');nfig+=1
figs.append(f);axs.append(a)

ff=np.linspace(-fs/2,fs/2,n)

Y_shift=scipy.fftpack.fftshift(Y)

Y_filt=Y_shift.copy()
Y_filt[np.abs(ff)>10]=0

axs[nfig].plot(ff,np.imag(Y_filt));fit(axs[nfig])
#%%

Y_filt=scipy.fftpack.ifftshift(Y_filt)
y_filt=scipy.fftpack.ifft(Y_filt)

f,a=getFig('Filtered Signal');nfig+=1
figs.append(f);axs.append(a)

axs[nfig].plot(x,y_filt);fit(axs[nfig])

