# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *

sim='4A15';save=True
dataset='clData/6blocks/pvClCd/'
n,t,clp,cdp=np.loadtxt(dataset+sim+'pClCd.dat',skiprows=1,unpack=True)
n,t,clv,cdv=np.loadtxt(dataset+sim+'vClCd.dat',skiprows=1,unpack=True)
fig,ax=getFig(sim+'pvCd.dat')
t-=t[0]
ax.plot(t,cdp,label='cdp')
ax.plot(t,cdv,label='cdv')
handle,labels=ax.get_legend_handles_labels()

if (save):
    path='pgfPlots/'
    name=fig.canvas.get_window_title()
    savePlotFile(path=path+name,vary=labels)

