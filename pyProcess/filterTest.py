# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from lib.myPlots import *
from lib.stats import *
from scipy.signal import butter, filtfilt
pi=np.pi
i=0
#%%
fig=plt.figure()
fig.canvas.set_window_title('Filter Test')
ax=fig.add_subplot(111)

fq=20;fs=50
x=np.linspace(0,1,2*fs+1)
s=np.sin(2*pi*fq*x)
ax.plot(x,s,'b',linewidth=2,label='unfiltered')
#%%
for i in range(6):
    cut=0.00001+i*0.1999
    
    b,a=butter(4,cut,analog=False)
    f=filtfilt(b,a,s)
    
    ax.plot(x,f,label='filtered '+str(cut))
    handle,labels=ax.get_legend_handles_labels()
    legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=1)
    axShow(ax)
