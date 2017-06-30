import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
from scipy import stats
pi=np.pi
plt.close('all')
#%%
dataset='rom/4A15W11AoA10/upperSurf/POD_time.dat'
t=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=range(1))
podt=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=range(1,11))
dataset='rom/4A15W11AoA10/upperSurf/POD_SngV.dat'
m,s=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=range(2))

aoa=10+5*(1-np.cos(pi*t))
#%%
fig,ax=getFig('POD Chronos')
figs,axs=getFig('POD Energy')

for n in range(5):
    ax.plot(t,podt[n,:],linewidth=2,label='m'+str(n))

handle,labels=ax.get_legend_handles_labels()
legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=6)
fit(ax)

axs.semilogy(m[:11],s[:11],'b-o',linewidth=2)
#%%
#line=ax.axvline(x=t[0],color='green',linewidth=2,linestyle='--')
#for i in range(640):
#    adum=[t[i],t[i]]
#    line.set_xdata(adum)
#    axShow(ax)
#    plt.pause(0.01)
#    fig.savefig('/home/rpt1g12/Desktop/post/4A15W11AoA10/rom/pics/img.'+'{:04d}'.format(i)+'.png')