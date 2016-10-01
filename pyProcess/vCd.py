# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *

sim='A00';name=sim+'vClCdPSD'
dataset='clData/6blocks/pvClCd/'
n,tin,clpin,cdpin=np.loadtxt(dataset+sim+'pClCd.dat',skiprows=1,unpack=True)
n,tin,clvin,cdvin=np.loadtxt(dataset+sim+'vClCd.dat',skiprows=1,unpack=True)

tin*=0.3
save=True;scale=False;sclg='density'
    
tmin=tmin=tin[-1]-(60.0*0.3);
ns=512
n0=np.where(tin>tmin)[0][0]
t0=tin[n0:]-tin[n0];clp0=clpin[n0:];cdp0=cdpin[n0:]
clv0=clvin[n0:];cdv0=cdvin[n0:]
    
clpn,tn,nsam,fsam=rsample(clp0,t0,verbose=True,nsample=ns)
cdpn,tn,nsam,fsam=rsample(cdp0,t0,nsample=ns)
clvn,tn,nsam,fsam=rsample(clv0,t0,nsample=ns)
cdvn,tn,nsam,fsam=rsample(cdv0,t0,nsample=ns)
    
#%%
nw=2;ovlp=0.5;sgnl=cdvn
    
nseg,novlp,ntt,fmax,fmin=defWin(tn,sgnl,nw,ovlp,verbose=False)

ff,Sxx=psdw(sgnl,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
if (scale):
    Sxx/=np.var(Sxx)
sin20=np.sin(np.deg2rad(20))
st=(sin20)*ff

#%%
plt.close('all')
fig,ax=getFig(name)
ax.loglog(st,Sxx,label='SCd')
    
ax.grid(b=True, which='major', color='k', linestyle='--')
ax.grid(b=True, which='minor', color='k', linestyle=':')
ax.set_xlabel(r'$St$')
ax.set_ylabel(r'$PSD$')
handle,labels=ax.get_legend_handles_labels()
legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=1)
fit(ax)
ax.set_xlim(0,1)


if (save):
    path='pgfPlots/'
    savePlotFile(path=path+name+'.dat',vary=labels)

