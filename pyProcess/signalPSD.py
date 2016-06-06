import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
#%%
sim='A00vClCd'
dataset='clData/6blocks/pvClCd/'+sim+'.dat';
n,tin,clin,cdin=np.loadtxt(dataset,skiprows=1,unpack=True)
save=True;scale=False;sclg='spectrum'

tmin=0.0;ns=512
n0=np.where(tin>tmin)[0][0]
t0=tin[n0:]-tin[n0];cl0=clin[n0:];cd0=cdin[n0:]

cln,tn,nsam,fsam=rsample(cl0,t0,verbose=True,nsample=ns)
cdn,tn,nsam,fsam=rsample(cd0,t0,nsample=ns)

#%%
nw=2;ovlp=0.5;sgnl=cdn

nseg,novlp,ntt,fmax,fmin=defWin(tn,sgnl,nw,ovlp,verbose=False)
#sgnl=myFilter(sgnl,0.25/(fmax))
ff,Sxx=psdw(sgnl,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
if (scale):
    Sxx/=np.var(Sxx)
sin20=np.sin(np.deg2rad(20))
st=(sin20/0.3)*ff
#%%
fig,ax=getFig(sim,211)
ax.loglog(st,Sxx,label='Sxx')
#ax.set_xlim(0,10)
ax.grid(b=True, which='major', color='k', linestyle='--')
ax.grid(b=True, which='minor', color='k', linestyle=':')
ax.set_xlabel(r'$St$')
ax.set_ylabel(r'$PSD$')
fit(ax)
ax2=addAx(fig,212)
ax2.plot(tn,cdn,label='signal')
fit(ax2)