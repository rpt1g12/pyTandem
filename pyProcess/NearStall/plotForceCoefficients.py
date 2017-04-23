import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
from scipy import stats
import getpass
user=getpass.getuser()
pi=np.pi
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')
figs=[];axs=[];hdls=[];lbls=[];lgds=[];nfig=-1;fs=18
#%%
folder='clData/aoa10/'
sim="4A15W11AoA10"
dataset='/home/'+user+'/anaconda3/pyTandem/'+folder+sim+'.dat';
ts,te=0,315
n,t0,clin,cdin,taoa,tmach=np.loadtxt(dataset,skiprows=1,unpack=True)
#%%
ns=np.where((t0>=ts))[0][0]
ne=np.where((t0<=te))[0][-1]
t=t0[ns:ne]
cl=clin[ns:ne];cd=cdin[ns:ne]
aoa=taoa[ns:ne]
#%%
f,a=getFig('t_cl');figs.append(f),axs.append(a);nfig+=1
axs[nfig].plot(t,cl,lw=2,color='blue',label=r'$C_L$')
axs[nfig].plot(t,cd,lw=2,color='red',label=r'$C_D$')
axs[nfig].axvline(x=t[0],color='black',linewidth=2,linestyle='--') 
axs[nfig].set_xlabel(r'$t^*$',fontsize=fs)
axs[nfig].set_ylabel(r'$C_L\, &\, C_D$',fontsize=fs)
#%%
f,a=getFig('t_aoa');figs.append(f),axs.append(a);nfig+=1
axs[nfig].plot(t,aoa,lw=2,color='blue',label=r'$\alpha$')
#%%
f,a=getFig('aoa_cl');figs.append(f),axs.append(a);nfig+=1
axs[nfig].plot(aoa,cl,lw=2,color='blue',label=r'$C_L$')
axs[nfig].plot(aoa,cd,lw=2,color='red',label=r'$C_D$')
#%%
nfigs=len(figs)
for i in range(nfigs):
    fit(axs[i])
    hdl,lbl,lgd=getLabels(ax=axs[i],ncol=2,fontsize=15)
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
    