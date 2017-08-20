import numpy as np
from scipy.signal import spectrogram 
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
from scipy import stats
import getpass
from matplotlib import colors,ticker
user=getpass.getuser()
pi=np.pi
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')
figs=[];axs=[];hdls=[];lbls=[];lgds=[];nfig=-1;fs=18
#%%
vmin=1e-5;vmax=1e-1
nsample=2**10
nw=32;ovlp=0.80;sclg='density'
save=True
prandtl=True #Apply Prandtl-Glauert correction?
AoA=10 #Angle of attack
A=0 #WLE amplitude
nwave=8 #Number of WLE wavelengths
if A>0:
    wavy=True
else:
    wavy=False
ts,te='min','max' #Time interval considered
#%%
folder='clData/aoa10/'
sim="{:1d}A{:02d}W11AoA{:02d}".format(nwave,A,AoA)
dataset='/home/'+user+'/anaconda3/pyTandem/'+folder+sim+'.dat';
spath='/home/rpt1g12/Documents/thesis/data/nearStall/clcdHistory/'
n,t0,clin,cdin,taoa,tmach=np.loadtxt(dataset,skiprows=1,unpack=True)
if ts=='min': ts=min(t0) 
if te=='max': te=max(t0) 
M_inf=tmach[-1]
if prandtl:
    beta=np.sqrt(1-M_inf**2)
else:
    beta=1
#%%
ns=np.where((t0>=ts))[0][0]
ne=np.where((t0<=te))[0][-1]
cl,t,nsam,fsam=rsample(clin[ns:ne]*beta,t0[ns:ne],verbose=False,force=True,nsample=nsample)
cd,t,nsam,fsam=rsample(cdin[ns:ne]*beta,t0[ns:ne],verbose=False,force=True,nsample=nsample)
aoa,t,nsam,fsam=rsample(taoa[ns:ne],t0[ns:ne],verbose=True,force=True,nsample=nsample)

var=cl
var_filt=fourierFilter(var,'high',fsam,0.07,15,True)
t-=t[0]
#%%
if wavy:
    fname='{:1d}WLE{:02d}'.format(nwave,AoA)
else:
    fname='{:1d}SLE{:02d}'.format(nwave,AoA)
f,a=getFig('t_cl_'+fname);figs.append(f),axs.append(a);nfig+=1
axs[nfig].plot(t,var,lw=2,color='blue',label='cl_full')
axs[nfig].plot(t,var_filt,lw=2,color='red',label='cl_filt')
axs[nfig].plot(t,var-var_filt,lw=2,color='cyan',label='cl-cl_filt')
axs[nfig].axvline(x=t[0],color='black',linewidth=2,linestyle='--') 
axs[nfig].set_xlabel(r'$t^*$',fontsize=fs)
axs[nfig].set_ylabel(r'$C_L\, &\, C_D$',fontsize=fs)
#%% Define windows

sgnl=var_filt
nseg,novlp,ntt,fmax,fmin=defWin(t,sgnl,nw,ovlp,verbose=True)

#%% Perform STFT
ff2,tt,psd=spectrogram(sgnl,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)


#%% Plot Spectogram
TT,FF=np.meshgrid(tt,ff2)
f,a=getFig('Spectogram'+fname);figs.append(f),axs.append(a);nfig+=1
im=a.pcolor(TT,FF,psd,norm=colors.LogNorm(vmin=psd.min(), vmax=psd.max()),cmap='jet',vmin=vmin,vmax=vmax)
cb=f.colorbar(im,orientation='vertical')
cb.set_label('PSD',fontsize=fs)
a.set_ylabel('f',fontsize=fs)
a.set_xlabel('t',fontsize=fs)
a.set_ylim(fmin,fmax)
a.set_yscale('log')
#%%
nfigs=len(figs)-1
for i in range(nfigs):
    fit(axs[i])
    hdl,lbl,lgd=getLabels(ax=axs[i],ncol=3,fontsize=15)
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
