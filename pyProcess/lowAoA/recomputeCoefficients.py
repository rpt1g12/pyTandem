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
folder='clData/lowAoA/'
sim="1A15W11AoA00"
dataset='/home/'+user+'/anaconda3/pyTandem/'+folder+sim+'.dat';
ts,te=0,315
n,tin,clin,cdin,aoain,tmach=np.loadtxt(dataset,skiprows=1,unpack=True)
#%%
aoa0=0;
sina,cosa=(np.sin(np.deg2rad(aoa0)),np.cos(np.deg2rad(aoa0)))

cl,cd=(clin*cosa+cdin*sina,cdin*cosa-clin*sina)
#%%
f,a=getFig('ClCd');figs.append(f),axs.append(a);nfig+=1

axs[nfig].plot(tin,cl,label='cl')
axs[nfig].plot(tin,cd,label='cd')

hdl,lbl,lgd=getLabels(axs[nfig],ncol=2,fontsize=15)
hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)

#savePlotFile(path=folder,ax=axs[nfig],vary=lbl)
#%%
f,a=getFig('ClCd_pv');figs.append(f),axs.append(a);nfig+=1
folder='clData/lowAoA/pvClCd/'
sim="1A15W11AoA00"
dataset='/home/'+user+'/anaconda3/pyTandem/'+folder+sim+'.dat';
n,tin,clpin,clvin,cdpin,cdvin=np.loadtxt(dataset,skiprows=1,unpack=True)

aoa0=0
sina,cosa=(np.sin(np.deg2rad(aoa0)),np.cos(np.deg2rad(aoa0)))

clp,clv,cdp,cdv=(clpin*cosa+cdpin*sina,clvin*cosa+cdvin*sina,cdpin*cosa-clpin*sina,cdvin*cosa-clvin*sina)

axs[nfig].plot(tin,clp,label='Clp')
axs[nfig].plot(tin,cdp,label='Cdp')
axs[nfig].plot(tin,clv,label='Clv')
axs[nfig].plot(tin,cdv,label='Cdv')

#hdl,lbl,lgd=getLabels(axs[nfig],ncol=4,fontsize=15)
#hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
#
#savePlotFile(path=folder,ax=axs[nfig],vary=lbl)

axs[nfig].plot(tin,clp+clv,label=r'$C_l$',linestyle='--')
axs[nfig].plot(tin,cdp+cdv,label=r'$C_d$',linestyle='--')
hdl,lbl,lgd=getLabels(axs[nfig],ncol=6,fontsize=15)
hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
fit(axs[nfig],(0,0.2))