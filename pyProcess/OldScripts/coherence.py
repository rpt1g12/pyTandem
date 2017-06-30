# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from lib.myPlots import *
from lib.stats import *
from scipy.signal import coherence as msc
from scipy.signal import welch as psdw
sin20=np.sin(np.deg2rad(20))
#%%
plt.close('all')
opt=0
tprob=31;   nbase=87;   nprob=31
freq=0.19
px=np.linspace(-0.5,0.48,tprob)
#path="signalData/6blocks/4A15W11AoA20/";lze=201;lz=0.44;lwle=0.11
#path="signalData/6blocks/test/";lze=201;lz=0.44;lwle=0.11
#path="signalData/6blocks/8A00W11AoA20/";lze=201;lz=0.88;lwle=0.11
path="signalData/6blocks/4A00W11AoA20/";lze=201;lz=0.88;lwle=0.11
#path="signalData/6blocks/A00W11AoA20/";lze=101;lz=0.22;lwle=0.11

figFreq,axFreq=getFig(str(nbase)+'PSD') 

CX=np.zeros((lze,nprob))
CY=CX.copy()
for m in range(nprob):

    cnum='{:03d}'.format(m+1)
    cname="circ"+cnum+".dat"
    t0=np.loadtxt(path+cname,skiprows=1,unpack=True,usecols=[0])
    ndat=len(t0)
    x0=np.loadtxt(path+cname,skiprows=1,unpack=True,usecols=range(1,lze+1))
    
    
    
    for k in range(lze):
        out,t,nsam,fsam=rsample(x0[k,:],t0,rmAvg=True)
        if (k==0):
            p=np.zeros((lze,len(out),nprob))
        p[k,:,m]=out
        
       
    #%%

    sgnl=p[nbase,:,m];sclg='density'

    nw=4;ovlp=0.50;
    
    nseg,novlp,ntt,fmax,fmin=defWin(t,sgnl,nw,ovlp,verbose=False)

    for k in range(lze):
        f,c=msc(sgnl,p[k,:,m],fs=fsam,window='hann',nperseg=nseg,noverlap=novlp)
        fdum,adum=psdw(p[k,:,m],fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
        if (k==0):
            CC=np.zeros((lze,f.shape[0]))
            CPSD=np.zeros((lze,fdum.shape[0]))
        CC[k,:]=c.copy()
        CPSD[k,:]=adum.copy()
    z=np.linspace(-lz/2,lz/2,lze);z/=lwle;z+=(lz/lwle)/2
    st=(sin20/0.3)*f

       


    #X,Y=np.meshgrid(st,z)
    
    stdif=abs(st-freq)
    sts=np.argmin(stdif)
    CX[:,m]=CC[:,sts]
    psdmax=CPSD[:,sts].max()
    CY[:,m]=CPSD[:,sts]/psdmax
    fdum,adum=psdw(sgnl,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
    if (m==0):
        pp=np.zeros((len(adum),nprob+1))
    pp[:,m]=adum[:].copy()

#%%

top=nprob
for m in range(top):    
    axFreq.loglog(st,pp[:,m],label=str(m))

axFreq.axvline(x=st[sts],color='black',linewidth=3,linestyle='--')
axFreq.grid(b=True, which='major', color='k', linestyle='--')
axFreq.grid(b=True, which='minor', color='k', linestyle=':')
axFreq.set_xlabel(r'$St$',fontsize=20)
axFreq.set_ylabel(r'$PSD$',fontsize=20)

print(st[0],st[-1],(st[-1]-st[0])/len(st))

#%%
X,Y=np.meshgrid(px[:nprob],z)
lvl=np.linspace(0,1,41)
cm='jet'
fig,ax=getFig('MSC_St='+'{:4.2f}'.format(st[sts]))
img=ax.contourf(X,Y,CX,levels=lvl,cmap=cm) 

path2='pgfPlots/'
ax.axes.get_yaxis().set_visible(False)
ax.axes.get_xaxis().set_visible(False)
fig.savefig(path2+fig.canvas.get_window_title()+'.pdf')
#ax.axhline(y=z[nbase],color='black',linewidth=3,linestyle='--')
#%%
figPSD,axPSD=getFig('PSD_St='+'{:4.2f}'.format(st[sts]))
#psdmax=CY.max()
#CY/=psdmax
axPSD.contourf(X,Y,CY,levels=lvl,cmap='jet')
figPSD.savefig(path+figPSD.canvas.get_window_title()+'.png')
#%%

#plt.close('all')
#ax.set_visible(False)
#cbar=fig.colorbar(img,orientation='horizontal')
#cbar.ax.yaxis.set_visible(False)
#cbar.ax.xaxis.set_visible(False)
#fig.savefig(path+cm+'.pdf')