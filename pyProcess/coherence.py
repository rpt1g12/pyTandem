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
px=np.linspace(-0.5,0.48,31)
#path="signalData/6blocks/4A15W11AoA20/";lze=201;lz=0.44;lwle=0.11
path="signalData/6blocks/A00W11AoA20/";lze=101;lz=0.22;lwle=0.11

CX=np.zeros((lze,31))
for m in range(31):

    cnum='{:03d}'.format(m+1)
    cname="circ"+cnum+".dat"
    t0=np.loadtxt(path+cname,skiprows=1,unpack=True,usecols=[0])
    ndat=len(t0)
    x0=np.loadtxt(path+cname,skiprows=1,unpack=True,usecols=range(1,lze+1))
    
    
    
    for k in range(lze):
        out,t,nsam,fsam=rsample(x0[k,:],t0,rmAvg=True)
        if (k==0):
            x=np.zeros((lze,len(out)))
        x[k,:]=out
      
    nw=2;ovlp=0.50;sgnl=x[000,:]
    
    nseg,novlp,ntt,fmax,fmin=defWin(t,sgnl,nw,ovlp,verbose=False)
    
    
    #%%
    flim=1
    for k in range(lze):
        f,c=msc(x[87,:],x[k,:],fs=fsam,window='hann',nperseg=nseg,noverlap=novlp)
        if (k==0):
            CC=np.zeros((lze,f.shape[0]))
        CC[k,:]=c
    z=np.linspace(-lz/2,lz/2,lze);z/=lwle;z+=(lz/lwle)/2
    st=(sin20/0.3)*f
    X,Y=np.meshgrid(st,z)
    
    stdif=abs(st-0.19)
    sts=np.argmin(stdif)
    CX[:,m]=CC[:,sts]
    

print(st[0],st[-1],(st[-1]-st[0])/len(st))

#%%
X,Y=np.meshgrid(px,z)
lvl=np.linspace(0,1,41)
cm='gray'
fig,ax=getFig('MSC at St_SLE='+'{:4.2f}'.format(st[sts]))
img=ax.contourf(X,Y,CX,levels=lvl,cmap=cm) 
fig.savefig(path+fig.canvas.get_window_title()+'.png')
ax.axes.get_yaxis().set_visible(False)
ax.axes.get_xaxis().set_visible(False)
#fig.savefig(path+fig.canvas.get_window_title()+'.pdf')
#%%

#plt.close('all')
#ax.set_visible(False)
#cbar=fig.colorbar(img,orientation='horizontal')
#cbar.ax.yaxis.set_visible(False)
#cbar.ax.xaxis.set_visible(False)
#fig.savefig(path+cm+'.pdf')