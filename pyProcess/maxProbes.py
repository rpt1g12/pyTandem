import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *

#%%
plt.close('all')
save=True; scale=False; step=False;sclg='density'
nxk=(9,8)
sim='A4A15W11AoA20'
path='maxProbes/'+sim+'/60-140-10/maxpln'
file=path+str(1)+'.dat'
t=np.loadtxt(file,skiprows=1,unpack=True,usecols=range(1))
nsam=len(t)
t*=0.3
mval=np.zeros((nxk[0],nsam,5,nxk[1]))

for k in range(nxk[1]):
    for q in range(5):
        file=path+str(k+1)+'.dat'
        mval[:,:,q,k]=np.loadtxt(file,skiprows=1,unpack=True,usecols=range(q+1,nxk[0]*5+1,5))

#%%
fTime=plt.figure()
fTime.canvas.set_window_title('Time '+sim)
axTime=fTime.add_subplot(111)

fFreq=plt.figure()
fFreq.canvas.set_window_title('Frequency '+sim)
axFreq=fFreq.add_subplot(111)
#%%
i=0;k=7;vflag=False
q=4;nw=16;ovlp=0.5
sgn=mval[i,:,q,k]
sgn,tn,nsam,fsam=rsample(sgn,t)
nseg,novlp,ntt,fmax,fmin=defWin(tn,sgn,nw,ovlp,verbose=False)

axTime.lines.clear()
axFreq.lines.clear()
for i in range(0,nxk[0]-2,1):
    if(vflag):
        sgn=(mval[i,:,1,k]**2+mval[i,:,2,k]**2+mval[i,:,3,k]**2)**0.5
    else:
        sgn=mval[i,:,q,k]
    sgn,tn,nsam,fsam=rsample(sgn,t,verbose=True,rmAvg=True)
    ff,psgn=psdw(sgn,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
    axTime.plot(tn,sgn+i*0.4)
    fit(axTime)
    
    plt.sca(axFreq)
    #cll();clt()
    if (step):
        psgn*=10**i
        if (scale):
            var=np.var(sgn)
            psgn/=var
    
    sin20=np.sin(np.deg2rad(20))
    #st=(sin20/0.3)*ff
    st=sin20*ff
    axFreq.loglog(st,psgn,label='i'+str(i),linewidth=2)
    
handle,labels=axFreq.get_legend_handles_labels()
legend=axFreq.legend(handle,labels,bbox_to_anchor=(0,0),ncol=2,loc=3)    
axFreq.grid(b=True, which='major', color='k', linestyle='--')
axFreq.grid(b=True, which='minor', color='k', linestyle=':')
axFreq.set_xlabel(r'$St$',fontsize=20)
axFreq.set_ylabel(r'$PSD$',fontsize=20)
fit(axFreq)
axFreq.set_xlim(0,fmax)
axFreq.figure.canvas.draw()

#%%
name=sim.split('W')[0]+'pS'+'k'+str(k)
if (save==True):
    plt.sca(axFreq)
    path='pgfPlots/'
    savePlotFile(path=path+name+'.dat',vary=labels,varx=['St'+str(i) for i in range(len(labels))])