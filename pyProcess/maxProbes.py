import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
#%%
fTime=plt.figure()
fTime.canvas.set_window_title('Time')
axTime=fTime.add_subplot(111)

fFreq=plt.figure()
fFreq.canvas.set_window_title('Frequency')
axFreq=fFreq.add_subplot(111)
#%%
nxk=(9,6)
path='/home/rpt1g12/anaconda3/pyTandem/maxProbes/3A15W11AoA20/60-140-10/maxpln'
file=path+str(1)+'.dat'
t=np.loadtxt(file,skiprows=1,unpack=True,usecols=range(1))
nsam=len(t)
mval=np.zeros((nxk[0],nsam,5,nxk[1]))

for k in range(nxk[1]):
    for q in range(5):
        file=path+str(k+1)+'.dat'
        mval[:,:,q,k]=np.loadtxt(file,skiprows=1,unpack=True,usecols=range(q+1,nxk[0]*5+1,5))


#%%
i=0;k=3;vflag=False
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
    ff,psgn=psdw(sgn,fs=fsam,nperseg=nseg,noverlap=novlp)
    axTime.plot(tn,sgn+i*0.4)
    fit(axTime)
    
    plt.sca(axFreq)
    #cll();clt()
    axFreq.loglog(ff,psgn,label='i'+str(i))
    
handle,labels=axFreq.get_legend_handles_labels()
legend=axFreq.legend(handle,labels,bbox_to_anchor=(0,0),ncol=2,loc=3)    
axFreq.grid(b=True, which='major', color='k', linestyle='--')
axFreq.grid(b=True, which='minor', color='k', linestyle=':')
axFreq.set_xlabel(r'$f^*$',fontsize=20)
axFreq.set_ylabel(r'$PSD$',fontsize=20)
fit(axFreq)
axFreq.set_xlim(0,fmax)
axFreq.figure.canvas.draw()
