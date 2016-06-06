import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
#from lib.matplotlib2tikz import save as tikz_save
#plt.close('all')
#%%
osim=4
if(osim==0):
    sim='A00'
elif(osim==2):
    sim='A15'
elif(osim==3):
    sim='3A15'
elif(osim==4):
    sim='A4A15'
sim+='W11AoA20'
dataset='clData/6blocks/'+sim+'.dat';
n,tin,clin,cdin=np.loadtxt(dataset,skiprows=1,unpack=True)
save=False;scale=False;sclg='spectrum'
#%%
#fOriginal=plt.figure()
#fOriginal.canvas.set_window_title('Original')
#axOriginal=fOriginal.add_subplot(111)

fTime=plt.figure()
fTime.canvas.set_window_title('Time '+sim)
axTime=fTime.add_subplot(111)

fFreq=plt.figure()
fFreq.canvas.set_window_title('Frequency '+sim)
axFreq=fFreq.add_subplot(111)


#%%
show=0
plt.sca(axTime)
axTime.lines.clear()
cll();clt();ns=512
clh=0.72;cdh=0.29;tmin=tin[-1]-60.0
n0=np.where(tin>tmin)[0][0]
t0=tin[n0:]-tin[n0];cl0=clin[n0:];cd0=cdin[n0:]
cln,tn,nsam,fsam=rsample(cl0,t0,verbose=True,nsample=ns)
cdn,tn,nsam,fsam=rsample(cd0,t0,nsample=ns)
s=np.sqrt(np.var(cln))
clM=cln.mean()
cds=np.sqrt(np.var(cdn))
cdM=cdn.mean()
if (show==1):
    axTime.plot(tin,clin,'b:');lClin=plt.gca().lines[-1]
axTime.plot(tn,cln,'b',label='Cl');lCln=plt.gca().lines[-1]
axTime.plot(tn,cdn,'b',label='Cd');lCln=plt.gca().lines[-1]
axTime.plot([tn[0],tn[-1]],[clh,clh],'k',linewidth=2,label='exp');lClh=plt.gca().lines[-1]
axTime.plot([tn[0],tn[-1]],[clM,clM],'g-.',linewidth=2,label='avg');lClnm=plt.gca().lines[-1]
axTime.plot([tn[0],tn[-1]],[clM-s,cln.mean()-s],'r-.',linewidth=2,label='s0');lCls0=plt.gca().lines[-1]
axTime.plot([tn[0],tn[-1]],[clM+s,clM+s],'r-.',linewidth=2,label='s1');lCls1=plt.gca().lines[-1]
axTime.plot([tn[0],tn[-1]],[cdh,cdh],'k',linewidth=2,label='exp');lCdh=plt.gca().lines[-1]
axTime.plot([tn[0],tn[-1]],[cdM,cdM],'g-.',linewidth=2,label='avg');lCdnm=plt.gca().lines[-1]
axTime.plot([tn[0],tn[-1]],[cdM-cds,cdn.mean()-cds],'r-.',linewidth=2,label='s0');lCds0=plt.gca().lines[-1]
axTime.plot([tn[0],tn[-1]],[cdM+cds,cdM+cds],'r-.',linewidth=2,label='s1');lCds1=plt.gca().lines[-1]
axTime.set_xlabel(r'$t^*$',fontsize=20)
axTime.set_ylabel(r'$C_l$',fontsize=20)
handle,labels=axTime.get_legend_handles_labels()
legend=axTime.legend(handle,labels,bbox_to_anchor=(1,1),ncol=5)
scl=r'$\sigma_L=$'+str(np.around(s,3))+'\t'+r'$\bar{C_L}=$'+str(np.around(clM,3))
scd=r'$\sigma_D=$'+str(np.around(cds,3))+'\t'+r'$\bar{C_D}=$'+str(np.around(cdM,3))
axTime.text(0.25,1.0,scl,ha='center',va='bottom',transform=axTime.transAxes)
axTime.text(0.75,1.0,scd,ha='center',va='bottom',transform=axTime.transAxes)
fit(axTime)

#%%

nw=2;ovlp=0.5;sgnl=cdn

nseg,novlp,ntt,fmax,fmin=defWin(tn,sgnl,nw,ovlp,verbose=False)
#sgnl=myFilter(sgnl,0.25/(fmax))
ff,pcl=psdw(sgnl,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
if (scale):
    pcl/=np.var(cln)
#%%
plt.sca(axFreq)
axFreq.lines.clear()
#cll();clt()
sin20=np.sin(np.deg2rad(20))
st=(sin20/0.3)*ff
axFreq.loglog(st,pcl,label='SCl')

axFreq.grid(b=True, which='major', color='k', linestyle='--')
axFreq.grid(b=True, which='minor', color='k', linestyle=':')
axFreq.set_xlabel(r'$St$',fontsize=20)
axFreq.set_ylabel(r'$PSD$',fontsize=20)
fit(axFreq)
axFreq.set_xlim(0,2)
#%%
name=sim.split('W')[0]+'SCd'
if (save==True):
    plt.sca(axFreq)
    path='pgfPlots/'
    savePlotFile(path=path+name+'.dat',vary=['SCd'],varx=['f'])
