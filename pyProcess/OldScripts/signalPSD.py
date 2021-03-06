import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
#%%
plt.close('all')
fig,ax=getFig('Grid',111)
fig2,ax2=getFig('Signals',111)
for ll in range(1):
    #sim='G'+str(ll)
    linreg=1
    sim='8A00vClCd';prdl=1
    ttotal=42
    #sim='4A15pClCd';prdl=0
    dataset='clData/6blocks/pvClCd/new/'+sim+'.dat';
    n,tin,clin,cdin=np.loadtxt(dataset,skiprows=1,unpack=True)
#    sim='8A00vClCd';prdl=1
#    dataset='clData/6blocks/pvClCd/'+sim+'.dat';
#    n,tin,clin2,cdin2=np.loadtxt(dataset,skiprows=1,unpack=True)
#    clin+=clin2;cdin+=cdin2


    save=False;
    name='8WLE_pClCdPSD'
    
#    tmin=0.0;ns=1024
#    n0=np.where(tin>tmin)[0][0]
#    t0=tin[n0:]-tin[n0];cl0=clin[n0:];cd0=cdin[n0:]
    
    ns=1024
    tmin=tin[-1]-(ttotal/0.3)
    n0=np.where(tin>tmin)[0][0]
    t0=tin[n0:]*0.3#-tin[n0];
    if prdl==1:
        fctr=np.sqrt(1-0.3**2)
    else:
        fctr=1
    cl0=clin[n0:]*fctr;
    cd0=cdin[n0:]*fctr
    
    cln,tn,nsam,fsam=rsample(cl0,t0,verbose=True,nsample=ns)
    cdn,tn,nsam,fsam=rsample(cd0,t0,nsample=ns)
    
    clM=cln.mean()
    cdM=cdn.mean()
    if linreg==1:
        #Linear regression
        x=tn.copy();y=cln.copy()
        z=np.polyfit(x,y,3)    
        p_cl=np.poly1d(z)
    
        y=cdn.copy()
        z=np.polyfit(x,y,3)
        p_cd=np.poly1d(z)
        lin_cl=p_cl(x)
        lin_cd=p_cd(x)
    else:
        lin_cl=np.ones(len(tn))*clM
        lin_cd=np.ones(len(tn))*cdM
    
    #%%
    scale=False;sclg='density'
    nw=6;ovlp=0.5;
    sgnl=cdn-lin_cd
    
    nseg,novlp,ntt,fmax,fmin=defWin(tn,sgnl,nw,ovlp,verbose=False)
    #sgnl=myFilter(sgnl,0.25/(fmax))
    ff,Sxx=psdw(sgnl,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
    if (scale):
        Sxx/=np.var(Sxx)
    
    sin20=np.sin(np.deg2rad(20))
    st=(sin20)*ff
    #st=(sin20/0.3)*ff
    #%%

    ax.loglog(st,Sxx,label=sim)
    
    ax.grid(b=True, which='major', color='k', linestyle='--')
    ax.grid(b=True, which='minor', color='k', linestyle=':')
    ax.set_xlabel(r'$St$')
    ax.set_ylabel(r'$PSD$')
    handle,labels=ax.get_legend_handles_labels()
    legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=1)
    fit(ax)
    ax.set_xlim(0,1)

    ax2.plot(tn,sgnl,label=sim)
    handle2,labels2=ax2.get_legend_handles_labels()
    ax2.set_xlabel(r'$t^*$')
    legend2=ax2.legend(handle2,labels2,bbox_to_anchor=(1,1),ncol=1)
    fit(ax2)
    
    print('Mean='+'{:6.2e}'.format(cdn.mean()))
    print('sigma='+'{:6.2e}'.format(np.sqrt(np.var(sgnl))))
#%%
#name=sim+'PSD'

if (save==True):
    plt.sca(ax)
    path='pgfPlots/'
    savePlotFile(path=path+name+'.dat',vary=['SCd'])