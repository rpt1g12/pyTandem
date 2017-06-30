import numpy as np
import matplotlib.pyplot as plt
from lib.myPlots import *
from lib.stats import *
from scipy.signal import butter, filtfilt

#%%
save=0
sims=('G1','G2','G3','G4')
nsecs=(201,201,201,161)
#sims=('G3','G2')
fig,ax=getFig('Grid')
fig2,ax2=getFig('Diff')
for ll in range(len(sims)):
    sim=sims[ll]
    path='sections/6blocks/gridDependence/'+sim+'/f08/'
    filename='SecCfCp';ext='0.dat'
    
    #%% 
    opt=0;
    nsec=nsecs[ll];
    direc=2;
    nsam=512
    secs=range(nsec)
    ncol=2;fname='all'
  
    
    for i in secs:
        if (direc==2):
            m=int(round(12.5*i))
            if (i%2==0):
                countm=i//2+1
                cs="M"
                cs+=str(countm)
            else:
                if((i//2)%2==0):
                    countp=i//4+1
                    cs="P"
                    cs+=str(countp)
                else:
                    countt=i//4+1
                    cs="T"
                    cs+=str(countt)
        else:
            cs=''
        dataset=path+filename+cs+ext
        cf=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=[0])
    
        x=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=[6,7,8])
       
        f=cf.copy();xn=x[0,:].copy()
        if(i==secs[0]):
            cfa=f.copy()
        else:
            cfa+=f
    cfa/=len(secs)
    ax.plot(xn,cfa,label=sim,linewidth=2)
    if(ll==0):
        dif=np.zeros((len(sims),nsam))
        dif[0],xdif=rsample(cfa,xn,nsam,force=True)[0:2]
        ax2.plot(xdif,dif[ll]-dif[0],label=sim,linewidth=2)
    else:
        dif[ll]=rsample(cfa,xn,nsam,force=True)[0]
        dif[ll]=(dif[ll]-dif[0])/dif[0]
        ax2.plot(xdif,dif[ll],label=sim,linewidth=2)
    
#%%
ax.set_xlabel(r'$x/L_c$',fontsize=20)
ax.set_ylabel(r'$<C_f>$',fontsize=20)
ax2.set_xlabel(r'$x/L_c$',fontsize=20)
ax2.set_ylabel(r'$\frac{<C_{f}>_{Gi}-<C_{f}>_{G1}}{<C_{f}>_{G1}}$',fontsize=20)
        
handle,labels=ax.get_legend_handles_labels()
legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=ncol)
ax.grid(True)
#ax.invert_yaxis()
ax.set_ylim(-0.005,0.005)
axShow(ax)

handle2,labels2=ax2.get_legend_handles_labels()
legend2=ax2.legend(handle2,labels2,bbox_to_anchor=(1,1),ncol=ncol)
ax2.grid(True)
ax2.set_ylim(-1,1)

axShow(ax2)

#%%
path='pgfPlots/'
fname='gridCp'
plt.sca(ax)
if(save==1): 
    savePlotFile(path=path+fname+'.dat',vary=labels)