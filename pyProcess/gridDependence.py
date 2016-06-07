import numpy as np
import matplotlib.pyplot as plt
from lib.myPlots import *
from lib.stats import *
from scipy.signal import butter, filtfilt

#%%
sims=('G2','G3','G4','4A15W11AoA20')
nsecs=(201,201,161,201)
#sims=('G3','G2')
fig,ax=getFig('Grid')
for ll in range(len(sims)):
    sim=sims[ll]
    path='sections/6blocks/'+sim+'/f08/'
    filename='SecCfCp';ext='0.dat'
    
    #%% 
    opt=0;
    nsec=nsecs[ll];
    direc=2;pdir=0;
    save=1;filt=0;nsam=100
    if (direc!=2):
        opt=4
    if (opt==0): #all
        secs=range(nsec)
        ncol=2;fname='all'
    elif (opt==1): #medium
        secs=range(0,nsec,2)
        ncol=2;fname='m'
    elif (opt==2): #peak
        secs=range(1,nsec,4)
        ncol=1;fname='p'
    elif (opt==3):  #though  
        secs=range(3,nsec,4)
        ncol=1;fname='t'
    elif (opt==4):    
        secs=[0]
        ncol=1; fname='sec'+str(direc)
    elif (opt==5): #peak&trough
        secs=range(1,nsec,2)
        ncol=2;fname='p&t'
    elif (opt==6): #custom
        secs=[7]
        ncol=2;fname='custom'   
    
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
        if (filt!=0):
            f,xn,nsam,fsam=rsample(cf,x[0,:],nsample=100)
            b,a=butter(4,filt,analog=False)
            f=filtfilt(b,a,f)
        else:
            f=cf.copy();xn=x[pdir,:].copy()
        #ax.plot(xn,f,label=cs)
        if(i==secs[0]):
            cfa=f.copy()
        else:
            cfa+=f
    cfa/=len(secs)
    ax.plot(xn,cfa,label=sim+'avg',linewidth=2)
    ax.set_xlabel(r'$x/L_c$')
    ax.set_ylabel(r'$<C_p>$')
    handle,labels=ax.get_legend_handles_labels()
    legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=ncol)
    ax.grid(True)
    ax.figure.show()