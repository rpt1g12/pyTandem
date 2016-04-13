import numpy as np
import matplotlib.pyplot as plt
from lib.myPlots import *

#%%
path='sections/6blocks/4A15W11AoA20/cfTwCp/'
filename='SecCfCp';ext='0.dat'
#tikzpath='/home/rpt1g12/Dropbox/phd/figures/wleResults/'

#%%
fig=plt.figure()
ax=fig.add_subplot(111)
opt=6;nsec=17;direc=2;pdir=0;save=0
if (direc!=2):
    opt=4
if (opt==0): #all
    secs=range(nsec)
    ncol=3;fname='all'
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
    cf=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=[1])
    if(i==secs[0]):
        cfa=cf.copy()
    else:
        cfa+=cf
    x=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=[6,7,8])
    ax.plot(x[pdir,:],cf,label=cs)
cfa/=len(secs)
ax.plot(x[pdir,:],cfa,label='avg')
ax.set_xlabel(r'$x/L_c$')
ax.set_ylabel(r'$<C_p>$')
handle,labels=ax.get_legend_handles_labels()
legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=ncol)
ax.grid(True)
ax.figure.show()

if(save==1): 
    savePlotFile(path=path+fname+'.dat',vary=labels)