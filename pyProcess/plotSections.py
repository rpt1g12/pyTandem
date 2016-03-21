import numpy as np
import matplotlib.pyplot as plt
from lib.myPlots import *

#%%
path='sections/6blocks/A00W11AoA20/cfTwCp/'
filename='cfSec';ext='0.dat'
#tikzpath='/home/rpt1g12/Dropbox/phd/figures/wleResults/'

#%%
fig=plt.figure()
ax=fig.add_subplot(111)
opt=0
if (opt==0): #all
    secs=range(9)
    ncol=3
elif (opt==1): #medium
    secs=[0,2,4,6,8]
    ncol=2
elif (opt==2): #peak
    secs=[1,5]
    ncol=1
elif (opt==3):  #though  
    secs=[3,7]
    ncol=1
elif (opt==4):    
    secs=[0]
    ncol=1
elif (opt==5): #peak&trough
    secs=[1,3,5,7]
    ncol=2
for i in secs:
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
    dataset=path+filename+cs+ext
    cf=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=[4])
    if(i==secs[0]):
        cfa=cf.copy()
    else:
        cfa+=cf
    x=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=[6,7,8])
    ax.plot(x[0,:],cf,label=cs)
cfa/=len(secs)
ax.plot(x[0,:],cfa,label='avg')
ax.set_xlabel(r'$x/L_c$')
ax.set_ylabel(r'$<C_f>$')
handle,labels=ax.get_legend_handles_labels()
legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=ncol)
ax.grid(True)
ax.figure.show()

savePlotFile(path=path+'allCf.dat',vary=labels)