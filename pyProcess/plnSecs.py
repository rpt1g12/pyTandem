import numpy as np
import matplotlib.pyplot as plt
from lib.myPlots import *

#%%
xyz=['x','y','z'];
path='sections/6blocks/3A15W11AoA20/cfTwCp/'
o=[-0.48,0,0];n=[1,0,0];d=[0.01,0,0];oo=o.copy()
filename="plnS"
ext='0.dat'
#tikzpath='/home/rpt1g12/Dropbox/phd/figures/wleResults/'

#%%
fig=plt.figure()
ax=fig.add_subplot(111)
opt=1;nsec=6;pdir=2;ncol=1
wl=0.11
for i in range(0,nsec):
    print("sec:"+str(i))
    for j in range(3):
        oo[j]=o[j]+i*n[j]*d[j]
        if (n[j]==1):
            cs=xyz[j]+str(round(oo[j],4))
    dataset=path+filename+cs+ext
    cf=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=[4])
    x=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=[6,7,8])
    x[2]/=wl
    cfs=[b for (a,b) in sorted(zip(x[pdir],cf))]
    xs=sorted(x[pdir])
    ax.plot(xs,cfs,label=cs)
ax.set_xlabel(r'$z/\lambda_{LE}$')
ax.set_ylabel(r'$<C_p>$')
handle,labels=ax.get_legend_handles_labels()
legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=ncol)
ax.grid(True)
ax.figure.show()

#savePlotFile(path=path+fname+'.dat',vary=labels)