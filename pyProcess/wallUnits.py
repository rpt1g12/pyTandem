import numpy as np
import matplotlib.pyplot as plt
from lib.myPlots import *
marker = ['s', 'o', 'd', 'v', '^']
color = ['b-', 'r-', 'g-', 'm-', 'c-']

#%%
lxi=315
def indx(i,k,lxi):
    """returns index on the wall arrays"""
    l=k*(lxi+1)+i
    return l 
    
#%%
dataset='wallUnits/clean/AoA20/wplus.dat'
x,y,z,xp,yp,zp=np.loadtxt(dataset,skiprows=1,unpack=True)

#%%
fig=plt.figure()
ax=fig.add_subplot(111)
wu=yp
opt=3;skip=5
if (opt==0):
    secs=range(9)
elif (opt==1):
    secs=[0,2,4,6,8]
elif (opt==2):
    secs=[1,5]
elif (opt==3):    
    secs=[3,7]
elif (opt==4):    
    secs=[0]
elif (opt==5): 
    secs=[1,3,5,7]
j=0
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
    k=m
    ls=indx(0,k,lxi);le=indx(lxi,k,lxi)
    ax.plot(x[ls:le],wu[ls:le],color[j]+marker[j],markevery=skip,label=cs)
    j+=1
    if(j>4): j=0
ax.set_xlabel(r'$x/L_c$')
ax.set_ylabel(r'$<y^+>$')
handle,labels=ax.get_legend_handles_labels()
legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=ncol)
ax.grid(True)
ax.figure.show()
path='wallUnits/clean/AoA20/'
#plt.savefig(savepath,dpi=300)
savePlotFile(path=path+'yPlus.dat',vary=labels)