import numpy as np
import matplotlib.pyplot as plt
from lib.myPlots import *

#%%
#plt.close('all')
path='aerofoils/'
filename='naca0020';ext='.dat'
#tikzpath='/home/rpt1g12/Dropbox/phd/figures/wleResults/'
dataset=path+filename+ext
pname='A15'
pext='Probes.dat'
nprob=2;
x=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=[0,1,2])
dataset=path+pname+pext
for i in range(nprob):
    ii=i+1
    var=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=range(3*i,3*ii))   
    l=len(var[0])    
    if (i==0):
        xp=np.zeros((l,2,nprob))
    xp[:,0,i]=var[0];xp[:,1,i]=var[1]



#%%
fig=plt.figure()
fig.canvas.set_window_title(pname+' Probes')
ax=fig.add_subplot(111)

#ax.plot(x[0],-x[1],label='BtSide')
#ax.plot(x[0],x[1],label='UpSide')
for i in range(nprob):
    ax.plot(xp[:,0,i],xp[:,1,i],'kx',label='Prob'+str(i))
ax.set_aspect('equal')
handle,labels=ax.get_legend_handles_labels() 
savePlotFile(path='aerofoils/'+pname+'PFig.dat',vary=labels)
