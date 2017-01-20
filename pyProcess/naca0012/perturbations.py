# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import fcbFD
from lib.myPlots import *
from lib.myPlot3dOperator import *
import lib.Plot3DClasses as p3d
from scipy import stats
from scipy.interpolate import griddata
pi=np.pi
import getpass
user=getpass.getuser()
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')
cosa,sina=np.cos(np.deg2rad(5)),np.sin(np.deg2rad(5))
import importlib
#%%
importlib.reload(p3d)
#%%
visible=False;save=True;rotate=True

path='/home/'+user+'/Desktop/post/nacaForcing/ss001b/'


files=p3d.getFileNames(path=path)
nt=len(files)

vnames=['r','u','v','w','e'];Re=125000;M=0.4

fl=p3d.flow(path,"grid.xyz",files[100])
fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames)
#%%
skip=15;j=15
nf=len(vnames)
xprobes=fl.blk[4].var['x'].getValues()[0::skip,j,0]
yprobes=fl.blk[4].var['y'].getValues()[0::skip,j,0]

nprobes=len(xprobes)
probes=np.zeros((nt,nprobes,nf))
#%%
f,a,im=fl.contourf(varname='u',vmax=1e-8,vmin=-1e-8,nlvl=10,bar=False)
a.set_aspect('equal')
a.set_xlim(-0.5,0.5)
a.set_ylim(-0.25,0.25)
#fl.drawMeshPlane(ax=a)
a.scatter(xprobes,yprobes,s=4,color='green')
f.tight_layout(pad=0)
#%%
for n in range(nt):
    fl.rdSol(vnames=vnames,sfile=files[n])
    for i in range(nf):
        name=vnames[i]
        probes[n,:,i]=fl.blk[4].var[name].getValues()[0::skip,j,0]
        
#%%
f2,a2=getFig('Perturbation in time')
time=np.linspace(0.005,5.0,nt)*0.4;dt=time[1]-time[0]
x0=xprobes[0]
fctr=[1.25e7,2.5e5]
lstyle=['--','-']
rprobes=range(1,nprobes)
switch=0.4
for i in rprobes:
    ii=min(1,int(i/(switch*nprobes)))
    dvdt=fcbFD(probes[:,i,4],dt)
    a2.plot(time,dvdt*fctr[ii]+i,lw=2,color='blue',ls=lstyle[ii])

a2.set_ylim(rprobes[0]-1,rprobes[-1]+1)
    