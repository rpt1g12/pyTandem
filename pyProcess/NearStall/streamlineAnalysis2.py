import os
import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
from lib.myPlot3dOperator import *
import lib.Plot3DClasses as p3d
from scipy import stats
from scipy.interpolate import griddata
from scipy.signal import argrelextrema as locmaxmin
pi=np.pi
import getpass
user=getpass.getuser()
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')
figs=[];axs=[];hdls=[];lbls=[];lgds=[];nfig=-1;fs=18
import importlib
#%%
importlib.reload(p3d)
#%% Options
animate=True
read=False
save=True;reconstruct=False
bar=True
option=1 #1 Scale with POD projection, 2 use leastDamped modes, 3 use DMD amplitude
plotOptions=2 #1 Only POD, 2 Only DMD, 3 ALL
reducedSpectrum=False
if option==1:
    scale='POD'
elif option==2:
    scale='LD'
elif option==3:
    scale='B0'
#%%
dataset='/home/rperezt/Desktop/paraviewExports/coordinatesExport.dat'
xyz=np.loadtxt(dataset,skiprows=2,unpack=True)

f=open(dataset)
f.readline()
line=f.readline()[:-1]
nx,nt=line.split('\t')
nx,nt=int(nx),int(nt)

xyz=np.reshape(xyz,(3,nt,nx))

#%%
dt=1/32.0;mach=0.3
maxT=90;xi0=2;maxE=0.99
t=np.linspace(0,maxT*dt,maxT)
D0=np.zeros((nx,nt))
D0[:,:]=np.transpose(xyz[2,:,:])
var=[]
D=D0.copy()
D[:,:]=D0[:,:]-D0[0,0]
for i in range(xi0,nx):
    var.append(sum((D[i,:])**2)/nt)
    if var[-1]>1e-10:
        D[i,:]=D[i,:]#/(var[-1]**0.5)
temp=D[xi0:,:];del D
D=temp.copy(); del temp; nx=D.shape[0]

#%%
U,s,V=oPOD(D[:,:maxT])
rs=s.copy()
for m,_ in enumerate(s):
    rs[m]=np.sum(sn[:m+1])
    
i=1;tot=0
while tot<maxE:
    tot=sum(sn[:i])
    i+=1
rank=i
print('Rank is {:d} contains {}% of the energy'.format(rank,maxE))

ev,Phi=oDMD(D[:,:maxT])

#%%

ntmode=U.shape[1]

sn=s/np.sum(s)


plotRange=[]
freqs = np.zeros((ntmode))
grate = np.zeros((ntmode))
for i in range(ntmode):
    temp = np.log(ev[i])/(dt*mach)
    freqs[i] = temp.imag/(2*pi)
    grate[i] = temp.real
    evNorm=np.abs(ev)
    if freqs[i] >= 0  and grate[i]>0:
        plotRange.append(i)
#%%
nmodes=rank
f,a=getFig('sigma');figs.append(f);axs.append(a);nfig+=1
axs[nfig].plot(rs,lw=2)

f,a=getFig('PODmodes');figs.append(f);axs.append(a);nfig+=1
for m in range(nmodes):
    axs[nfig].plot(U[:,m],lw=2,label=str(m))

f,a=getFig('PODtime');figs.append(f);axs.append(a);nfig+=1
for m in range(nmodes):
    axs[nfig].plot(t,s[m]*V[:,m],lw=2,label=str(m))

f,a=getFig('DMDspectrum');figs.append(f);axs.append(a);nfig+=1
axs[nfig].scatter(freqs,grate,c=evNorm,s=50,cmap=plt.cm.jet)
axs[nfig].scatter(freqs[plotRange],grate[plotRange],c='yellow',s=500,cmap=plt.cm.jet,marker='*')
axs[nfig].axhline(y=0,lw=2,color='black')

f,a=getFig('DMDmodes');figs.append(f);axs.append(a);nfig+=1
for m in plotRange:
    axs[nfig].plot(Phi[:,m],lw=2,label=str(m))
    
#%%
    
if animate:
    Dr=np.dot(U[:,:rank],np.dot(np.diag(s[:rank]),np.transpose(V[:,:rank])))
    f,a=getFig('Streamline');nfig+=1
    axs.append(a);figs.append(f)
    xrange=(np.min(xyz[0,:,xi0:]),np.max(xyz[0,:,xi0:]))
    yrange=(np.min(D[:,:maxT]),np.max(D[:,:maxT]))
    a.set_xlim(xrange)
    a.set_ylim(yrange)
    for tt in range(maxT):
        a.plot(xyz[0,tt,xi0:],D[:,tt],lw=2,color='blue')
        a.plot(xyz[0,tt,xi0:],Dr[:,tt],lw=2,color='red')
        plt.pause(0.1)
        cll(a)