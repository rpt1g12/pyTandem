# -*- coding: utf-8 -*-
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

import importlib
#%%
importlib.reload(p3d)

#%%
comp=True;block=4;kk=62
#path='/home/rperezt/Desktop/post/1A15W11AoA15.8/ss001/'
path='/home/rperezt/Desktop/post/1A15W11AoA15.8/fullDomain/'
#path='/media/rperezt/082F-63FE/post/naca0012G1/ss003/'
path='/home/rperezt/Desktop/post/4A15W11AoA20/vsmallDomain/'
files=p3d.getFileNames(path=path)

vnames=['r','u','v','w','p'];
#vnames2=['Q','Wx','Wy','Wz','Cp'];

fl=p3d.flow(path,"grid.xyz",'solTA.qa')


fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames)
#fl.rdSol(vnames2,'solTQ+W+Cp.qa')
flInfo=fl.rdFlowInfo()
cosa,sina=np.cos(np.deg2rad(flInfo[1])),np.sin(np.deg2rad(flInfo[1]))
up=fl.blk[block]
del fl
#%%

up.getMetrics()

nx,ny,xo,yo=(up.mets[3].getValues()[:,0,kk],up.mets[4].getValues()[:,0,kk],
            up.var['x'].getValues()[:,0,kk],up.var['y'].getValues()[:,0,kk])
#%%
Wz=up.derive('v','x',True)-up.derive('u','y',True)

up.setData('Wz',Wz)

f,a,im=up.contourf(varname='Wz',vmin=-10,vmax=10,k=kk,nlvl=50,cmap=plt.cm.bwr,bar=False)
a.set_aspect('equal')
#%%
f2,a2=getFig('Profiles')
rakes=[];
dstr=[];theta=[];xbl=[];ybl=[];delta=[]
nn=101;l=0.008;dl=l/(nn-1)
rl=np.array([i*dl for i in range(nn)])
u=np.zeros_like(rl)
dudy=np.zeros(nn)
i=4
while i<60:#dudy[0]>=0:
    print('i={:}'.format(i))
    rakes.append(p3d.rake(xo[i],yo[i],nx[i],ny[i],nn,l))
    x,y=rakes[-1].x,rakes[-1].y  
    tx,ty=rakes[-1].tx,rakes[-1].ty
    wz=up.interpolate2dk('Wz',x,y,kk,mode='values')[:,0,0]
    u=-myIntegral(wz,rl);u[np.isnan(u)]=max(u)
    if comp:
        rr=up.interpolate2dk('r',x,y,kk,mode='values')[:,0,0]
    rakes[-1].var=u
    dudy=fcbFD(rakes[-1].var,dl)
    if True:
        j99=np.where(u<0.99*u[-1])[0][-1]+1;jmax=-1
        delta.append(j99*dl)
        um=u[jmax]
        rakes[-1].var/=um
        tmp=rakes[-1].var[:j99]
        tmp2=rl[:j99]
        a2.plot(tmp,tmp2/tmp2[-1])
        a.plot(x,y,lw=1,color='orange') 
        if comp:        
            rm=rr[jmax]
            rr/=rm
            tmp=1-(rr[:jmax])*rakes[-1].var[:jmax]
            tmp2=rl[:jmax]
            dstr.append(np.trapz(tmp,tmp2))
            tmp=((rr[:jmax])*rakes[-1].var[:jmax])*(1-rakes[-1].var[:jmax])
            theta.append(np.trapz(tmp,tmp2))  
        else:
            tmp=1-rakes[-1].var[:jmax]
            tmp2=rl[:jmax]
            dstr.append(np.trapz(tmp,tmp2))
            tmp=rakes[-1].var[:jmax]*(1-rakes[-1].var[:jmax])
            theta.append(np.trapz(tmp,tmp2))
        xbl.append(rakes[-1].x[j99])
        ybl.append(rakes[-1].y[j99])
    i+=1

#rakes.pop()
dudy=fcbFD(rakes[-1].var,dl)
theta=np.asarray(theta)
dstr=np.asarray(dstr)
delta=np.asarray(delta)
H=dstr/theta
nr=len(rakes)

cout='theta={:},Ue={:},theta/Ue={:}'.format(theta[-1],um,theta[-1]/um)
print(cout)   

#%%  
fH,aH=getFig('Shape Factor')
aH.plot(H)
fd,ad=getFig('d*')
ad.plot(dstr,'b')
ad.plot(delta,'r')
ft,at=getFig('Theta')
at.plot(theta)
a.plot(xbl,ybl,lw=2,color='black',linestyle='--')
#%%
#def myIntegral (u,l):
#    n=len(u)
#    U=np.zeros_like(u)
#    for i in range(1,n):
#        U[i]=0.5*(l[i]-l[i-1])*(u[i]+u[i-1])+U[i-1]
#    return U