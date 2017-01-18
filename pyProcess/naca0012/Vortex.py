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
pi=np.pi
import getpass
user=getpass.getuser()
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')
#cosa,sina=np.cos(np.deg2rad(5)),np.sin(np.deg2rad(5))
import importlib
#%%
importlib.reload(p3d)
#%%
Gamma=1;R=0.5
nx,ny=201,201
hnx,hny=int((nx-1)/2),int((ny-1)/2)
x,y=np.linspace(-1,1,nx),np.linspace(-1,1,ny)
X,Y=np.meshgrid(x,y)
u,v=np.zeros_like(X),np.zeros_like(Y)
p=np.zeros_like(u)
U=np.zeros_like(u)
k=np.log(1e-6)/(R**2)
x0=np.sqrt(-1/(2*k))
f0=x0*np.exp(-0.5)
Gamma=3e-8/f0
for j in range(ny):
    for i in range(nx):
        r=x[i]*x[i]+y[j]*y[j]
        A=Gamma*np.exp(k*r)
        u[j,i],v[j,i]=-A*Y[j,i],A*X[j,i]
        U[j,i]=(u[j,i]**2+v[j,i]**2)**0.5
        p[j,i]=-(U[j,i]**2)*0.5
        
#%%
L=x0
var=U
vmin,vmax=var.min(),var.max()
print(vmin,vmax,k,x0,f0)
f,a=getFig()
a.set_aspect('equal')
a.contourf(X,Y,var,extend='both',vmin=vmin,vmax=vmax,levels=np.linspace(vmin,vmax,256))
a.quiver(X,Y,u,v)
a.streamplot(X,Y,u,v)
a.set_xlim(-L,L)
a.set_ylim(-L,L)
#%%
f2,a2=getFig()
a2.plot(X[hny,:],v[hny,:])
a2.plot(X[hny,:],p[hny,:])