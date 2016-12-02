# -*- coding: utf-8 -*-
import time
import numpy as np
import matplotlib.pyplot as plt
from lib.myPlots import *
from lib.stats import *
from scipy.signal import butter, filtfilt
pi=np.pi
#%%
n=101;r=0.5;x0,y0=(-0.0,-0.0)
x=np.linspace(-1,1,n);y=x.copy()
X,Y=np.meshgrid(x,y)
g=X.copy();R2=X.copy()
k0=(np.log(0.0001)/(r**2))
s=-2*(2*k0*r*np.exp(k0*r**2))
T=np.linspace(0,1,51)
#fig0,sur=getFig('Contour')
fig,ax=getFig('Slice')
plt.ion()
for t in T:
    sin=np.sin(2*pi*t)
    for i in range(n):
        for j in range(n):
            R2[i,j]=(X[i,j]-x0)**2+(Y[i,j]-y0)**2
            g[i,j]=np.exp(k0*R2[i,j])*sin
            #g[i,j]=(2*k0*(2*k0*R2[i,j]+1)*np.exp(k0*R2[i,j]))*s*sin
    
    
    #sur.contourf(X,Y,g)

    
    #mx,my=(np.argmax(g,0)[0],np.argmax(g,1)[0])
    
    cll(ax)
    ax.plot(X[50,:],g[50,:],'b')
    ax.figure.canvas.draw()
    ax.set_ylim(-1,1)
    st='{:8.5f}'.format(t)
    ssin='{:8.5f}'.format(sin)
    print(st+' '+ssin)
    plt.pause(0.1)