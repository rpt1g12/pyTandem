import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *# -*- coding: utf-8 -*-
pi=np.pi
plt.close('all')
#%%
f =lambda P,Q: (1/3)*P*(Q-(2/9)*P**2)-(2/27)*(P**2-3*Q)**1.5
g =lambda P,Q: (1/3)*P*(Q-(2/9)*P**2)+(2/27)*(P**2-3*Q)**1.5
r =lambda P,R: R/P

#%%

P=[0.0207463,-0.183008,-0.116876,-0.0276233,-0.0284057,-0.126413,-0.0101807,-0.037246,0.0146103,0.0131832,
   0.0232029,-0.0102549,-0.0369838,-0.028263,-0.227014]
Q=[-0.023682,-0.364932,-2.87091,0.0465961,0.0752722,-0.317069,-0.658934,-0.268602,-0.116358,-0.332687,
   -0.0655772,-0.659182,0.062802,0.0383862,-0.318144]
R=[-0.00127918,-3.93063,1.05153,0.00692368,0.00972644,-0.979061,0.164935,-0.0382219,0.00557037,0.0580009,
   -0.0030519,0.164999,0.00258864,-0.00018194,-0.165668]



for n in [-1]:
    Q2=np.linspace(-3,P[n]**2/3,512)
    Q3=np.linspace(-3,P[n]**2/3,512)
    
    fig,ax=getFig('Q-R({:2d}); P,Q,R=({:06.2f},{:06.2f},{:06.2f})'.format(n,P[n],Q[n],R[n]))

    s2=f(P[n],Q2)
    s3=g(P[n],Q3)
    
    R4=np.linspace(-1,1)
    s4=r(P[n],R4)
    
    ax.plot(s2,Q2,c='black')
    ax.plot(s3,Q3,c='black')
    ax.plot(R4,s4,c='black',ls='--')
    
    ax.plot(R[n],Q[n],c='red',marker='x',ms=16)
    ax.grid(True)
    fct=1.1
    ax.set_xlim(-fct*abs(R[n]),fct*abs(R[n]))
    ax.set_ylim(-fct*abs(Q[n]),fct*abs(Q[n]))
    ax.axvline(c='black')
    axShow(ax)

