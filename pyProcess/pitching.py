import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
pi=np.pi
plt.close('all')
#%%
mach=0.3
aoa0=10*pi/180;aoa1=80*pi/180
u0=np.array([mach*np.cos(aoa0),mach*np.sin(aoa0)])
up=np.array([mach*np.cos(aoa1),mach*np.sin(aoa1)])


t0=0;t1=11;tps=2;tp=4;tpe=tps+tp
nt=(t1-t0-1)*10+1
t=np.linspace(t0,t1-1,nt)
ut=np.array([u0 for i in range(nt)])
taoa=np.array([aoa0*180/pi for i in range(nt)])
taoaf=taoa.copy()
tmach=np.array([mach for i in range(nt)])
tmachf=tmach.copy()
g=np.zeros(nt)
dts=pi/tp

k=(1-np.cos(aoa1-aoa0))
for i in range(nt):
    ts=max(min(pi/tp*(t[i]-tps),pi),0);
    ut[i]=u0+(up-u0)*0.5*(1-np.cos(ts))
    taoa[i]=180*np.arctan(ut[i,1]/ut[i,0])/pi
    taoaf[i]=180*(aoa0+(aoa1-aoa0)*0.5*(1-np.cos(ts)))/pi
    tmach[i]=np.linalg.norm(ut[i])
    g[i]=((np.cos(ts))**2-1)/4
    tmachf[i]=((1+2*k*g[i])*mach**2)**0.5

fig,ax=getFig('Pitching')
ax2=ax.twinx()
ax.plot(t[:],ut[:,0],'b')
ax.plot(t[:],ut[:,1],'r')
ax2.plot(t[:],taoa[:],'g')
ax2.plot(t[:],taoaf[:],'og')
ax.plot(t[:],tmach[:],'m')
ax.plot(t[:],tmachf[:],'om')
ax.plot(t[:],g[:],'--k')
