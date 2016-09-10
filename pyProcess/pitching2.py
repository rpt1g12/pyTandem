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

daoa=(aoa1-aoa0)

for i in range(nt):
    ts=max(min(pi/tp*(t[i]-tps),pi),0);
    taoa[i]=(aoa0+(aoa1-aoa0)*0.5*(1-np.cos(ts)))
    ut[i]=mach*np.array([np.cos(taoa[i]),np.sin(taoa[i])])
    taoaf[i]=180*np.arctan(ut[i,1]/ut[i,0])/pi
    tmach[i]=np.linalg.norm(ut[i])
    tmachf[i]=mach

fig,ax=getFig('Pitching')
ax2=ax.twinx()
ax.plot(t[:],ut[:,0],'b')
ax.plot(t[:],ut[:,1],'r')
ax2.plot(t[:],180*taoa[:]/pi,'g')
ax2.plot(t[:],taoaf[:],'og')
ax.plot(t[:],tmach[:],'m')
ax.plot(t[:],tmachf[:],'om')
