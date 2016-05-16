from lib.plots2d import *
import matplotlib.pyplot as plt
import numpy as np
import lib.stats as stats
from lib.indx import indx

#%%
p1=plot2d('cfData/Cf_forced50.dat')
p2=plot2d('cfData/Cf_inflow75.dat')
p3=plot2d('cfData/Cf_grid2force50.dat')
p4=plot2d('cfData/Cf_grid3force100.dat')
p5=plot2d('cfData/Cf_grid4force100.dat')
p1.x+=0.5
p2.x+=0.5
p3.x+=0.5
p4.x+=0.5
p5.x+=0.5
span=0.2
h=span/50
h2=span/75
h3=0.22/100.0
p3y=p3.y[:]*0.4

#%%
z=np.array([i*h-0.1 for i in range(0,51)])
z2=np.array([i*h2-0.1 for i in range(0,76)])
z3=np.array([i*h3-0.1 for i in range(0,101)])
r=z[:]*0.0
r2=z2[:]*0.0
r3=z3[:]*0.0
u=np.array([0.0 for i in range(0,151)])
v=np.array([0.0 for i in range(0,151)])
w=np.array([0.0 for i in range(0,201)])
w3=np.array([0.0 for i in range(0,301)])
w2=np.array([0.0 for i in range(0,201)])
w4=np.array([0.0 for i in range(0,321)])

#%%
for i in range(0,151):
    for k in range(0,51):
        l=indx(i,k,150)
        r[k]=p1.y[l]
    u[i]=stats.avg(r,z)
for i in range(0,151):
    for k in range(0,51):
        l=indx(i,k,150)
        r[k]=p2.y[l]
    v[i]=stats.avg(r,z)
for i in range(0,201):
    for k in range(0,51):
        l=indx(i,k,200)
        r[k]=p3.y[l]
    w[i]=stats.avg(r,z)
for i in range(0,201):
    for k in range(0,76):
        l=indx(i,k,200)
        r2[k]=p4.y[l]
    w2[i]=stats.avg(r,z)
for i in range(0,301):
    for k in range(0,76):
        l=indx(i,k,300)
        r2[k]=p5.y[l]
    w3[i]=stats.avg(r,z)
    
#%%
px,py=np.loadtxt('cfData/Jones2008.dat',skiprows=1,unpack=True)

#%%
plt.rc('text', usetex=False)
plt.rc('font', family='serif')
plt.subplot(111)
plt.plot(p1.x[0:151],u,'b-s',label='Forced',linewidth=3.0,markevery=5)
plt.plot(p2.x[0:151],v,'r-d',label='Inflow',linewidth=3.0,markevery=5)
plt.plot(p3.x[0:201],w,'g-^',label='G2',linewidth=3.0,markevery=5)
plt.plot(p4.x[0:201],w2,'r-s',label='G3',linewidth=3.0,markevery=5)
plt.plot(p5.x[0:301],w3,'c-d',label='G4',linewidth=3.0,markevery=5)
plt.plot(px,py,'ko',label='3DU',markersize=5.0)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),
           loc=3,ncol=3,
           mode="expand", borderaxespad=0.)
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.grid(True)
#plt.gca().invert_yaxis()
plt.gca().yaxis.set_ticks(np.arange(-0.02, 0.09, 0.01))
plt.xlabel(r'$x/L_c$',fontsize=16)
plt.ylim([-0.02,0.02])
plt.ylabel(r'$C_f$',fontsize=16)
plt.show()
