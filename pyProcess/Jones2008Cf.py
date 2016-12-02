from lib.plots2d import *
import matplotlib.pyplot as plt
import numpy as np
import lib.stats as stats
from lib.indx import indx
from lib.myPlots import *
plt.close('all')
cosa=np.cos(np.deg2rad(5))
sina=np.sin(np.deg2rad(5))

#%%
p1=plot2d('cfData/Cf_forced50.dat')
p2=plot2d('cfData/Cf_inflow75.dat')
p3=plot2d('cfData/Cf_grid2force50.dat')
p4=plot2d('cfData/Cf_grid3force100.dat')
p5=plot2d('cfData/Cf_grid4force100.dat')
x1=p1.x;y1=-np.tan(np.deg2rad(5))*(x1+0.5)
x=x1*cosa-y1*sina;x1=x.copy()
x2=p2.x;y2=-np.tan(np.deg2rad(5))*(x2+0.5)
x=x2*cosa-y2*sina;x2=x.copy()
x3=p3.x;y3=-np.tan(np.deg2rad(5))*(x3+0.5)
x=x3*cosa-y3*sina;x3=x.copy()
x4=p4.x;y4=-np.tan(np.deg2rad(5))*(x4+0.5)
x=x4*cosa-y4*sina;x4=x.copy()
x5=p5.x;y5=-np.tan(np.deg2rad(5))*(x5+0.5)
x=x5*cosa-y5*sina;x5=x.copy()
span=0.2
h=span/50
h2=span/75
h3=0.22/100.0

#%%
z=np.array([i*h-0.1 for i in range(0,51)])
z2=np.array([i*h2-0.1 for i in range(0,76)])
r=z.copy()
r2=z2.copy()
u=np.zeros(151)
v=np.zeros(151)
w=np.zeros(201)
w3=np.zeros(301)
w2=np.zeros(201)

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
    w2[i]=stats.avg(r2,z2)
for i in range(0,301):
    for k in range(0,76):
        l=indx(i,k,300)
        r2[k]=p5.y[l]
    w3[i]=stats.avg(r2,z2)
    
#%%
px,py=np.loadtxt('cfData/Jones2008.dat',skiprows=1,unpack=True)
px+=-0.5
#%%
f,a=getFig('cfComparison')
a.plot(x1[0:151],u,label='Forced',lw=3)
a.plot(x2[0:151],v,label='Inflow',lw=3)
a.plot(x3[0:201],w,label='G2',lw=3)
a.plot(x4[0:201],w2,label='G3',lw=3)
a.plot(x5[0:301],w3,label='G4',lw=3)
a.plot(px,py,label='3DU',lw=3)
handle,labels,legend=getLabels(ax=a)
axShow(a)
#%%
save=True
if save:
    savePlotFile(ax=a)