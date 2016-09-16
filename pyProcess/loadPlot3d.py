# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
import lib.Plot3DClasses as p3d
from scipy import stats
pi=np.pi
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')

import importlib
importlib.reload(p3d)

#%%
#fl=p3d.flow("/home/rpt1g12/Desktop/post/A4A15W11AoA20/f64/","grid.xyz","solTA.qa")
fl=p3d.flow("/home/rperezt/Desktop/post/A4A15W11AoA20/vsmallDomain/f64/","grid.xyz","solTA.qa")

fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=['r','u','v','w','p'])


fl.blk[4].data
#%%
xlim=range(0,321);ylim=range(1);zlim=range(1)

surf=fl.blk[4].getSubset(ylim=ylim)
x=surf.var['x'].getValues()
print(x[0,0,0])
z=surf.var['z'].getValues()
print(z[0,0,0])
p=surf.var['p'].getValues()
print(p[0,0,0])


plt.contourf(x[:,0,:],z[:,0,:],p[:,0,:])
plt.axis('equal')
#print(x.getValues().shape)

