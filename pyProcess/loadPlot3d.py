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
fl=p3d.flow("/home/rperezt/Desktop/post/8A00W11AoA10/vsmallDomain/","grid.xyz","solT348.0002.q")

fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=['r','u','v','w','p'])


fl.blk[4].data
#%%
surf=fl.blk[4].getSubset(ylim=[0])
x=surf.var['x'].getValues()

z=surf.var['z'].getValues()

p=surf.var['r'].getValues()

pavg=surf.var['p'].avgDir()[:,0]
xavg=surf.var['x'].avgDir()[:,0]

plt.plot(xavg,pavg)

#plt.contourf(x[:,0,:],z[:,0,:],p[:,0,:])
#plt.axis('equal')
#print(x.shape)

