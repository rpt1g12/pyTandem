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
fl=p3d.flow("/home/rpt1g12/Desktop/post/A4A15W11AoA20/f64/","grid.xyz","solTA.qa")

fl.rdHdr()
fl.rdGrid()
fl.rdSol()

#%%
xlim=range(0,321);ylim=range(0,0);zlim=range(0,0)

surf=fl.blk[4].getSubset(xlim,ylim,zlim)
x=surf.var['x']
print(x.getValues())
