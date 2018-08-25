# -*- coding: utf-8 -*-

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

import importlib
#%%
importlib.reload(p3d)
#%%
path ="/media/rperezt/My Passport/results/6Blocks/post/A00W11AoA20/fullDomain/"
path ="/media/rperezt/My Passport/results/6Blocks/post/A00W11AoA20/fullDomain/"
fl = p3d.flow(path,"grid.xyz","solTQ+W000.qa")

vnames=["sp0","sp1"]
fl.rdHdr()
fl.rdGrid()
fl.rdFun(vnames=vnames,ffile="sponge000.f")
#%%
fl.contourf(varname="sp0",vmin=0,vmax=12,k=1)