# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
from lib.myPlot3dOperator import *
import lib.Plot3DClasses as p3d
import lib.Plot3DShapes as shapes
from scipy import stats
from scipy.interpolate import griddata
from scipy.signal import argrelextrema
import getpass
user='/home/{:}/'.format(getpass.getuser())
pi=np.pi
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')

import importlib
#%%
importlib.reload(p3d)
importlib.reload(shapes)
#%%
path=user+'Desktop/post/A4A15W11AoA20/f64/'
gfile='grid.xyz';sfile='solTA.qa'

fl=p3d.flow(path,gfile,sfile)
fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=['r','u','v','w','p'])

fl.blk[4].wrVarASCII(varnames=['x','y','u','v'],k=87,fpath=path)

#%%
path=user+'Desktop/post/8A15W11AoA20/aerofoil/'
gfile='grid.xyz';sfile='solT180.0002.q'

fl=p3d.flow(path,gfile,sfile)
fl.rdHdr()
fl.rdGrid()

mfl=fl.mergeBlocks(bkx=1,bky=2)
mfl.blk[0].wrVarASCII3D(varnames=['x','y','z'],fpath=path)