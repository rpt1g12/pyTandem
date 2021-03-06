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
path='/media/'+user+'/My Passport/results/6Blocks/post/8A15W11AoA20/vsmallDomain/'
gfile='grid.xyz'
sfile='solTA.qa'
vnames=['uu','vv','ww','uv','uw','vw']

fl=p3d.flow(path,gfile,sfile)
fl.rdHdr()
fl.rdGrid()