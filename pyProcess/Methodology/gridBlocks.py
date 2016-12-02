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
path='/home/'+user+'/Desktop/post/8A15W11AoA20/fullDomain/T1_plane/'
fl=p3d.flow(path,"grid.xyz","Q+W+Cp/solTQ+W+Cp0000.q")


vnames=['r','u','v','w','p']

fl.rdHdr()
fl.rdGrid()

#%%
f,a=getFig('mesh')
for bk in range(6):
    fl.blk[bk].drawMeshPlane(direction=2,pln=0,skp=3,ax=a,showBlock=True)
#%%
a.set_ylim(-0.1,0.1)
a.set_xlim(-0.6,-0.4)
a.set_aspect('equal')
axShow(a)
f.savefig('pgfPlots/leZoom.pdf', bbox_inches='tight', pad_inches=0)