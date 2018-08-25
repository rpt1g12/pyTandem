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
#path='/media/'+user+'//post/8A15W11AoA20/fullDomain/T1_plane/'
path='/media/'+user+'/My Passport/results/6Blocks/post/8A15W11AoA20/fullDomain/T1_plane/'
fl=p3d.flow(path,"grid.xyz","Q+W+Cp/solTQ+W+Cp0000.q")


vnames=['r','u','v','w','p']

fl.rdHdr()
fl.rdGrid()
y3 = fl.blk[3].var['y'].getValues().copy()*-1
y4 = fl.blk[4].var['y'].getValues().copy()*-1
#%% Create a second flow object
fl2 = p3d.flow(path=path,gfile="grid2.xyx",sfile="test.q")
fl2.setHdr(nbks=[2,2,1],lxib=[121,321],letb=[361,361],lzeb=[1])
#%%
fl2.blk[2]=fl.blk[3].getSubset()
fl2.blk[3]=fl.blk[4].getSubset()
fl2.blk[0]=fl.blk[3].getSubset();fl2.blk[0].setData(vname='y',val=y3)
fl2.blk[1]=fl.blk[4].getSubset();fl2.blk[1].setData(vname='y',val=y4)

#%%
f,a=getFig('mesh')
for bk in range(4):
    fl2.blk[bk].drawMeshPlane(direction=2,pln=0,skp=[3],ax=a,showBlock=True)
#%%
a.set_ylim(-0.1,0.1)
a.set_xlim(-0.6,-0.4)
a.set_aspect('equal')
axShow(a)
savePath = '/media/rperezt/082F-63FE/phd/thesis/figures/Methodology/BoundaryConditions/'
saveFigOnly(path=savePath,fig=f,ax=a,name='leZoom2',ext='.pdf')
#f.savefig('pgfPlots/leZoom.pdf', bbox_inches='tight', pad_inches=0)