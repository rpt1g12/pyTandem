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
save=False
path='/home/'+user+'/Desktop/post/naca0012G1/ss003/'

files=p3d.getFileNames(path=path)

vnames=['r','u','v','w','p'];Re=125000;M=0.4


fl=p3d.flow(path,"grid.xyz",'solTA.qa')
fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames)
up=fl.blk[1]
up.getMetrics()
up.derive('u','y')
up.getViscosity()

dudy=up.var['dudy'].avgDir()[:,-1]
nu=up.var['nu'].avgDir()[:,-1]
tw=nu*dudy/Re
cf=tw*2/(M**2)

x=up.var['x'].getValues()[:,-1,0]

jx,jcf=np.loadtxt('cfData/Jones2008_up.dat',skiprows=1,unpack=True)

f,a=getFig('Cf_dudy')

a.plot(x,cf,lw=2,color="blue")
a.plot(jx-0.5,jcf,lw=2,color="red")
fit(a)

if save:
    savePlotFile(ax=a,vary=['cfdudy'])