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
save=True
path='/home/'+user+'/Desktop/post/naca0012G3/ss003/'

files=p3d.getFileNames(path=path)

vnames=['r','u','v','w','p'];Re=125000;M=0.4


fl=p3d.flow(path,"grid.xyz",'solTA.qa')
fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames)
up=fl.blk[4]
up.getMetrics()
up.derive('u','y')
up.getViscosity()
bot=fl.blk[1]
bot.getMetrics()
bot.derive('u','y')
bot.getViscosity()

dudy_up=up.var['dudy'].avgDir()[:,0]
nu_up=up.var['nu'].avgDir()[:,0]
tw_up=nu_up*dudy_up/Re
cf_up=tw_up*2/(M**2)
x_up=up.var['x'].getValues()[:,0,0]

dudy_bot=bot.var['dudy'].avgDir()[:,-1]
nu_bot=bot.var['nu'].avgDir()[:,-1]
tw_bot=nu_bot*dudy_bot/Re
cf_bot=-tw_bot*2/(M**2)
x_bot=bot.var['x'].getValues()[:,-1,0]

jx,jcf=np.loadtxt('cfData/Jones2008_bot.dat',skiprows=1,unpack=True)

f,a=getFig('Cf_dudy_G3')

a.plot(x_up,cf_up,lw=2,color="blue",label='up')
a.plot(x_bot,cf_bot,lw=2,color="blue",label='bot')
#a.plot(jx-0.5,jcf,lw=2,color="red")
fit(a)

if save:
    savePlotFile(ax=a)