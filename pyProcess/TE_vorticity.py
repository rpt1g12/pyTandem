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
path='/home/'+user+'/Desktop/post/8A00W11AoA20/vsmallDomain/'
#fl0=p3d.flow(path,"grid.xyz","solTavgQ+W+D.qa")
fl0=p3d.flow(path,"grid.xyz","solTA.qa")
vnames=['Q','wx','wy','wz','D']
fl0.rdHdr()
fl0.rdGrid()
fl0.rdSol(vnames=vnames)
botom_0=((fl0.blk[2].var['wz'].avgDir(2))**2+(fl0.blk[2].var['wy'].avgDir(2))**2+(fl0.blk[2].var['wx'].avgDir(2))**2)**0.5
#botom_0=fl0.blk[2].var['D'].avgDir(2)
botom_0_x=fl0.blk[2].var['x'].avgDir(2)
botom_0_y=fl0.blk[2].var['y'].avgDir(2)
top_0=((fl0.blk[5].var['wz'].avgDir(2))**2+(fl0.blk[5].var['wy'].avgDir(2))**2+(fl0.blk[5].var['wx'].avgDir(2))**2)**0.5
#top_0=fl0.blk[5].var['D'].avgDir(2)
top_0_x=fl0.blk[5].var['x'].avgDir(2)
top_0_y=fl0.blk[5].var['y'].avgDir(2)
#%%
#f,a=getFig('wz')
#lvl=np.linspace(-30,30,16)
#a.contourf(botom_0_x,botom_0_y,botom_0,levels=lvl,exted='both')
#a.contourf(top_0_x,top_0_y,top_0,levels=lvl,exted='both')
#axShow(a) 
#%%
path='/home/'+user+'/Desktop/post/8A15W11AoA20/vsmallDomain/'
#fl1=p3d.flow(path,"grid.xyz","solTavgQ+W+D.qa")
fl1=p3d.flow(path,"grid.xyz","solTA.qa")
vnames=['Q','wx','wy','wz','D']
fl1.rdHdr()
fl1.rdGrid()
fl1.rdSol(vnames=vnames)
botom_1=((fl1.blk[2].var['wz'].avgDir(2))**2+(fl1.blk[2].var['wy'].avgDir(2))**2+(fl1.blk[2].var['wx'].avgDir(2))**2)**0.5
#botom_1=fl1.blk[2].var['D'].avgDir(2)
botom_1_x=fl1.blk[2].var['x'].avgDir(2)
botom_1_y=fl1.blk[2].var['y'].avgDir(2)
top_1=((fl1.blk[5].var['wz'].avgDir(2))**2+(fl1.blk[5].var['wy'].avgDir(2))**2+(fl1.blk[5].var['wx'].avgDir(2))**2)**0.5
#top_1=fl1.blk[5].var['D'].avgDir(2)
top_1_x=fl1.blk[5].var['x'].avgDir(2)
top_1_y=fl1.blk[5].var['y'].avgDir(2)
#%%
f1,a1=getFig(r'$\Delta{U}$')
lvl=np.linspace(-0.2,0.2,11)
a1.contourf(botom_1_x,botom_1_y,botom_1-botom_0,levels=lvl,exted='both')
im=a1.contourf(top_1_x,top_1_y,top_1-top_0,levels=lvl,exted='both')
a1.plot(fl0.blk[1].var['x'].val[:,-1,0],fl0.blk[1].var['y'].val[:,-1,0],lw=2,color='black')
a1.plot(fl0.blk[4].var['x'].val[:,0,0],fl0.blk[4].var['y'].val[:,0,0],lw=2,color='black')
a1.set_xlim(0.3,0.8)
a1.set_ylim(-0.1,0.4)


cb=f1.colorbar(im, orientation='vertical')
cb.set_label(r'$U|_{WLE}-U|_{SLE}$',fontsize=20)
a1.set_xlabel(r'$x$',fontsize=20)
a1.set_ylabel(r'$y$',fontsize=20)
axShow(a1)