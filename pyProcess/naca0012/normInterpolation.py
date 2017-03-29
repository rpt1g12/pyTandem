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
kk=19
save=False
#path='/media/rperezt/082F-63FE/post/naca0012G1/ss003/'
#path='/media/rperezt/082F-63FE/post/8A15W11AoA20/vsmallDomain/'

files=p3d.getFileNames(path=path)

vnames=['r','u','v','w','p'];

fl=p3d.flow(path,"grid.xyz",'solTA.qa')
fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames)
flInfo=fl.rdFlowInfo()
f,a,im=fl.contourf(varname='p',vmin=0.6,vmax=0.8,k=kk,nlvl=256,cmap=plt.cm.hot)
a.set_aspect('equal')

#%%
#fl.getMetrics()
#up=fl.blk[4]
#etax,etay,etaz=up.mets[3].getValues(),up.mets[4].getValues(),up.mets[5].getValues()
#
#X,Y=up.var['x'].getValues(),up.var['y'].getValues()
#p=up.var['p'].getValues()[:,0,kk]
#p0=(1/1.4)+0.5*0.4**2
#u0=np.sqrt(2*(p0-p))
#nx,ny,xo,yo=etax[:,0,kk],etay[:,0,kk],X[:,0,kk],Y[:,0,kk]
##a.quiver(xo,yo,nx,ny)
#           
##%%
#flJoin=fl.mergeBlocks(3,2)
#U=flJoin.blk[0].var['u'].getValues()**2+flJoin.blk[0].var['v'].getValues()**2+flJoin.blk[0].var['w'].getValues()**2
#U=np.sqrt(U)
#flJoin.blk[0].setData('U',val=U)
#pt=flJoin.blk[0].var['p'].getValues()+0.5*flJoin.blk[0].var['r'].getValues()*U**2
#flJoin.blk[0].setData('pt',val=pt)
#nxi=xo.shape[0]
#
#rakes=[]
#nn=151;l=0.05;dl=l/(nn-1)
#for i in range(50,140,15):
#    rakes.append(p3d.rake(xo[i],yo[i],nx[i],ny[i],nn,l))
#    x,y=rakes[-1].x,rakes[-1].y  
#    a.scatter(x,y,lw=2,color='blue')
#    rakes[-1].var=flJoin.blk[0].interpolate2dk('u',x,y,kk,mode='values')[:,0,0]/u0[i]
#
##%%    
#f1,a1=getFig()
#nr=len(rakes)
#rl=[i*dl for i in range(nn)]
#a1.axvline(x=1,color='black',linestyle='--',lw=2)
#for i in range(nr):
#    a1.plot(rakes[i].var,rl)
#
