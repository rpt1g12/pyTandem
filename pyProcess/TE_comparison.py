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
fl0=p3d.flow(path,"grid.xyz","solTA.qa")
vnames=['r','u','v','w','p']
fl0.rdHdr()
fl0.rdGrid()
fl0.rdSol(vnames=vnames)

path='/home/'+user+'/Desktop/post/8A15W11AoA20/vsmallDomain/'
fl1=p3d.flow(path,"grid.xyz","solTA.qa")
fl1.rdHdr()
fl1.rdGrid()
fl1.rdSol(vnames=vnames)



path='/home/'+user+'/Desktop/post/8Difference/vsmallDomain/'
dfl=p3d.flow(path,'grid.xyz','solTDiff.qa')



#%%
lxis=[321,81]
lxiss=[321,81,321,81]
lets=[21,121]
letss=[21,21,121,121]


for bk in [1,2,4,5]:
    u0=(fl0.blk[bk].var['u'].getValues())**2+(fl0.blk[bk].var['v'].getValues())**2+(fl0.blk[bk].var['w'].getValues())**2
    u1=(fl1.blk[bk].var['u'].getValues())**2+(fl1.blk[bk].var['v'].getValues())**2+(fl1.blk[bk].var['w'].getValues())**2
    fl0.blk[bk].setData(vname='U',val=u0)
    fl1.blk[bk].setData(vname='U',val=u1)
#%%
dfl.setHdr([2,2,1],lxis,lets,[1])
#%%
varname='p'
for varname in ['p','U']:
    count=0
    for bk in [1,2,4,5]:    
        xi=fl0.blk[bk].var['x'].getValues()
        yi=fl0.blk[bk].var['y'].getValues()
        zi=fl0.blk[bk].var['z'].getValues()
        
        iblk1=fl1.blk[bk].interpolate(varname,xi,yi,zi)
    
        x=np.zeros((lxiss[count],letss[count],1))
        y=x.copy()
        z=x.copy()
        p=x.copy()
        p0=x.copy()
        x[:,:,0]=iblk1.var['x'].avgDir(2)
        y[:,:,0]=iblk1.var['y'].avgDir(2)
        z[:,:,0]=iblk1.var['z'].avgDir(2)
        p[:,:,0]=iblk1.var[varname].avgDir(2)
        p0[:,:,0]=fl0.blk[bk].var[varname].avgDir(2)*np.sqrt(1-0.3**2)
        
        dfl.blk[count].setData(vname='x',val=x)
        dfl.blk[count].setData(vname='y',val=y)
        dfl.blk[count].setData(vname='z',val=z)
        dfl.blk[count].setData(vname=varname,val=p-p0)
        count+=1

#%%
plt.close('all')
varname='U';bval=0.1
f,ax=getFig('Difference_'+varname)
ax.plot(fl0.blk[1].var['x'].val[:,-1,0],fl0.blk[1].var['y'].val[:,-1,0],lw=2,color='black')
ax.plot(fl0.blk[4].var['x'].val[:,0,0],fl0.blk[4].var['y'].val[:,0,0],lw=2,color='black')
bars=[False,False,False,False]
for count in range(4):
    dfl.blk[count].contourf(varname=varname,vmin=-bval,vmax=bval,ax=ax,bar=bars[count],nlvl=21)

ymin=-0.05
ax.set_xlim(0.25,0.75)
ax.set_ylim(ymin,ymin+0.5)
axShow(ax)
ax.axes.get_xaxis().set_visible(False)
ax.axes.get_yaxis().set_visible(False)


