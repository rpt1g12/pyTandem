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
path=user+'Desktop/post/8A15W11AoA20/fullDomain/'
gfile='grid.xyz';sfile='solTRMS.qa'
#c=shapes.circle(user+'Desktop/circle/','circle.xyz','solCircle.q',0,0,0,1,0,270,45)

L_wAvg=0;rngK=range(0,200,10)

f,a=getFig('Wake_Contour')
f1,a1=getFig('Wake_Line')

fl=p3d.flow(path,gfile,sfile)
fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=['r','u','v','w','p'])
#%%
for k in rngK:
    print('k={:}'.format(k))
    pl19=fl.blk[5].getSubset(zlim=[k])
    
    
    #%%
#    if k==0:
#        pl19.contourf(varname='p',vmin=-0.01,vmax=0.01,nlvl=256,ax=a)
#    else:
#        pl19.contourf(varname='p',vmin=-0.01,vmax=0.01,nlvl=256,ax=a,bar=False)
    
    #%%
    z0=pl19.var['z'].getValues()[0,0,0]
#    pTE=np.asarray([0.5,0,z0])
#    dr=np.asarray([np.cos(np.deg2rad(20)),np.sin(np.deg2rad(20)),0])
#    dn=np.asarray([-np.sin(np.deg2rad(20)),np.cos(np.deg2rad(20)),0])
#    p0=pTE+(1.25*dr)-(0.5*dn)
#    
#    p1=p0+(1.5*dn)
    
    p0=[.75,0,z0]
    p1=[.75,1,z0]
    
    l=shapes.line(user+'Desktop/circle/','line.xyz','solLine.q',list(p0),list(p1),101)
    l.blk[0].drawMeshPlane(ax=a)
    
    #%%

    varname='p'
    x=l.blk[0].var['x'].getValues()
    y=l.blk[0].var['y'].getValues()
    z=l.blk[0].var['z'].getValues()
    bk_dummy=pl19.interpolate(varname,x,y,z)
    
    l_len=l.blk[0].var['l'].getValues(True)
    l_val=bk_dummy.var[varname].getValues(True)
    
    
    a1.plot(l_len,l_val)
    l_val_s=smooth(l_val,4)
    a1.plot(l_len,l_val_s)
    
    l_max=argrelextrema(l_val_s,np.greater)
    a1.scatter(l_len[l_max],l_val_s[l_max],color='red',marker='o')
    
    axShow(a1)
    
    #%%
    L_w=l_len[l_max][-1]-l_len[l_max][0]
    print(L_w)
    L_wAvg+=L_w

nk=len(rngK)
L_wAvg/nk
del fl