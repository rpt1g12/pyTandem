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

path='/home/'+user+'/fortran/tandemProject/out/ss001/'

files=p3d.getFileNames(path=path)

nfiles=len(files)

#nfiles=1

vnames=['r','u','v','w','p']

ix=85;jy=25

f,a=getFig('Velocity Profile i={:d}'.format(ix))

p_inf=1/1.4
u_inf=0.4
Re=125000

fl=p3d.flow(path,"grid.xyz",files[0])
fl.rdHdr()
fl.rdGrid()
up=fl.blk[4]
y1=up.var['y'].getValues()[:,1,0]-up.var['y'].getValues()[:,0,0]
up.getMetrics()
for n in range(nfiles):    

    fl.rdSol(vnames=vnames,sfile=files[n])
    

    
    #for ix in range(0,90,5):
#    y=up.var['y'].getValues()[ix,:jy,0]
#    u=up.var['u'].getValues()[ix,:jy,0]
#    v=up.var['v'].getValues()[ix,:jy,0]
#    w=up.var['w'].getValues()[ix,:jy,0]
#    U=np.sqrt(u**2+v**2+w**2)
#    a.plot(u,y-y[0],lw=2,marker='o')
    if n==0:
        x=up.var['x'].getValues()[:,0,0]
        p=up.var['p'].getValues()[:,0,0]
        dudy=up.derive('u','y',True)[:,0,0]
    else:
        p+=up.var['p'].getValues()[:,0,0]
        dudy+=up.derive('u','y',True)[:,0,0]

p_avg=p/nfiles

cp=(p_avg-p_inf)/(0.5*u_inf**2)

dudy_avg=dudy/nfiles

cf=(dudy_avg)/(Re*0.5*u_inf**2)

ustar=np.sqrt(2*u_inf**2*np.abs(cf))

yplus=ustar*y1*Re

a.plot(x,cf,lw=2,marker='o')