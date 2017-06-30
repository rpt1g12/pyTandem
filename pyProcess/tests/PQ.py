# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import welch as psdw
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
from lib.myPlot3dOperator import *
import lib.Plot3DClasses as p3d
from scipy import stats
from scipy.interpolate import griddata
from scipy.interpolate import interp2d
from scipy.optimize import fmin_bfgs as fmin
pi=np.pi
import getpass
user=getpass.getuser()
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')

import importlib
#%%
importlib.reload(p3d)
#%%
A=15 #WLE Amplitude, if SLE A=0
AoA=20 #Angle of Attack
nwave=4 #Number of LE wavelengths
gfile='grid.xyz' #Grid filename
sfile='solTA.qa' #Solution filename
svnames=['r','u','v','w','p']

#%%
varname='P'; vmin,vmax=0,1

#%% Paths set-up
if A>0:
    sfolder='{}WLE'.format(nwave)
else:
    sfolder='{}SLE'.format(nwave)
subpath='average/'
simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/dellHDD/post/{}/{}".format(user,simfolder,subpath)

fl=p3d.flow(path,gfile,sfile)
fl.rdHdr()
fl.rdGrid()
fl.rdSol(svnames,sfile)
#%%
fl.blk[4].getStrain(invariants=True)
#P=fl.blk[4].var['dudx'].getValues()
#P+=fl.blk[4].var['dvdy'].getValues()
#P+=fl.blk[4].var['dwdz'].getValues()
#fl.blk[4].setData(vname='P',val=P)
#%%
f,a,i=fl.blk[4].contourf(varname='R',vmin=vmin,vmax=1,plane=1,k=0,nlvl=21,cmap=plt.cm.jet,bar=False)