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
A=15 #WLE Amplitude, if SLE A=0
AoA=20 #Angle of Attack
kk=125 #Number of Spanwise planes to shift
nwave=8 #Number of LE wavelengths
gfile='grid.xyz' #Grid filename
ffile='Rij000.f' #Function filename
sfile='solTCf+tw+Cp.qa' #Solution filename
svnames=['r','u','v','w','p']
fvnames=['uu','vv','ww','uv','uw','vw']
varname=svnames[4]; vmin,vmax=-1,1
sfileExt=sfile.split('.')[-1]
sfile_shift=sfile.split('.'+sfileExt)[0]+'_shifted.'+sfileExt
ffile_shift=ffile.split('.f')[0]+'_shifted.f'
gfile_shift=gfile.split('.xyz')[0]+'_shifted.xyz'
#%% Paths set-up
if A>0:
    sfolder='{}WLE'.format(nwave)
else:
    sfolder='{}SLE'.format(nwave)
subpath='average/vsmallDomain/'
simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/dellHDD/post/{}/{}".format(user,simfolder,subpath)

fl=p3d.flow(path,gfile,sfile)
fl.rdHdr()
fl.rdGrid()
#%%
fl.rdSol(svnames,sfile)
#fl.rdFun(fvnames,ffile)

fl.contourf(varname=varname,vmin=vmin,vmax=vmax,k=kk,nlvl=256,cmap=plt.cm.jet,bar=False)

#sfl=fl.shiftK(kk,['x','y'])
sfl=fl.shiftK(kk,svnames)
#sfl=fl.shiftK(kk,fvnames)


sfl.contourf(varname=varname,vmin=vmin,vmax=vmax,k=0,nlvl=256,cmap=plt.cm.jet,bar=False)

#sfl.wrGrid(gfile_shift)
sfl.wrSol(svnames,sfile_shift)
#sfl.wrFun(fvnames,ffile_shift)