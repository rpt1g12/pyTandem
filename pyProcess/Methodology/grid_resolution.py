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
plt.close('all');
figs=[];axs=[];hdls=[];lbls=[];lgds=[];nfig=-1;fs=18
import importlib
#%%
importlib.reload(p3d)
#%% Options
save=False
A=0 #WLE Amplitude, if SLE A=0
AoA=10 #Angle of Attack
nwave=8 #Number of WLE wavelenghts
topBlk=1 # Suction side block
botBlk=0 # Pressure side block
if A>0:
    sim='WLE'
    slices=[0,12,24,36]
else:
    sim='SLE'
    slices=np.asarray(range(0,48,1)) # Array of slice indices

nslice=len(slices) # Number of slices
#%% Paths set-up
if sim=='SLE':
    subpath='average/ss003/'
else:
    subpath='average/ss005/'
simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/sonyHDD/post/{}/{}".format(user,simfolder,subpath)
spath='/media/rperezt/082F-63FE/phd/thesis/data/Methodology/WallUnits/'
print('Reading data from:\n {}'.format(path))
#%%
vnames=['r','u','v','w','p']; # Variable names
gfile='grid.xyz' # Grid file name
files=p3d.getFileNames(path=path,pattern='solTA*')
sfile=files[0]#'solTavgCf+tw+Cp.q'# Solution file name

fl=p3d.flow(path,"grid.xyz",sfile) # Check-out flow object
fl.rdHdr() # Read header from grid file and set-up flow object properties
fl.rdGrid() # Read Grid
fl.rdSol(vnames=vnames) # Read solution file
flInfo=fl.rdFlowInfo() # Read solution file info, i.e. Mach, AoA, Re and time
Re = flInfo[2]
#%%
bk = fl.blk[4]
bk.getViscosity(Re = Re)
#%%
x_coord = fl.blk[4].var['x'].getValues()
y_coord = fl.blk[4].var['y'].getValues()
z_coord = fl.blk[4].var['z'].getValues()
nx,ny,nz = bk.size
#%%
y1 = np.sqrt( (y_coord[:,0,:]-y_coord[:,1,:])**2 )
y1 = y1/2
x1 = np.zeros_like(y1)
z1 = np.ones_like(y1)*(z_coord[0,0,0] - z_coord[0,0,1])
for i in range(nx-1):
    x1[i,:] = np.sqrt( (x_coord[i,0,:]-x_coord[i+1,0,:])**2 )
x1[-1,:] = x1[-2,:]
print(x1.max(),y1.max(),z1.max())
flag = 0
#%%
vnames=['cf','twx','twy','twz','cp']; # Variable names
gfile='grid.xyz' # Grid file name
files=p3d.getFileNames(path=path,pattern='solTavgCf+tw*')
sfile=files[0]#'solTavgCf+tw+Cp.q'# Solution file name

fl=p3d.flow(path,"grid.xyz",sfile) # Check-out flow object
fl.rdHdr() # Read header from grid file and set-up flow object properties
fl.rdGrid() # Read Grid
fl.rdSol(vnames=vnames) # Read solution file
bk2 = fl.blk[4]
#%%
plt.close('all')
f,a = getFig("Wall_units_{}".format(simfolder))
name = '{:1d}{}{:d}.dat'.format(nwave,sim,AoA)
x_plot = np.sum(x_coord[:,0,:],axis=1)/nz
for flag in range(3):
    nu = bk.var['nu'].getValues()[:,0,:]
    rho = bk.var['r'].getValues()[:,0,:]
    tw =  (np.sqrt(bk2.var['twx'].getValues()**2+bk2.var['twy'].getValues()**2+bk2.var['twz'].getValues()**2))[:,0,:]
    u_star = np.sqrt(tw/rho)
    if flag == 2:
        label="z+"
        d1 = z1
    elif flag == 1:
        label="y+"
        d1 = y1
    else:
        label="x+"
        d1 = x1
    y_plus = ( u_star * d1 )/nu
    y_plus = np.sum(y_plus,axis=1)/nz
    print(y_plus.max(),y_plus.min(),y_plus.mean())
    a.plot(x_plot,y_plus,lw=2,label=label)

fit(a)
handle,labels,legend = getLabels(a)
axShow(a)
savePlotFile(path=spath,ax=a,varx=['x'],vary=labels,name=name,sameX=True)