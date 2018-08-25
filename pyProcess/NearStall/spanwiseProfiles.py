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
side='both'
var='p' # Variable to plot
A=15 #WLE Amplitude, if SLE A=0
AoA=17 #Angle of Attack
nwave=8 #Number of WLE wavelenghts
topBlk=4 # Suction side block
botBlk=1 # Pressure side block
xpos = -0.45
#%% Paths set-up
media = "082F-63FE"
simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/{}/post/{}/ss001/".format(user,media,simfolder)
spath='/media/{}/{}/phd/thesis/figures/nearStall/{}spanForced/'.format(user,media,var)
if not os.path.exists(spath) and save:
    os.makedirs(spath)
print('Reading data from:\n {}'.format(path))
#%%
vnames=['r','u','v','w','p']; # Variable names
gfile='grid.xyz' # Grid file name
files=p3d.getFileNames(path=path,pattern='solT*')
sfile=files[0]#'solTavgCf+tw+Cp.q'# Solution file name
nt = len(files)
fl=p3d.flow(path,"grid.xyz",sfile) # Check-out flow object
fl.rdHdr() # Read header from grid file and set-up flow object properties
fl.rdGrid() # Read Grid
fl.rdSol(vnames=vnames) # Read solution file
flInfo=fl.rdFlowInfo() # Read solution file info, i.e. Mach, AoA, Re and time
mach = flInfo[0]
#%% Extract slices
t=[]
iv_top=np.zeros((nt,nze))
nze=fl.blk[botBlk].size[2]
nxi=fl.blk[botBlk].size[0]
zi = fl.blk[topBlk].var['z'].val[0,0,:]
for i in range(nt):
    fl.rdSol(vnames=vnames,sfile=files[i])
    flInfo=fl.rdFlowInfo()
    t.append(flInfo[-1])
    cp = fl.blk[topBlk].var['p'].getValues()[:,0,:]
    cp =  (cp - (1/1.4))/(0.5*mach**2)
    xi = np.ones(nze) * xpos
    size=nxi*nze
    x=fl.blk[topBlk].var['x'].val[:,0,:];x=np.reshape(x,x.size)
    z=fl.blk[topBlk].var['z'].val[:,0,:];z=np.reshape(z,z.size)
    points=np.array(np.transpose([x,z]))
    values=cp;values=np.reshape(values,size)
    ipoints=np.array(np.transpose([np.reshape(xi,xi.size),np.reshape(zi,zi.size)]))
    iv_top[i,:]=griddata(points,values,ipoints,method="cubic")

t = np.asarray(t)
t = t - t[0]
#%%
fname='{}{}'.format(simfolder,var)
f,a=getFig(fname);nfig+=1 # Get figure and axis objects
figs.append(f);axs.append(a) # Append them to the figures and axes arrays
data = np.zeros_like(iv_top)
for i in range(nt):
    corr = acorrp(iv_top[i,:])
    data[i,:] = corr
    a.plot(zi,corr,lw=2,color="blue")
    a.set_ylim(-1,1)
    plt.pause(0.05)
    cll(a)
#%%
f,a = getFig("pCorrForced")
a.contourf(t,np.linspace(0,8,nze),data.transpose(),255,cmap=plt.cm.bwr,vmin=-1,vmax=1) 
#%%
saveFigOnly(path=spath)