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

#%% Sumulation name and set-up
#hName='8SLE';kk=0;block=0
#path="/media/rpt1g12/dell\'s_hdd/post/8A00W11AoA20/subsets/ss004/"
#tpath='/home/rpt1g12/Documents/thesis/data/StalledWLE/fig24_8SLE/'
#hName='8SLE_150';kk=0;block=0
#path="/media/rpt1g12/dell\'s_hdd/post/8A00W11AoA20/subsets/ss006/"
#tpath='/home/rpt1g12/Documents/thesis/data/StalledWLE/fig24_8SLE/'
#hName='8WLE_Central';kk=0;block=0
#path="/media/rpt1g12/dell\'s_hdd/post/8A15W11AoA20/subsets/ss003/ss001/"
#tpath='/home/rpt1g12/Documents/thesis/data/StalledWLE/fig25_8WLE/'
hName='8WLE_SSL';kk=0;block=0
path="/media/rpt1g12/dell\'s_hdd/post/8A15W11AoA20/subsets/ss007/ss001/"
tpath='/home/rpt1g12/Documents/thesis/data/StalledWLE/fig25_8WLE/'
#%% Options
save=False
#%%

files=p3d.getFileNames(path=path)

vnames=['r','u','v','w','p'];

fl=p3d.flow(path,"grid.xyz",files[0])
fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames)
flInfo=fl.rdFlowInfo()
#%%
flavg=fl.getAvg(files,vnames)
flavg.wrSol(sfile='solTA.qa',vnames=vnames)
#%%
f,a,im=flavg.contourf(varname='uu',vmin=0,vmax=0.02,k=0,nlvl=256,cmap=plt.cm.hot_r);nfig+=1
figs.append(f);axs.append(a)
axs[nfig].set_aspect('equal')
#%%

up=flavg.blk[block]
#%%
var='uu';xprob,yprob=[],[]
xiprob,etprob=[],[]
for i in range(10,180,10):
    print('i={:}'.format(i))
    x,y=up.var['x'].getValues()[i,:,kk],up.var['y'].getValues()[i,:,kk]
    val=up.var[var].getValues()[i,:,kk]
    ip=np.argmax(val)
    xiprob.append(i);etprob.append(ip)
    xprob.append(x[ip]);yprob.append(y[ip])
    axs[nfig].plot(xprob[-1],yprob[-1],marker='x',ms=10,color='blue')
    axs[nfig].plot(x[:],y[:],color='green',lw=2)
    
xprob=np.asarray(xprob)
yprob=np.asarray(yprob)

axShow(axs[nfig])
#%%
ivar=['p','u','v']
nvar=len(ivar)
nprob=len(xprob)
nt=len(files)
#nt=64*60+1
t=np.zeros(nt)
data=np.zeros((nt,nprob,nvar))
for n in range(nt):
    t[n]=fl.rdFlowInfo(sfile=files[n])[3]
    fl.rdSol(sfile=files[n])
    for i in range(nvar):
        var=ivar[i]
        data[n,:,i]=fl.blk[block].var[var].getValues()[xiprob,etprob,kk]
#%% Plot and save histories
for i in range(nvar):
    var=ivar[i]
    f,a=getFig(hName+'_Hist_{}'.format(var));nfig+=1
    figs.append(f);axs.append(a)
    for j in range(nprob):
        axs[nfig].plot(t[:]-t[0],data[:,j,i])
    fit(axs[nfig])
    savePlotFile(ax=axs[nfig],sameX=True)
#%% Plot and Save probe locations

f,a=getFig(hName+'_probes');nfig+=1
figs.append(f);axs.append(a)
axs[nfig].plot(xprob[:],yprob[:],lw=2,color='blue')
axs[nfig].set_aspect('equal')
savePlotFile(ax=axs[nfig],path=tpath)