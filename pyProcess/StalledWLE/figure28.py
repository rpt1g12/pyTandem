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
#%%Options
avg=False
#%% Sumulation name and set-up
hName='8SLE_T2';kk=0;block=0;nrange=[185,205,230,260]
path="/media/rpt1g12/dell\'s_hdd/post/8A00W11AoA20/subsets/ss004/"
tpath='/home/rpt1g12/Documents/thesis/figures/StalledWLE/fig26_8SLE/'
hName='8SLE_150';kk=0;block=0;nrange=[185,205,230,260]
path="/media/rpt1g12/dell\'s_hdd/post/8A00W11AoA20/subsets/ss006/"
tpath='/home/rpt1g12/Documents/thesis/figures/StalledWLE/fig26_8SLE/'
#hName='8WLE_T4';kk=0;block=0;nn=125
#path="/media/rpt1g12/dell\'s_hdd/post/8A15W11AoA20/subsets/ss003/ss001/"
#tpath='/home/rpt1g12/Documents/thesis/figures/StalledWLE/fig26_8WLE/'
#hName='8WLE_T8';kk=0;block=0;nrange=[70,85,100,125]
#path="/media/rpt1g12/dell\'s_hdd/post/8A15W11AoA20/subsets/ss007/ss001/"
#tpath='/home/rpt1g12/Documents/thesis/figures/StalledWLE/fig26_8WLE/'
#%%
files=p3d.getFileNames(path=path)

vnames=['r','u','v','w','p'];

fl=p3d.flow(path,"grid.xyz",files[0])
fl.rdHdr()
fl.rdGrid()
for nn in nrange:
    fl.rdSol(vnames=vnames,sfile=files[nn])
    flInfo=fl.rdFlowInfo()
    #%%
    var='p'
    for block in range(len(fl.blk)):
        A=2*(fl.blk[block].var[var].getValues()-(1/1.4))/0.3**2
        fl.blk[block].setData(vname='Cp',val=A)
    #%%
    f,a,im=fl.contourf(varname='Cp',vmin=-1.75,vmax=0,k=0,nlvl=1024,cmap=plt.cm.hot,bar=False);nfig+=1
    figs.append(f);axs.append(a)
    a.set_xlim(-0.6,-0.1)
    a.set_ylim(0,0.25)
    axs[nfig].set_aspect('equal')
    #%%
    saveFigOnly(path=tpath,fig=f,ax=a,name='{}_{:06d}'.format(hName,nn))