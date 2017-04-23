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
figs=[];axs=[];hdls=[];lbls=[];lgds=[];nfig=-1;fs=18
import importlib
#%%
importlib.reload(p3d)

#%%
nxi,neta=180,100
save=False
path="/media/rpt1g12/dell's_hdd/post/8A15W11AoA20/subsets/ss003/"

files=p3d.getFileNames(path=path)
nt=len(files)
vnames=['r','u','v','w','p'];

fl=p3d.flow(path,"grid.xyz",files[0])
fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames,sfile=files[0]);flInfo=fl.rdFlowInfo()

sspath="/media/rpt1g12/dell's_hdd/post/8A15W11AoA20/subsets/ss003/ss001/"
ss=p3d.flow(sspath,"grid.xyz",files[0])
#%%
ss.setHdr([1,1,1],[nxi],[neta],[1])
sub=fl.blk[4].getSubset(xlim=range(nxi),ylim=range(neta),link=True)
ss.blk[0]=sub
#%%
ss.wrGrid()
for n in range(nt):
    fl.rdSol(vnames=vnames,sfile=files[n]);flInfo=fl.rdFlowInfo()
    sub=fl.blk[4].getSubset(xlim=range(nxi),ylim=range(neta),link=True)
    ss.blk[0]=sub
    ss.wrSol(vnames=vnames,flInfo=flInfo,sfile=files[n])
    ssInfo=ss.rdFlowInfo(sfile=files[n])
    print(ssInfo)