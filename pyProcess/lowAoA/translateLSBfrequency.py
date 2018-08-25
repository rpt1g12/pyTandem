import os
import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
from lib.myPlot3dOperator import *
import lib.Plot3DClasses as p3d
from scipy import stats
from scipy.interpolate import griddata
from scipy.signal import argrelextrema as locmaxmin
pi=np.pi
import getpass
user=getpass.getuser()
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')
figs=[];axs=[];hdls=[];lbls=[];lgds=[];nfig=-1;fs=18
import importlib
#%%
sim="8SLE10"
path= "/media/rperezt/082F-63FE/phd/thesis/data/lowAoA/blAnalysis/"
dataset=path+sim+".dat"
x0,delta,delta_star,theta,H,ue,uDelta,CpDelta=np.loadtxt(dataset,skiprows=1,unpack=True)
#%%
if sim=="1WLE06":
    xsep=-0.4350
elif sim=="1WLE10":
    xsep=-0.4465   
elif sim=="8SLE10":
    xsep=-0.4162
elif sim=="1SLE06":
    xsep=-0.3594
ans,dum,dum1,dum2=rsample(theta,x0,tnew=[xsep,0])
theta_sep = ans[0]
ans,dum,dum1,dum2=rsample(ue,x0,tnew=[xsep,0])
u_sep =ans[0]

f_theta=28*0.3*theta_sep/u_sep

print(f_theta,theta_sep,u_sep)

