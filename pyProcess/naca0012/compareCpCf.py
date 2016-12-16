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
opt=1; plus=False; save=False; dudy=True
path='/home/'+user+'/Desktop/post/naca0012G1/ss003/'

files=p3d.getFileNames(path=path)

vnames=['cf','twx','twy','twz','cp']

fl=p3d.flow(path,"grid.xyz",'solTavgCfTwCp.qa')
fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames)
up=fl.blk[4]
bot=fl.blk[1]

cp_up=up.var['cp'].avgDir()[:,0]
cf_up=up.var['cf'].avgDir()[:,0]
#cf_up=2*(up.var['twx'].avgDir()[:,0])/(0.4**2)
x_up=up.var['x'].getValues()[:,0,0]

if plus:    
    y_up=up.var['y'].getValues()[:,0,0]
    y1_up=up.var['y'].getValues()[:,1,0]
    x1_up=up.var['x'].getValues()[:,1,0]
    dy_up=np.sqrt((x1_up-x_up)**2+(y1_up-y_up)**2)

cp_bot=bot.var['cp'].avgDir()[:,-1]
cf_bot=bot.var['cf'].avgDir()[:,-1]
#cf_bot=2*(bot.var['twx'].avgDir()[:,-1])/(0.4**2)
x_bot=bot.var['x'].getValues()[:,-1,0]


#Load Jones2008 data
xcp_up,jcp_up=np.loadtxt('cpData/Jones2008_up.dat',skiprows=1,unpack=True)
xcp_bot,jcp_bot=np.loadtxt('cpData/Jones2008_bot.dat',skiprows=1,unpack=True)
xcf_up,jcf_up=np.loadtxt('cfData/Jones2008_up.dat',skiprows=1,unpack=True)
xcf_bot,jcf_bot=np.loadtxt('cfData/Jones2008_bot.dat',skiprows=1,unpack=True)
sxcf_up,scf_up=np.loadtxt('pgfPlots/Cf_dudy.dat',skiprows=1,unpack=True)


if opt==0:    
    cvar=r'$C_p$'
    var_up=cp_up
    var_bot=cp_bot
    jvar_up=jcp_up
    jx_up=xcp_up-0.5
    jvar_bot=jcp_bot
    jx_bot=xcp_bot-0.5
    title='Cp_comparison'
else:
    cvar=r'$C_f$'
    var_up=cf_up
    var_bot=cf_bot
    jvar_up=jcf_up
    jx_up=xcf_up-0.5
    jvar_bot=jcf_bot
    jx_bot=xcf_bot-0.5
    title='Cf_comparison'
    ns=np.where(var_up<0)[0][0]
    nr=np.where(var_up<0)[0][-1]
    xs=x_up[ns];xr=x_up[nr]

f,a=getFig(title)
    
a.plot(x_up,var_up,color='blue',lw=2,marker='o',label='LES',markevery=5)
a.plot(jx_up,jvar_up,color='red',lw=2,label='DNS')
#a.plot(x_bot,var_bot,color='blue',lw=2,marker='o',markevery=5)
#a.plot(jx_bot,jvar_bot,color='red',lw=2)



a.set_xlabel(r'$x$',fontsize=20)
a.set_ylabel(cvar,fontsize=20)


if opt==1:
    a.axhline(y=0,lw=2,color='black')
    s=r'$x_{sep}=$'+'{:-5.3f}'.format(xs)+'\t'+r'$x_{rea}=$'+'{:-5.3f}'.format(xr)
    a.text(0.25,-0.15,s,ha='center',va='bottom',transform=a.transAxes)
    a.set_xlim(-0.5,0.5)
    a.set_ylim(-0.01,0.02)
    if plus:
        f2,a2=getFig('Wall Units')
        y_plus=125000*0.5*0.4*np.sqrt(2*abs(var_up))*dy_up
        a2.plot(x_up,y_plus)
    if dudy:
        a.plot(sxcf_up,scf_up,lw=2,color='green',label='2ndOrder')
else:    
    fit(a,(0,0.05))
    a.invert_yaxis()
    
h,lbl,lgd=getLabels(ax=a)
f.tight_layout()

if save:
    savePlotFile(ax=a,vary=['g1_up','dns_up','g1_bot','dns_bot'])