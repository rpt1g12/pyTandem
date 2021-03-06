# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from lib.myPlots import *
pi=np.pi
sina=np.sin(np.deg2rad(20))
cosa=np.cos(np.deg2rad(20))
plt.close('all')
#%%
path='/home/rpt1g12/anaconda3/pyTandem/pgfPlots/'
#%%
fname='8WLE_pint.dat'
xm,xl,cl_wle,cd_wle=np.loadtxt(path+fname,skiprows=1,unpack=True)
cl_wle/=xl;cd_wle/=xl
#cd=cd_wle*cosa+cl_wle*sina
#cl=cl_wle*cosa-cd_wle*sina
#cd_wle=cd.copy();cl_wle=cl.copy()
for i in range(3,203,7):
    cl_wle[i-3:i+4]=cl_wle[i-3:i+4].mean()
    cd_wle[i-3:i+4]=cd_wle[i-3:i+4].mean()

#f0,a0=getFig('8WLE_pint_avg')
#a0.grid(True)
#a0.plot(xm,cl_wle,color='blue',lw=2,label='cl_wle')
#a0.plot(xm,cd_wle,color='red',lw=2,label='cd_wle')
#ax=plt.gca()
#handle,labels=ax.get_legend_handles_labels()
#legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=1)
#a0.set_xlabel(r'$x$',fontsize=18)
#a0.set_ylabel(r'$\frac{{C_{L,D}}}{{x}}$',fontsize=18)
#ax.set_ylim(-20,50)#fit(a0)
#axShow(ax)
#%%
fname='8SLE_pint.dat'
xm,xl,cl_sle,cd_sle=np.loadtxt(path+fname,skiprows=1,unpack=True)
cl_sle/=xl;cd_sle/=xl
cl_sle*=np.sqrt(1-0.3**2)
cd_sle*=np.sqrt(1-0.3**2)
#cd=cd_sle*cosa+cl_sle*sina
#cl=cl_sle*cosa-cd_sle*sina
#cd_sle=cd.copy();cl_sle=cl.copy()
for i in range(3,203,7):
    cl_sle[i-3:i+4]=cl_sle[i-3:i+4].mean()
    cd_sle[i-3:i+4]=cd_sle[i-3:i+4].mean()

#f1,a1=getFig('8SLE_pint_avg')
#ax=plt.gca()
#ax.grid(True)
#ax.plot(xm,cl_sle,color='blue',lw=2,label='cl_sle')
#ax.plot(xm,cd_sle,color='red',lw=2,label='cd_sle')
#ax=plt.gca()
#handle,labels=ax.get_legend_handles_labels()
#legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=1)
#ax.set_xlabel(r'$x$',fontsize=18)
#ax.set_ylabel(r'$\frac{{C_{L,D}}}{{x}}$',fontsize=18)
#ax.set_ylim(-20,50)#fit(ax)
#axShow(ax)
Delta_CL=cl_wle*xl-cl_sle*xl
Delta_CD=cd_wle*xl-cd_sle*xl
#%%
fname='8WLE_circulation.dat'
z0,c_wle=np.loadtxt(path+fname,skiprows=1,unpack=True)
fname='8SLE_circulation.dat'
z1,c_sle=np.loadtxt(path+fname,skiprows=1,unpack=True)
c_sle*=np.sqrt(1-0.3**2)
c_sle*=(2/(0.3*0.88));c_wle*=(2/(0.3*0.88))
#f2,a2=getFig(fname)
#ax=plt.gca()
#ax.grid(True)
#ax.plot(z0,c_wle,color='blue',lw=2,label='wle')
#ax.plot(z0,c_sle,color='red',lw=2,label='sle')
#ax.axhline(y=c_wle.mean(),color='blue',lw=1,label='wle_avg')
#ax.axhline(y=c_sle.mean(),color='red',lw=1,label='sle_avg')
#ax.set_xlabel(r'$z/L_{LE}$',fontsize=18)
#ax.set_ylabel(r'$2\Gamma/U_\infty$',fontsize=18)
#handle,labels=ax.get_legend_handles_labels()
#legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=1)
#fit(ax)
#%%
f3,a3=getFig('Delta_p')
ax=plt.gca()
ax.grid(True)
delta_cl=cl_wle-cl_sle
delta_cd=cd_wle-cd_sle

xp=np.linspace(-0.515,0.5,len(xm))
ax.plot(xp,delta_cl,color='blue',lw=2,label='cl')
ax.plot(xp,delta_cd,color='red',lw=2,label='cd')
ax.set_xlim(-0.5,0.5)
ax=plt.gca()
handle,labels=ax.get_legend_handles_labels()
legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=1)
ax.set_xlabel(r'$x$',fontsize=18)
ax.set_ylabel(r'$\Delta$',fontsize=18)
axShow(ax)
#%%
f4,a4=getFig('Delta_Gamma')
ax=plt.gca()
ax.grid(True)
delta_g=c_wle-c_sle
ax.plot(z0,delta_g,color='blue',lw=2,label='Gamma')
ax=plt.gca()
handle,labels=ax.get_legend_handles_labels()
legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=1)
ax.set_xlabel(r'$z/L_{LE}$',fontsize=18)
ax.set_ylabel(r'$\Delta\Gamma$',fontsize=18)
fit(ax)
#%%
save=False

if (save==True):
    savePlotFile(ax=a4)