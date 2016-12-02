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
from scipy.interpolate import splrep,splev
pii=np.pi
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')
cosaoa=np.cos(np.deg2rad(20))
sinaoa=np.sin(np.deg2rad(20))
import importlib
#%%
importlib.reload(p3d)
#%%
user='rpt1g12'
path='/home/'+user+'/Desktop/post/8A15W11AoA20/aerofoil/'
prdtl=1#np.sqrt(1-0.3**2)
fl=p3d.flow(path,"grid.xyz","solTa+n+p.q")

#vnames=['Cf','twx','twy','twz','Cp']
vnames=['a','nx','ny','nz','p']

fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames=vnames)
#%%
lvl=np.linspace(-1,1,51)
#upfig,upax=getFig('Upper_Surface')
bl0=0 #Block number
#upax.contourf(fl.blk[bl].var['x'].val[:,0,:],fl.blk[bl].var['z'].val[:,0,:],fl.blk[bl].var['ny'].val[:,0,:],levels=lvl,cmap='bwr')
#%%
if bl0==0:
    lst=[0,1]
else:
    lst=[1,0]
fsum=0
for bl in lst:
    x=fl.blk[bl].var['x'].getValues()
    y=fl.blk[bl].var['y'].getValues()
    z=fl.blk[bl].var['z'].getValues()
    p=fl.blk[bl].var['p'].getValues()
    nx=fl.blk[bl].var['nx'].getValues()
    ny=fl.blk[bl].var['ny'].getValues()
    a=fl.blk[bl].var['a'].getValues()
    #%%
    #x=fl.blk[bl].var['x'].getValues(True)
    #y=fl.blk[bl].var['y'].getValues(True)
    #z=fl.blk[bl].var['z'].getValues(True)
    #p=fl.blk[bl].var['p'].getValues(True)
    #nx=fl.blk[bl].var['nx'].getValues(True)
    #ny=fl.blk[bl].var['ny'].getValues(True)
    #a=fl.blk[bl].var['a'].getValues(True)
    #%%
#    xr=cosaoa*x+sinaoa*y
#    yr=cosaoa*y-sinaoa*x
#    nxr=cosaoa*nx+sinaoa*ny
#    nyr=cosaoa*ny-sinaoa*nx
#    x=xr.copy()
#    y=yr.copy()
    #%%
    fy=p*ny
    fx=p*nx
    f=(fy*cosaoa-fx*sinaoa)*prdtl
    #f=(fx*cosaoa+fy*sinaoa)*prdtl
    #f=fy
    pdyn=0.88*0.5*0.3**2
    fsum+=(f*a).sum()
    #%%
    #xlim=(x.min(),x.max());zlim=(z.min(),z.max());
    #xi,zi=np.mgrid[xlim[0]:xlim[1]:413j,zlim[0]:zlim[1]:201j]
    #isize=list(xi.shape)
    #%%
#    if bl==bl0:
#        ipoints=x[:,0,6]
    #ipoints=-0.45*np.cos(np.linspace(0,pii,801))
    #ipoints=np.linspace(-0.5149,-0.485,11)
    ipoints=np.linspace(-0.485,0.5,301)
    #print(ipoints.min())
    pi=np.zeros((ipoints.shape[0],201))
    nxi=pi.copy()
    nyi=pi.copy()
    ai=pi.copy()
    fi=pi.copy()
    xi=pi.copy()
    zi=pi.copy()
    method='linear';fval=np.nan
    for kk in range(x.shape[2]):
        points=x[:,0,kk]
        values=f[:,0,kk]
        ival=griddata(points,values,ipoints,method=method,fill_value=fval)
        fi[:,kk]=ival.copy()
#        values=nx[:,0,kk]
#        ival=griddata(points,values,ipoints,method=method,fill_value=fval)
#        nxi[:,kk]=ival.copy()
#        values=ny[:,0,kk]
#        ival=griddata(points,values,ipoints,method=method,fill_value=fval)
#        nyi[:,kk]=ival.copy()
#        values=a[:,0,kk]
#        ival=griddata(points,values,ipoints,method=method,fill_value=fval)
#        ai[:,kk]=ival.copy()
#        fi[:,kk]=(pi[:,kk]*nyi[:,kk]*ai[:,kk]*cosaoa-pi[:,kk]*nxi[:,kk]*ai[:,kk]*sinaoa)
        xi[:,kk]=ipoints
        zi[:,kk]=z[0,0,kk]
        
    #%%
    kk=6
    figs,axs=getFig('block {:d}'.format(bl))
    axs.plot(x[:,0,kk],f[:,0,kk],color='blue',marker='o')
    axs.plot(ipoints,fi[:,kk],color='red',marker='s')
    #%%
    #method='linear';fval=0
    #points=np.array(np.transpose([x,z])) 
    #values=a
    #ipoints=np.array(np.transpose([np.reshape(xi,xi.size),np.reshape(zi,zi.size)]))
    #ival=griddata(points,values,ipoints,method=method,fill_value=fval)
    #pi=np.reshape(ival,isize)
    #%%
#    
#    iupfig,iupax=getFig('iUpper_Surface')
#    iupax.contourf(xi,zi,pi)
    #%%
    if (bl==bl0):
        fi_avg=np.zeros((xi.shape[0],2))
        ncount=fi_avg.copy()
    for i in range(xi.shape[0]):
        for k in range(xi.shape[1]):
            rdum=fi[i,k]
            if (not np.isnan(rdum)):
                ncount[i,bl]+=1
                fi_avg[i,bl]+=fi[i,k]
    
    
fi_avg_t=fi_avg[:,0]+fi_avg[:,1]
fctr=0.88/200
print(fctr*np.trapz(fi_avg_t,xi[:,0])/pdyn)
print(fsum/pdyn)
    
#%%
figp,axp=getFig('plot')
axp.plot(xi[:,0],fi_avg_t,color='black',lw=2)
#axp.plot(xi[:,0],pi_avg[:,0],color='red',lw=2)
#axp.plot(xi[:,0],pi_avg[:,1],color='blue',lw=2)
#axp.plot(xi,pi)
#%%
#sle=pi_avg_t.copy()
#wle=pi_avg_t.copy()
#delta=wle-sle
#fdelta,adelta=getFig('delta_cd')
#adelta.grid(True)
#adelta.plot(ipoints,delta,color='black',lw=2,label='delta')
#adelta.plot(ipoints,wle,color='blue',lw=2,label='wle')
#adelta.plot(ipoints,sle,color='red',lw=2,label='sle')
#handle,labels=adelta.get_legend_handles_labels()
#legend=adelta.legend(handle,labels,bbox_to_anchor=(1,1),ncol=1)
##%%
#save=True
#name=fdelta.canvas.get_window_title()
#ax=adelta
#if (save==True):
#    plt.sca(ax)
#    path='pgfPlots/'
#    savePlotFile(path=path+name+'.dat',vary=labels)