import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.stats import rsample
from lib.myPlots import *
from scipy import stats
import getpass
user=getpass.getuser()
pi=np.pi
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')
#%%
sin20=np.sin(np.deg2rad(20))

folder='naca0012/';sclcd='scl'
osim=2;aoa=5;iaoah=20;freq=0.20;ttotal=20*0.3;ns=1024
nw=4;ovlp=0.5;

sina=np.sin(np.deg2rad(aoa))
cosa=np.cos(np.deg2rad(aoa))

#Expected values
clh=0.621;cdh=0.0358

if(osim==1):
    sim='G1'
    prdl=0
    linreg=0
    name='G1'
elif(osim==2):
    sim='G2'
    prdl=0
    linreg=0
    name='G2'
    
dataset='/home/'+user+'/anaconda3/pyTandem/clData/'+folder+sim+'.dat';
   
n,tin,clin,cdin,taoa,tmach=np.loadtxt(dataset,skiprows=1,unpack=True)

#%%

fTime,axTime=getFig('Time '+sim)

show=False;showAoA=False;getAvg=True


tmin=tin[-1]-(ttotal/0.3)       #in jfm paper its -(25/0.3)
n0=np.where(tin>tmin)[0][0]
t0=tin[n0:]#*0.3#-tin[n0];
if prdl==1:
    fctr=np.sqrt(1-0.3**2)
else:
    fctr=1
cl0=clin[n0:]*fctr;
cd0=cdin[n0:]*fctr
taoa0=taoa[n0:]
cln,tn,nsam,fsam=rsample(x=cl0,t=t0,verbose=True,nsample=ns)
cdn,tn,nsam,fsam=rsample(x=cd0,t=t0,nsample=ns)
dtaoa,tn,nsam,fsam=rsample(x=taoa0,t=t0,nsample=ns)

clM=cln.mean()
cdM=cdn.mean()

if linreg==1:
    #Linear regression
    x=tn.copy();y=cln.copy()
    z=np.polyfit(x,y,3)    
    p_cl=np.poly1d(z)

    y=cdn.copy()
    z=np.polyfit(x,y,3)
    p_cd=np.poly1d(z)
    lin_cl=p_cl(x)
    lin_cd=p_cd(x)
else:
    lin_cl=np.ones(len(tn))*clM
    lin_cd=np.ones(len(tn))*cdM

s=np.sqrt(np.var(cln-lin_cl))
cds=np.sqrt(np.var(cdn-lin_cd))

#%%
if (show):
    axTime.plot(tin[:n0],clin[:n0],'b');lClin=plt.gca().lines[-1]
    axTime.plot(tin[:n0],cdin[:n0],'r');lCdin=plt.gca().lines[-1]
    axTime.axvline(x=tmin,color='green',linewidth=2,linestyle='--')
    
if (showAoA):
    ax2=axTime.twinx()
    ax2.plot(tn,dtaoa,color='lime',linewidth=2,label=r'$\alpha$')
    ax2.set_ylabel(r'$\alpha(t^*)$',fontsize=20)
    
axTime.plot(tn,cln,'b',linewidth=2,label='Cl'+name)
axTime.plot(tn,cdn,'r',linewidth=2,label='Cd'+name)
axTime.plot([tn[0],tn[-1]],[clh,clh],'k',linewidth=2,label=r'DNS $\alpha=$'+str(iaoah))
axTime.plot([tn[0],tn[-1]],[cdh,cdh],'k',linewidth=2)

if (getAvg):
    axTime.plot(tn,lin_cl,'g-.',linewidth=2,label='avg');lClnm=plt.gca().lines[-1]
    axTime.plot(tn,lin_cl-s,'r-.',linewidth=2,label='s0');lCls0=plt.gca().lines[-1]
    axTime.plot(tn,lin_cl+s,'r-.',linewidth=2,label='s1');lCls1=plt.gca().lines[-1]
    axTime.plot(tn,lin_cd,'g-.',linewidth=2,label='avg');lCdnm=plt.gca().lines[-1]
    axTime.plot(tn,lin_cd-cds,'r-.',linewidth=2,label='s0');lCds0=plt.gca().lines[-1]
    axTime.plot(tn,lin_cd+cds,'r-.',linewidth=2,label='s1');lCds1=plt.gca().lines[-1]


axTime.set_xlabel(r'$t^*$',fontsize=20)
axTime.set_ylabel(r'$C_l$',fontsize=20)
handle,labels=axTime.get_legend_handles_labels()
legend=axTime.legend(handle,labels,bbox_to_anchor=(1,1),ncol=3)
scl=r'$\sigma_L=$'+str(np.around(s,3))+'\t'+r'$\bar{C_L}=$'+str(np.around(clM,3))
scd=r'$\sigma_D=$'+str(np.around(cds,3))+'\t'+r'$\bar{C_D}=$'+str(np.around(cdM,3))
axTime.text(0.25,-0.15,scl,ha='center',va='bottom',transform=axTime.transAxes)
axTime.text(0.75,-0.15,scd,ha='center',va='bottom',transform=axTime.transAxes)

fit(axTime,spg=(0,0.4))

if (showAoA):
    fit(ax2,spg=(0,0))
    f2,a2=getFig('Cl & Cd vs alpha')
    a2.plot(dtaoa,cln,lw=2,color='blue',label='Cl')
    a2.plot(dtaoa,cdn,lw=2,color='red',label='Cd')
    a2.set_xlabel(r'$\alpha$',fontsize=20)
    fit(a2)