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

folder='pitchUp/';sclcd='scl'
osim=10;aoa=12;iaoah=aoa;freq=0.20;
ttotal=10*0.3;ns=1024
nw=4;ovlp=0.5;

sina=np.sin(np.deg2rad(aoa))
cosa=np.cos(np.deg2rad(aoa))

#Expected values
dataset='/home/'+user+'/anaconda3/pyTandem/clData/HansenClCd.dat';
aoah,hcl,hcd=np.loadtxt(dataset,skiprows=1,unpack=True)
clh=np.interp(aoa,aoah,hcl)
cdh=np.interp(aoa,aoah,hcd)
print(clh,cdh)
#clh=0.35;cdh=0.035

if(osim==1):
    sim='G1'
    prdl=0
    name='G1'
    folder='naca0012/'
if(osim==10):
    sim='G1'
    prdl=0
    name='G1'
    folder='pitchUp/'
elif(osim==2):
    sim='G2'
    prdl=0
    name='G2'
    folder='naca0012/'
elif(osim==3):
    sim='G3'
    prdl=0
    name='G3'
    folder='naca0012/'
    
dataset='/home/'+user+'/anaconda3/pyTandem/clData/'+folder+'clcdPresVisc_'+sim+'.dat';
   
n,tin,clinp,clinv,cdinp,cdinv=np.loadtxt(dataset,skiprows=1,unpack=True)

#%%

fTime,axTime=getFig('Time '+sim)

tmin=tin[-1]-(ttotal/0.3)       #in jfm paper its -(25/0.3)
n0=np.where(tin>tmin)[0][0]
t0=tin[n0:]#*0.3#-tin[n0];
if prdl==1:
    fctr=np.sqrt(1-0.3**2)
else:
    fctr=1

cl0p=clinp[n0:]*fctr;
cl0v=clinv[n0:]*fctr;
cd0p=cdinp[n0:]*fctr;
cd0v=cdinv[n0:]*fctr;


clnp,tn,nsam,fsam=rsample(x=cl0p,t=t0,verbose=True,nsample=ns)
cdnp,tn,nsam,fsam=rsample(x=cd0p,t=t0,nsample=ns)
clnv,tn,nsam,fsam=rsample(x=cl0v,t=t0,nsample=ns)
cdnv,tn,nsam,fsam=rsample(x=cd0v,t=t0,nsample=ns)


clMp=clnp.mean()
cdMp=cdnp.mean()
clMv=clnv.mean()
cdMv=cdnv.mean()

cln=clnp+clnv;cdn=cdnp+cdnv
clM,cdM=cln.mean(),cdn.mean()
cls=np.sqrt(np.var(cln))
cds=np.sqrt(np.var(cdn))

#%%
    
axTime.plot(tn,cln,'b',linewidth=2,label='Cl'+name)
axTime.plot(tn,cdn,'r',linewidth=2,label='Cd'+name)
axTime.plot([tn[0],tn[-1]],[clh,clh],'k',linewidth=2,label=r'DNS $\alpha=$'+str(iaoah))
axTime.plot([tn[0],tn[-1]],[cdh,cdh],'k',linewidth=2)


axTime.set_xlabel(r'$t^*$',fontsize=20)
axTime.set_ylabel(r'$C_l$',fontsize=20)
handle,labels=axTime.get_legend_handles_labels()
legend=axTime.legend(handle,labels,bbox_to_anchor=(1,1),ncol=3)
scl=r'$\sigma_L=$'+str(np.around(cls,3))+'\t'+r'$\bar{C_L}=$'+str(np.around(clM,3))
scd=r'$\sigma_D=$'+str(np.around(cds,3))+'\t'+r'$\bar{C_D}=$'+str(np.around(cdM,3))
axTime.text(0.25,-0.15,scl,ha='center',va='bottom',transform=axTime.transAxes)
axTime.text(0.75,-0.15,scd,ha='center',va='bottom',transform=axTime.transAxes)

fit(axTime,spg=(0,0.4))