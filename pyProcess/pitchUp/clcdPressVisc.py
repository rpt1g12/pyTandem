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
osim=10;aoa=13;iaoah=aoa;freq=0.20;
ttotal=150*0.3;ns=1024
nw=4;ovlp=0.5;

if(osim==10):
    sim='G1'
    prdl=0
    name='G1'
    folder='pitchUp/'

sina=np.sin(np.deg2rad(aoa))
cosa=np.cos(np.deg2rad(aoa))

#Expected values
dataset='/home/'+user+'/anaconda3/pyTandem/clData/HansenClCd.dat';
aoah,hcl,hcd=np.loadtxt(dataset,skiprows=1,unpack=True)
clh=np.interp(aoa,aoah,hcl)
cdh=np.interp(aoa,aoah,hcd)

print(clh,cdh)

    
dataset='/home/'+user+'/anaconda3/pyTandem/clData/'+folder+'clcdPresVisc_'+sim+'.dat';
   
n,tin,clinp,clinv,cdinp,cdinv=np.loadtxt(dataset,skiprows=1,unpack=True)

dataset='/home/'+user+'/anaconda3/pyTandem/clData/'+folder+'clcdAoA_'+sim+'.dat';
   
_,_,clin,cdin,AoAin,_=np.loadtxt(dataset,skiprows=1,unpack=True)

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


clnp,tn,nsam,fsam=rsample(x=cl0p,t=t0,verbose=True,force=True)
cdnp,tn,nsam,fsam=rsample(x=cd0p,t=t0,force=True)
clnv,tn,nsam,fsam=rsample(x=cl0v,t=t0,force=True)
cdnv,tn,nsam,fsam=rsample(x=cd0v,t=t0,force=True)
AoA,_,_,_=rsample(x=AoAin[n0:],t=t0,force=True)

clnh,cdnh=np.zeros_like(AoA),np.zeros_like(AoA)
for i in range(len(AoA)):
    a=AoA[i]
    clnh[i]=np.interp(a,aoah,hcl)
    cdnh[i]=np.interp(a,aoah,hcd)



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
axTime.plot(tn,clnh,color='blue',lw=2,linestyle='--',label=r'Exp')
axTime.plot(tn,cdnh,color='red',lw=2,linestyle='--')

axTime.set_xlabel(r'$t^*$',fontsize=20)
axTime.set_ylabel(r'$C_l$',fontsize=20)
handle,labels=axTime.get_legend_handles_labels()
legend=axTime.legend(handle,labels,bbox_to_anchor=(1,1),ncol=3)
scl=r'$\sigma_L=$'+str(np.around(cls,3))+'\t'+r'$\bar{C_L}=$'+str(np.around(clM,3))
scd=r'$\sigma_D=$'+str(np.around(cds,3))+'\t'+r'$\bar{C_D}=$'+str(np.around(cdM,3))
axTime.text(0.25,-0.15,scl,ha='center',va='bottom',transform=axTime.transAxes)
axTime.text(0.75,-0.15,scd,ha='center',va='bottom',transform=axTime.transAxes)

fit(axTime,spg=(0,0.4))
#%%

AoA0=15.6;AoA1=16.0
ns=np.where(AoA>AoA0)[0][0]-1
ne=np.where(AoA<AoA1)[0][-1]+1
claoa,aoa_2,nsam_aoa,faoa=rsample(cln[ns:ne],AoA[ns:ne],nsample=16)
dcldaoa=fcbFD(claoa,1/faoa)

f,a=getFig('dCl/dAoA')
a.plot(aoa_2,dcldaoa)