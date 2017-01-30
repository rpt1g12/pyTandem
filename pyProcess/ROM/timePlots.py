import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
from scipy import stats
import getpass
user=getpass.getuser()
pi=np.pi
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')
#%%

# Font Size
fs=18
folder='rom/4A15W11AoA10/upperSurf/'
sim="POD_time"
dataset='/home/'+user+'/anaconda3/pyTandem/'+folder+sim+'.dat';
t0=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=(range(1)))
modes0=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=(range(1,6)))

ts,te=t0[0],t0[-2]

sim="POD_SngV"
dataset='/home/'+user+'/anaconda3/pyTandem/'+folder+sim+'.dat';
nmode,sgma=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=[0,1])

fsgma,asgma=getFig('SingularValues')
asgma.semilogy(nmode[:5],sgma[:5],lw=2,marker='o',label='sgma')
asgma.set_ylabel(r'$\sigma_i$',fontsize=fs)
asgma.set_xlabel(r'mode #',fontsize=fs)
savePlotFile(ax=asgma)

sizes=modes0.shape

ns=np.where((t0>=ts))[0][0]
ne=np.where((t0>te))[0][0]
trom=t0[ns:ne]
modes=modes0[:,ns:ne]

fmod,amod=getFig('Modes')
fmod.set_size_inches(19.2/2,10.8/2,forward=True)
amod.set_xlabel(r'$t^*$',fontsize=fs)
fmod.tight_layout()

for m in range(1,sizes[0]):
    amod.plot(trom,sgma[m]*modes[m,:],lw=2,label='mode{:d}'.format(m))

fit(amod,spg=(0,0.1))

handle1,labels1,legend1=getLabels(ax=amod,ncol=5,fontsize=12,loc='upper right',hspace=0.2)
   
line1=amod.axvline(x=trom[0],color='black',linewidth=2,linestyle='--')

#amod.set_xlabel(r'$t^*$',fontsize=20)
amod.set_ylabel(r'POD Coefficient',fontsize=fs)


    
#%%
folder='clData/aoa10/'
sim="4A15W11AoA10"
dataset='/home/'+user+'/anaconda3/pyTandem/'+folder+sim+'.dat';
ts,te=45,314.98
n,t0,clin,cdin,taoa,tmach=np.loadtxt(dataset,skiprows=1,unpack=True)

ns=np.where((t0>=ts))[0][0]
ne=np.where((t0>te))[0][0]
tcl=t0[ns:ne]
cl=clin[ns:ne];cd=cdin[ns:ne]
aoa=taoa[ns:ne]

fcl,acl=getFig('cl')
#acl=addAx(fmod,(212))
acl.plot(tcl,cl,lw=2,color='blue',label=r'$C_L$')
acl.plot(tcl,cd,lw=2,color='red',label=r'$C_D$')

fit(acl)

handle2,labels2,legend2=getLabels(ax=acl,ncol=2,fontsize=15)

line2=acl.axvline(x=tcl[0],color='black',linewidth=2,linestyle='--') 
acl.set_xlabel(r'$t^*$',fontsize=fs)
acl.set_ylabel(r'$C_L\, &\, C_D$',fontsize=fs)


#faoa,aaoa=getFig('aoa')
#aaoa.plot(tcl,aoa,lw=2,color='blue',label=r'$\alpha$')


#line3=aaoa.axvline(x=tcl[0],color='black',linewidth=2,linestyle='--') 
#aaoa.set_xlabel(r'$t^*$',fontsize=20)
#aaoa.set_ylabel(r'$\alpha$',fontsize=20)

for a in [amod,acl]:
    axShow(a)
    a.set_xlim(ts,te)

#for i in range(int((te-ts)/0.125)):
#    print(i)
#    tt=ts+0.125*i
#    adum=[tt,tt]
#    line1.set_xdata(adum)
#    line2.set_xdata(adum)
##    line3.set_xdata(adum)
#    axShow(amod)
##    axShow(acl)
#    plt.pause(0.01)
#    #fmod.savefig('/home/rpt1g12/kdenlive/rom/mod/mod.'+'{:04d}'.format(i)+'.png',dpi=300)
##    fcl.savefig('/home/rpt1g12/kdenlive/rom/cl/cl.'+'{:04d}'.format(i)+'.png')
##    faoa.savefig('/home/rpt1g12/kdenlive/rom/aoa/aoa.'+'{:04d}'.format(i)+'.png')
