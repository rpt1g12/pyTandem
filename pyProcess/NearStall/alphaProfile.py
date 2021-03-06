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
figs=[];axs=[];hdls=[];lbls=[];lgds=[];nfig=-1;fs=18
#%%
save=True
spath='/home/rpt1g12/Documents/thesis/data/nearStall/aoaTime/'
#%%
sin=np.sin
deltaAoA=10
T=120
tup=np.linspace(0,120,240)
tfix=np.linspace(120,150,60)
tdown=np.linspace(150,270,240)
aup=10+10*sin(0.5*pi*tup/T)*sin(0.5*pi*tup/T)
afix=np.ones(60)*20
adown=20-10*sin(0.5*pi*(tdown-150)/T)*sin(0.5*pi*(tdown-150)/T)
#%%
f,a=getFig('aoaTime');figs.append(f),axs.append(a);nfig+=1
daup=fcbFD(aup,tup[1]-tup[0])
dadown=fcbFD(adown,tdown[1]-tdown[0])
dafix=fcbFD(afix,tfix[1]-tfix[0])
a.plot(tup,aup,color='blue',label='up')
a.plot(tfix,afix,color='red',label='fix')
a.plot(tdown,adown,color='green',label='down')
a.plot(tup,daup,color='blue',linestyle='--',label='dup')
a.plot(tfix,dafix,color='red',linestyle='--',label='dfix')
a.plot(tdown,dadown,color='green',linestyle='--',label='ddown')
#%%
nfigs=len(figs)
for i in range(nfigs):
    fit(axs[i])
    hdl,lbl,lgd=getLabels(ax=axs[i],ncol=3,fontsize=15)
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
#%% Saving
    if save:
        savePlotFile(path=spath,ax=axs[i])