import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.stats import rsample
from lib.myPlots import *
from lib.myclasses import *
from scipy import stats
import getpass
user=getpass.getuser()
pi=np.pi
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')
#%%
AoA=10
#Expected values
dataset='/home/'+user+'/anaconda3/pyTandem/clData/HansenClCd.dat';
aoa,cl,cd,cl0,cd0=np.loadtxt(dataset,skiprows=1,unpack=True)
#%%
f,a=getFig('Hansen Data')

a.plot(aoa,cl,label=r'$C_{L_WLE}$',color='blue',lw=2,marker='o')
a.plot(aoa,cd,label=r'$C_{D_WLE}$',color='red',lw=2,marker='o')
a.plot(aoa,cl0,label=r'$C_{L_SLE}$',color='blue',lw=2,marker='s',linestyle='--')
a.plot(aoa,cd0,label=r'$C_{D_SLE}$',color='red',lw=2,marker='s',linestyle='--')
a.axvline(x=AoA,color='black',linestyle=':',lw=2)
#%%
fcl,fcd,fcl0,fcd0=interf(aoa,cl),interf(aoa,cd),interf(aoa,cl0),interf(aoa,cd0)
print(fcl.f(AoA),fcd.f(AoA),fcl0.f(AoA),fcd0.f(AoA))