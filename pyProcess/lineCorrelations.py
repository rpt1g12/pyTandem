import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
#%%
save=True; scale=False; step=True;mod='full'

for case in range(3):
    if (case==0): sim='A4A15W11AoA20'
    if (case==1):sim='A00W11AoA20'
    if (case==2):sim='A15W11AoA20'
    fig=plt.figure()
    fig.canvas.set_window_title('Correlation '+sim)
    ax=fig.add_subplot(111)
    fig2=plt.figure()
    fig2.canvas.set_window_title('Correlation Contour '+sim)
    ax2=fig2.add_subplot(111)
    path='lines/'+sim+'/'
    lst=range(1,11);n=-1;lvl=np.linspace(0,1,21)
    for i in lst:
        n+=1
        file=path+str(i)+'+0.15.dat'
        uu,vv,ww=np.loadtxt(file,skiprows=1,unpack=True,usecols=range(5,8))
        p=np.loadtxt(file,skiprows=1,unpack=True,usecols=[4])
        nsam=len(uu)
        if (i==lst[0]):
            muu=np.zeros((nsam,len(lst)))
            mvv=muu.copy()
            mww=mvv.copy()
        p=p-p.mean()
        x=np.loadtxt(file,skiprows=1,unpack=True,usecols=range(11,14))
        z=(x[2]-x[2][0])/0.11
        muu[:,n]=acorrp(uu)
        mvv[:,n]=acorrp(vv)
        mww[:,n]=acorrp(ww)
        ax.plot(z,muu[:,n],'b')
        ax.plot(z,mvv[:,n],'r')
        ax.plot(z,mww[:,n],'g')
        #ax.plot(z,acorrp(p),'k')
    

    #%%
    ax.set_ylim(0,1)
    l=np.linspace(-0.4,0.5,10)
    X,Y=np.meshgrid(l,z)
    ax2.contourf(X,Y,mww,levels=lvl,cmap='Greys_r')
