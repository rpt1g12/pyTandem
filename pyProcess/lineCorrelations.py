import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
from lib.myIO import *
#%%
save=True; scale=False; step=True;mod='full';rmMean=False
plt.close('all')
for case in [0,1,2]:
    if (case==0): sim='A4A15W11AoA20'
    if (case==1):sim='A00W11AoA20'
    if (case==2):sim='A15W11AoA20'
    fig=plt.figure()
    fig.canvas.set_window_title('Correlation '+sim)
    ax=fig.add_subplot(111)
    fig2=plt.figure()
    fig2.canvas.set_window_title('Rxx_'+sim)
    ax2=fig2.add_subplot(111)
    fig3=plt.figure()
    fig3.canvas.set_window_title('Ryy_'+sim)
    ax3=fig3.add_subplot(111)
    fig4=plt.figure()
    fig4.canvas.set_window_title('Rzz_'+sim)
    ax4=fig4.add_subplot(111)
    path='lines/'+sim+'/'
    lst=range(0,31);n=-1;
    for i in lst:
        n+=1
        file=path+str(i)+'+0.25.dat'
        uu,vv,ww=np.loadtxt(file,skiprows=1,unpack=True,usecols=range(1,4))
        p=np.loadtxt(file,skiprows=1,unpack=True,usecols=[4])
        nsam=len(uu)
        if (rmMean):
            p=p-p.mean()
            uu=uu-uu.mean()
            vv=vv-vv.mean()
            ww=ww-ww.mean()
        x=np.loadtxt(file,skiprows=1,unpack=True,usecols=range(14,17))
        z=(x[2]-x[2][0])/0.11
        if (i==lst[0]):
            muu=np.zeros((nsam,len(lst)))
            mvv=muu.copy()
            mww=mvv.copy()
            x0=x[0][0]
        #plt.figure();plt.plot(z,uu);plt.show()
        muu[:,n]=acorrp(uu)
        mvv[:,n]=acorrp(vv)
        mww[:,n]=acorrp(ww)
        ax.plot(z,muu[:,n],'b')
        ax.plot(z,mvv[:,n],'r')
        ax.plot(z,mww[:,n],'g')
        #ax.plot(z,acorrp(p),'k')
    xn=x[0][-1]

    #%%
    ax.set_ylim(-1,1)
    if ((muu.min()<0) or (mvv.min()<0) or (mww.min()<0)):
        lvl=np.linspace(-1,1,41)
        cm='bwr'
    else:
        lvl=np.linspace(0,1,41)
        cm='inferno'
    l=np.linspace(x0,xn,len(lst));
    X,Y=np.meshgrid(l,z)
    ax2.contourf(X,Y,muu,levels=lvl,cmap=cm)
    ax3.contourf(X,Y,mvv,levels=lvl,cmap=cm)
    img=ax4.contourf(X,Y,mww,levels=lvl,cmap=cm)

    for a in [ax2,ax3,ax4]:    
        a.axes.get_yaxis().set_visible(False)
        a.axes.get_xaxis().set_visible(False)
    for f in [fig2,fig3,fig4]:
        f.savefig(path+f.canvas.get_window_title()+'.pdf')
#%%
    wrtContour(path+'muu.dat',X,Y,muu,['x','z','Ruu'])
    wrtContour(path+'mvv.dat',X,Y,mvv,['x','z','Rvv'])
    wrtContour(path+'mww.dat',X,Y,mww,['x','z','Rww'])
    
#%%
#plt.close('all')
#ax4.set_visible(False)
#cbar=fig4.colorbar(img,orientation='horizontal')
#cbar.ax.yaxis.set_visible(False)
#cbar.ax.xaxis.set_visible(False)
#fig4.savefig(path+cm+'.pdf')