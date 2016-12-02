from scipy.signal import coherence as msc
from scipy.signal import welch as psdw
from scipy import signal as sgl
import numpy as np
from lib.stats import *
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm
from lib.matplotlib2tikz import save as tikz_save
from os import system
from os import getcwd
#%%
cwd=getcwd()
for nprob in range(1,11):
    plt.clf()
    cpre='strprb'
    cprob=str(nprob).zfill(2)
    savepath='/home/rpt1g12/Desktop/circ/'+cpre+cprob+'.png'
    x=np.loadtxt(cwd+'/signalData/clean/AoA20/circLine/straight/circ0'+cprob+'.dat',skiprows=1,unpack=True,usecols=range(102))
    lzee=x.shape[0];lzes=1;lze=lzee-lzes
    n=150
    for i in range(lzes,lzee):
        x[i,:]=x[i,:]-avg(x[i,:],x[0,:])
        if (i==lzes):
            out1,out2,nsam,fsam=rsample(x[i,:],x[0,:],n,verbose=True)
            z=np.zeros((lzee,len(out1)))
            z[i,:]=out1[:];z[0,:]=out2[:]
        else:
            z[i,:],out2,nsam,fsam=rsample(x[i,:],x[0,:],n)
    lzee=z.shape[0];lzes=1;lze=lzee-lzes

    plot=0
    ovlp=0.3;nwin=2;

    nseg,novlp,ntot,fmax,fmin=defWindows(z,nwin,ovlp,plot)

    psd=1
    if(psd==0):
        for m in range(9):
            j=int(round(12.5*m))+1
            for i in range(lzes,lzee):
                ff,cpp=msc(z[j,0:ntot],z[i,0:ntot],fs=fsam,nperseg=nseg,noverlap=novlp)
                if (i==lzes and j==lzes):
                    Cpp=np.zeros((9,lze,ff.shape[0]))
                Cpp[m,i-1,:]=cpp


    else:
        for i in range(lzes,lzee):
            ff,cpp=psdw(z[i,0:ntot],fs=fsam,nperseg=nseg,noverlap=novlp,scaling='spectrum')
            if (i==lzes):
                Cpp=np.zeros((lze,ff.shape[0]))
            vscl=(np.var(z[i,0:ntot]))
            Cpp[i-1,:]=cpp/vscl

    zz=np.linspace(-0.11,0.11,101)#[i*(2.0/100) for i in range(101)]
    zz/=0.11
    X, Y = np.meshgrid(ff, zz)
    i=7;fmax=1
    m=i

    if (psd==1):
        v = []
        exp0 = -3;exp1=0
        nlev = 50
        for E in range(exp0,exp1):
            v = np.concatenate((v[:-1],np.linspace(10**E,10**(E+1),nlev)))
        Z = np.clip(Cpp, v[0],v[-1])
        CS=plt.contourf(X,Y,Z,norm = LogNorm(vmin=v[0],clip=False),levels=v)#,norm=LogNorm())#,vmin=cmin)
        def fmt(x, pos):
            a, b = '{:.2e}'.format(x).split('e')
            b = int(b)
            return r'${} \times 10^{{{}}}$'.format(a, b)
        v = []
        nlev = 2
        for E in range(exp0,exp1):
            v = np.concatenate((v[:-1],np.linspace(10**E,10**(E+1),nlev)))
        cbar=plt.colorbar(CS, format=ticker.FuncFormatter(fmt), spacing='proportional',ticks=v)
        cbar.ax.set_ylabel('PSD')
    else:
        if (i%2==0):
            countm=i//2+1
            cs="Middle "
            cs+=str(countm)
        else:
            if((i//2)%2==0):
                countp=i//4+1
                cs="Peak "
                cs+=str(countp)
            else:
                countt=i//4+1
                cs="Trough "
                cs+=str(countt)
        v = np.linspace(0.0, 1.0, 11, endpoint=True)
        CS=plt.contourf(X[:,:],Y[:,:],Cpp[m,:,:],v)
        plt.title(cs+' used as base Signal')
        cbar=plt.colorbar(CS)
        cbar.ax.set_ylabel('MSC')

    plt.ylabel(r'$\frac{z}{L_c \lambda_{LE}}$',fontsize=25)
    plt.xlabel(r'$f^*$',fontsize=16)
    plt.xlim(0,fmax)
    plt.show()    
    #plt.savefig(savepath,dpi=300)
    
    
#system('montage -mode concatenate -tile 4x3 /home/rpt1g12/Desktop/circ/'+cpre+
 #      '*.png /home/rpt1g12/Desktop/circ/out'+cpre+'.png')

