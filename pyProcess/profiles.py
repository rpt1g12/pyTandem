import numpy as np
import matplotlib.pyplot as plt

#%%
path='profiles/SIT/AoA10/trough/profile'
path2='/home/rpt1g12/Dropbox/phd/figures/wleResults/'
for n in range(0,10):
    dataset=path+str(n)+'.dat'
    #prf=[]#np.zeros()
    cp=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=range(1))
    cp = cp[~np.isnan(cp)]
    uin=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=range(1,4))
    rho=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=range(4,5))
    rho = rho[~np.isnan(rho)]
    xin=np.loadtxt(dataset,skiprows=1,unpack=True,usecols=range(9,12))
    for i in range(uin.shape[0]):    
        x1 = xin[i,:][~np.isnan(uin[i,:])]
        if (i==0): x=np.zeros((3,len(x1)));u=x.copy()
        x[i,:]=x1
        x1 = uin[i,:][~np.isnan(uin[i,:])]
        u[i,:]=x1
    var=np.sqrt(u[0,:]*u[0,:]+u[1,:]*u[1,:]+u[2,:]*u[2,:])
    if (n==0): k=0.005/np.absolute(var.max())
    varp=var*k+x[0,0]
    plt.plot(varp[:40],range(40),'b',label=r'$x/L_c=$'+str(x[0,0]))
plt.xlabel(r'$x/L_c$')
plt.ylabel(r'$\eta$')
plt.title(r'$|U|$ profiles')
plt.show()
