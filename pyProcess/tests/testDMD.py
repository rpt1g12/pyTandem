import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import dot, multiply, diag, power
from numpy import pi, exp, sin, cos, cosh, tanh, real, imag
from numpy.linalg import inv, eig, pinv
from scipy.linalg import svd, svdvals
from scipy.integrate import odeint, ode, complex_ode
from warnings import warn

from lib.stats import *
from lib.myPlots import *

plt.close('all')
figs=[];axs=[];hdls=[];lbls=[];lgds=[];nfig=-1;fs=18
#%%
# define time and space domains
x = np.linspace(-10, 10, 60)
t = np.linspace(0, 12*np.pi, 1025)
dt = t[2] - t[1]
Xm,Tm = np.meshgrid(x, t)

#%% create three spatiotemporal patterns
f1 = multiply(20-0.2*power(Xm, 2), exp((2.3j)*Tm))
f2 = multiply(Xm, exp(0.6j*Tm))
f3 = multiply(5*multiply(1/cosh(Xm/2), tanh(Xm/2)), 2*exp((0.1+2.8j)*Tm))

#%% combine signals and make data matrix
D = (f1 + f2 + f3).T
r=4
U,s,V,eV,eVec,Phi,Psi,PhiPOD=ROM(D,r=r)

#%% compute DMD reconstruction
rDMD = dot(Phi, Psi)
np.allclose(D[:,:-1], rDMD) # True

#%%

if r<0:
    r=len(s)
f,a=getFig('DMDmode');figs.append(f),axs.append(a);nfig+=1

plotRange=range(r)
for i in plotRange:
    axs[nfig].plot(Phi.real[:,i],lw=2,label='m{:02d}'.format(i))
#%%
#nmodes=10
f,a=getFig('DMDtime');figs.append(f),axs.append(a);nfig+=1
for i in plotRange:
    axs[nfig].plot(Psi.real[i,:],lw=2,label='m{:02d}'.format(i))

#%%
n=len(t)-2
f,a=getFig('streamline');figs.append(f),axs.append(a);nfig+=1
a.plot(D[:,n],lw=2,color='red',label='original')
a.plot(rDMD[:,n],lw=2,color='blue',linestyle='--',label='DMD')
fit(a)