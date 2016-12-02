import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
#%%
x=[83200,176000,735000,1490000]
y=[0.56,0.78,1.15,1.21]
#%%
n=len(x);nf=100
tmp0=0
tmp1=0
tmp2=0
tmp3=0
tmp4=0
for i in range(n):
    tmp0+=y[i]*np.log(x[i])
    tmp1+=y[i]
    tmp2+=np.log(x[i])
    tmp3+=(np.log(x[i]))**2
    
b=(n*tmp0-tmp1*tmp2)/(n*tmp3-tmp2**2)
a=(tmp1-b*tmp2)/n

f=np.zeros(nf);xf=np.linspace(x[0],x[-1],nf)
for i in range(nf):
    f[i]=a+b*np.log(xf[i])
 
plt.close('all')   
fig=plt.figure()
fig.canvas.set_window_title('Fit')
ax=fig.add_subplot(111)
#ax.plot(x,y,'bo')
ax.plot(xf,f,'k--')

savePlotFile(path='clFit.dat')
print('(a,b)= ',a,b)