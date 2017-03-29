import numpy as np
import matplotlib.pyplot as plt
from lib.myPlots import *
   
#%%
cs='y'
dataset='wallUnits/6blocks/4A15W11AoA20/'+cs+'+.dat'
x,wu=np.loadtxt(dataset,skiprows=1,unpack=True)

#%%
fig=plt.figure()
ax=fig.add_subplot(111)
skip=5;ncol=1
ax.plot(x,wu,markevery=skip,label=cs+'+')
ax.set_xlabel(r'$x/L_c$')
ax.set_ylabel(r'$<'+cs+'^+>$')
handle,labels=ax.get_legend_handles_labels()
legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=ncol)
ax.grid(True)
ax.figure.show()
path='wallUnits/clean/AoA20/'
#plt.savefig(savepath,dpi=300)
#savePlotFile(path=path+cs+'Plus.dat',vary=labels)