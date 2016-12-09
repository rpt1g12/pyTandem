import numpy as np
import matplotlib.pyplot as plt

from libGrid.gridFunctions import *
import libGrid.gridClasses as gcl
from lib.myPlots import *

import importlib
#%%
importlib.reload(gcl)

plt.close('all')
#%%

p0=gcl.point(dyx=[0,1])
p1=gcl.point(x=1.0,y=1)
p2=gcl.point(x=10.0,y=1)

ln=gcl.line(119,[p0,p1,p2],[50],0)
print(ln.io,ln.nseg)

ln.linePlot()