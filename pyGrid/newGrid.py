from lib.myGrid import *
from lib.gtools import *

plt.close('all')
fig,ax=getFig('lines')
#%%
p0=corner(0.0,0.0);p1=corner(2.0,1.0)
l=line(101,p0,p1,0)
#%%
l.c1.y=-2
l.c0.dydx=0
l.c1.dydx=0;
l.c1.dxdxi=0.01
l.c0.d2xdxi2=0.0002

l.linePlot(fig)