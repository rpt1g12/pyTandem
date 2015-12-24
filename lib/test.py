'''from plots2d import *

p1=plotwall()

p1.getData(0,2)

for i in range(0,51):
	p1.addPlot(i)

p1.showPlot()'''

import numpy as np
import matplotlib.pyplot as plt
from stats import *
y,x = np.loadtxt('signal.dat',skiprows=1,unpack=True,usecols=(0,1))
r=acorr(y)
s=np.array([t-25.0 for t in x])
plt.plot(s,r)
plt.show()

