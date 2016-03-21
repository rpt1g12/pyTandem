# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt

re=[]
re=(83200,176000,735000,1490000)
cl=[]
cl=(0.56,0.78,1.15,1.22)

plt.loglog(re,cl)
plt.show()