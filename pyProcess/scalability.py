import xlrd
import matplotlib.pyplot as plt
import numpy as np
#%%
book=xlrd.open_workbook('/home/rpt1g12/Dropbox/phd/excel/GridPoints.xls')
sheet=book.sheet_by_index(3)

p=sheet.col_values(0)
p.pop(0)
s=sheet.col_values(5)
s.pop(0)
r=sheet.col_values(6)
r.pop(0)
e=sheet.col_values(7)
e.pop(0)
#%%
fig1=plt.figure()
plt.plot(p,s,'b-s',label='Iridis4',linewidth=2.0,ms=10)
plt.plot(p,r,'r--',label='Linear',linewidth=2.0)
plt.ylabel(r'$S_u$',fontsize=20)
plt.xlabel(r'$n_p$',fontsize=20)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    ncol=2, mode="expand", borderaxespad=0.)       
plt.grid(True)
plt.show()
#%%
fig2=plt.figure()
plt.plot(p,e,'b-s',label='Iridis4',linewidth=2.0,ms=10)
plt.ylabel(r'$P_{eff}$',fontsize=20)
plt.xlabel(r'$n_p$',fontsize=20)
plt.grid(True)
plt.show()