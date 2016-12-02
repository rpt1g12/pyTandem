try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()
import os
import numpy as np
import getpass
user=getpass.getuser()

p0=FindSource('from')
p1=FindSource('to')
cl=[];cd=[];xm=[];xl=[]
frm=-0.488
print(frm)
to=-0.485*np.cos(np.deg2rad(20))
xs=np.linspace(frm,to,7)
lx0=len(xs)-1
for i in range(lx0):
        x0=frm
        x1=to
        xl.append(xs[i+1]-xs[i])
        xm.append(0.5*(xs[i]+xs[i+1]))
        p0.ClipType.Origin=[x0,0,0]
        p1.ClipType.Origin=[x1,0,0]
        cd.append((((servermanager.Fetch(FindSource('IntegrateVariables1'))).GetPointData()).GetArray('pn')).GetValue(0)/lx0)
        cl.append((((servermanager.Fetch(FindSource('IntegrateVariables1'))).GetPointData()).GetArray('pn')).GetValue(1)/lx0)
        print(xm[-1],xl[-1],cl[-1],cd[-1])

frm=to
to=0.5*np.cos(np.deg2rad(20))
x=np.linspace(frm,to,198)
lx=len(x)-1
for i in range(lx):
        x0=float(x[i]);x1=float(x[i+1])
        xl.append(x1-x0)
        xm.append(0.5*(x1+x0))
        p0.ClipType.Origin=[x0,0,0]
        p1.ClipType.Origin=[x1,0,0]
        cd.append((((servermanager.Fetch(FindSource('IntegrateVariables1'))).GetPointData()).GetArray('pn')).GetValue(0))
        cl.append((((servermanager.Fetch(FindSource('IntegrateVariables1'))).GetPointData()).GetArray('pn')).GetValue(1))
        print(xm[-1],xl[-1],cl[-1],cd[-1])

print(to)
data=[xm,xl,cl,cd]
n=len(cl);m=len(data)
s='xm \t xl \t cl \t cd \t'
s+='\n'
for i in range(n):
        for j in range(m):
                s+='{:12.5e} \t'.format(data[j][i])
        s+='\n'

path='/home/'+user+'/anaconda3/pyTandem/pgfPlots/8SLE_pint.dat'
print (path)

fh=open(path,'w')
fh.write(s)
fh.close()

del xm,xl,cl,cd,data,s
