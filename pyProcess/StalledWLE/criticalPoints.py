# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
from lib.myPlot3dOperator import *
import lib.Plot3DClasses as p3d
from scipy import stats
from scipy.interpolate import interp2d
from scipy.integrate import odeint
pi=np.pi
import getpass
user=getpass.getuser()
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')

import importlib
#%%
importlib.reload(p3d)
#%%
A=15 #WLE Amplitude, if SLE A=0
AoA=20 #Angle of Attack
nwave=2 #Number of LE wavelengths
gfile='grid.xyz' #Grid filename

sfile='solTA.qa' #Solution filename
vnames=['r','u','v','w','p'];


#%% Paths set-up
if A>0:
    sfolder='{}WLE'.format(nwave)
else:
    sfolder='{}SLE'.format(nwave)
subpath='average/'
simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/dellHDD/post/{}/{}".format(user,simfolder,subpath)

fl=p3d.flow(path,gfile,sfile)
fl.rdHdr()
fl.rdGrid()
fl.rdSol(vnames,sfile)
#%%
ss=fl.blk[4].getSubset(ylim=[1])
x=ss.var['x'].getValues()[:,0,:]
z=ss.var['z'].getValues()[:,0,:]
u=ss.var['u'].getValues()[:,0,:]
v=ss.var['v'].getValues()[:,0,:]
w=ss.var['w'].getValues()[:,0,:]

#%%
dudxi=ss.derive('u','xi',True)[:,0,:]
dudze=ss.derive('u','ze',True)[:,0,:]
dwdxi=ss.derive('w','xi',True)[:,0,:]
dwdze=ss.derive('w','ze',True)[:,0,:]

#%%
nxi,nze=u.shape
dummy=np.zeros((nxi,1,nze))
dummy[:,0,:]=np.sqrt(u**2+v**2+w**2)
ss.setData('U',val=dummy)
#%%
f,a,im=ss.contourf(varname='U',k=0,plane=1,vmin=0,vmax=0.01,cmap=plt.cm.Reds_r)
a.invert_yaxis()
a.set_aspect('equal')
ss.drawMeshPlane(direction=1,skp=(5,1),showBlock=True,ax=a)

#%%
xi,ze=np.linspace(0,nxi-1,nxi),np.linspace(0,nze-1,nze)
x_xize=interp2d(ze,xi,x)
z_xize=interp2d(ze,xi,z)
u_xize=interp2d(ze,xi,u)
w_xize=interp2d(ze,xi,w)

dudxi_xize=interp2d(ze,xi,dudxi)
dudze_xize=interp2d(ze,xi,dudze)
dwdxi_xize=interp2d(ze,xi,dwdxi)
dwdze_xize=interp2d(ze,xi,dwdze)

uw = lambda p,t: [u_xize(p[1],p[0])[0],-w_xize(p[1],p[0])[0]]
#%%
nlines=0
nbins=0
dt = 100; t0 = 0; nt=100
t1 = nt*dt
t = np.linspace(t0,t1,nt)
xz=np.zeros((2,nt))

kRange=np.linspace(0,nze-1,4*nwave+1)
#%%
minHits=1
nrelease=0
sinks={}
sources={}
#%%Clean all lines
for i in range(nrelease):
    a.lines.pop(-1)
axShow(a);nrelease=0
#%%
iRange=np.linspace(5,(nxi-1)/2,10)
for fctr in [1,-1]:
    nlines=0
    if fctr>0:
        lcolor='blue'
    else:
        lcolor='green'
    for kk in kRange:
        for ii in iRange: 
            nlines+=1
            nrelease+=1
            print('Streamline {:06d} out of {:06d}'.format(nlines,len(kRange)*len(iRange)))
            streamline=odeint(uw,(ii,kk),fctr*t)
            for i in range(nt):
                xz[0,i]=x_xize(streamline[i,1],streamline[i,0])
                xz[1,i]=z_xize(streamline[i,1],streamline[i,0])      
            a.plot(xz[0,:],xz[1,:],lw=2,color=lcolor)
            ibin,kbin = streamline[-1,:]
            xi=np.round(ibin,1);ze=np.round(kbin,1)
            critPt_key='{:3.1f},{:3.1f}'.format(xi,ze)
            ibin=min(max(0,ibin),nxi-1);kbin=min(max(0,kbin),nze-1)
            ibin=int(np.floor(ibin));kbin=int(np.floor(kbin))
            if fctr>0:
                if critPt_key in sinks.keys():
                    sinks[critPt_key][2]+=1
                    
                else:
                    sinks.update({critPt_key:[xi,ze,1.0]})
            else:
                if critPt_key in sources.keys():
                    sources[critPt_key][2]+=1
                    
                else:
                    sources.update({critPt_key:[xi,ze,1.0]})
            
toPop=[]
sinkIndices=[]
for key,item in sinks.items():
    if item[2]>minHits:
        print(item[2])
        xi,ze=item[:2]
        sinkIndices.append([xi,ze])
        x,z=x_xize(ze,xi),z_xize(ze,xi)
        a.scatter(x,z,c='yellow',s=150)
    else:
        toPop.append(key)            
for key in toPop:
    sinks.pop(key)                

toPop=[]
sourceIndices=[]
for key,item in sources.items():
    if item[2]>minHits:
        print(item[2])
        xi,ze=item[:2]
        sourceIndices.append([xi,ze])
        x,z=x_xize(ze,xi),z_xize(ze,xi)
        a.scatter(x,z,c='magenta',s=150)
    else:
        toPop.append(key)            
for key in toPop:
    sources.pop(key)                 
        
axShow(a) 

#%%
nRays=8;radii=[[1,1]]
dAlpha=2*pi/nRays
phase=0.0
for fctr in [1,-1]:
    nlines=0
    phase=0.0
    if fctr>0:
        lcolor='blue'
        critIndices=sourceIndices
    else:
        lcolor='green'
        critIndices=sinkIndices
    for index in critIndices:
        xi0,ze0=index[:]
        for radius in radii:
            for nn in range(nRays):
                nlines+=1
                nrelease+=1 
                dxi=radius[0]*np.cos(nn*dAlpha+phase)
                dze=radius[1]*np.sin(nn*dAlpha+phase)
                print('Streamline {:06d} out of {:06d}'.format(nlines,(nRays)*len(radii)*len(critIndices)))
                ii=xi0+dxi;kk=ze0+dze
                streamline=odeint(uw,(ii,kk),fctr*t)
                for i in range(nt):
                    xz[0,i]=x_xize(streamline[i,1],streamline[i,0])
                    xz[1,i]=z_xize(streamline[i,1],streamline[i,0])      
                a.plot(xz[0,:],xz[1,:],lw=2,color=lcolor)
                ibin,kbin = streamline[-1,:]
                xi=np.round(ibin,1);ze=np.round(kbin,1)
                critPt_key='{:3.1f},{:3.1f}'.format(xi,ze)
                ibin=min(max(0,ibin),nxi-1);kbin=min(max(0,kbin),nze-1)
                ibin=int(np.floor(ibin));kbin=int(np.floor(kbin))

                if fctr>0:
                    if critPt_key in sinks.keys():
                        sinks[critPt_key][2]+=1
                        
                    else:
                        sinks.update({critPt_key:[xi,ze,1.0]})
                else:
                    if critPt_key in sources.keys():
                        sources[critPt_key][2]+=1
                        
                    else:
                        sources.update({critPt_key:[xi,ze,1.0]})
            phase+=dAlpha/2

toPop=[]
sourceIndices=[]
for key,item in sources.items():
    if item[2]>minHits:
        print(item[2])
        xi,ze=item[:2]
        sourceIndices.append([xi,ze])
        x,z=x_xize(ze,xi),z_xize(ze,xi)
        a.scatter(x,z,c='magenta',s=150)
    else:
        toPop.append(key)            
for key in toPop:
    sources.pop(key)                
                
toPop=[]
sinkIndices=[]
for key,item in sinks.items():
    if item[2]>minHits:
        print(item[2])
        xi,ze=item[:2]
        sinkIndices.append([xi,ze])
        x,z=x_xize(ze,xi),z_xize(ze,xi)
        a.scatter(x,z,c='yellow',s=150)
    else:
        toPop.append(key)            
for key in toPop:
    sinks.pop(key) 
        
axShow(a) 

#%%
fPQ,aPQ=getFig('PQ')
sourcePQ=[]
sinkPQ=[]
for indices in sourceIndices:
    xi,ze=indices[:]
    dudx=dudxi_xize(ze,xi)[0]
    dudz=-dudze_xize(ze,xi)[0]
    dwdx=dwdxi_xize(ze,xi)[0]
    dwdz=-dwdze_xize(ze,xi)[0]
    P=-(dudx+dwdz)
    Q=dudx*dwdz-dudz*dwdx
    aPQ.scatter(P,Q,c='magenta',s=100)
    sourcePQ.append([P,Q])

for indices in sinkIndices:
    xi,ze=indices[:]
    dudx=dudxi_xize(ze,xi)[0]
    dudz=-dudze_xize(ze,xi)[0]
    dwdx=dwdxi_xize(ze,xi)[0]
    dwdz=-dwdze_xize(ze,xi)[0]
    P=-(dudx+dwdz)
    Q=dudx*dwdz-dudz*dwdx
    aPQ.scatter(P,Q,c='yellow',s=100)
    sinkPQ.append([P,Q])  
    
sourcePQ=np.asarray(sourcePQ)
sinkPQ=np.asarray(sinkPQ)
minP=min(np.concatenate((sinkPQ[:,0],sourcePQ[:,0])))
maxP=max(np.concatenate((sinkPQ[:,0],sourcePQ[:,0])))
minQ=min(np.concatenate((sinkPQ[:,1],sourcePQ[:,1])))
maxQ=max(np.concatenate((sinkPQ[:,1],sourcePQ[:,1])))
widthP=max(np.abs(minP),np.abs(maxP))
widthQ=max(np.abs(minQ),np.abs(maxQ))
widthMax=max(widthP,widthQ)
p=np.linspace(-widthMax*1.1,1.1*widthMax,151)
q=(p**2)/4
aPQ.plot(p,q,lw=2,color='black')
aPQ.axvline(x=0,lw=2,color='black')
aPQ.axhline(y=0,lw=2,color='black')
aPQ.set_xlim(-widthP*1.1,widthP*1.1)
aPQ.set_ylim(-widthQ*1.1,widthQ*1.1)
#%%
print('\nStable Foci (x,z):\n')
for i in range(len(sinkIndices)):
    xi,ze=sinkIndices[i][:]
    P,Q=sinkPQ[i,:]
    x,z=x_xize(ze,xi)[0],z_xize(ze,xi)[0]
    if Q>(P**2)/4:
        print('({:3.5f},{:3.5f})'.format(x,z))

print('\nUnstable Foci (x,z):\n')
for i in range(len(sourceIndices)):
    xi,ze=sourceIndices[i][:]
    P,Q=sinkPQ[i,:]
    x,z=x_xize(ze,xi)[0],z_xize(ze,xi)[0]
    if Q>(P**2)/4:
        print('({:3.5f},{:3.5f})'.format(x,z))
        