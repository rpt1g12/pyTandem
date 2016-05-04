import numpy as np
import matplotlib.pyplot as plt
from lib.myPlots import *
import struct
import os

#%%

def rdType(file,t,verbose=False):
    f=file
    if(t=='i'): l=4
    elif(t=='f'): l=4
    elif(t=='c'): l=1
    elif(t=='d'): l=8
    v=struct.unpack(t,f.read(l))[0]
    if (verbose):
        print(v)
    return v
    
def rdString(file,verbose=False):
    f=file    
    eflag=1
    s=''
    lh=f.tell()
    #s=f.read(1)
    while (eflag!=0):
        f.seek(lh)
        c=rdType(f,'c')
        c=chr(c[0]);lh+=4;
        f.seek(lh)
        eflag=rdType(f,'i')
        s+=c
    if(verbose):
        print(s)
    return s

def rdVar(f,lh,nvar,lxi,let,lze):
    ltomb=(lxi)*(let)*(lze)
    k8=np.float32; lk8=np.dtype(k8).itemsize
    lh+=nvar*ltomb*lk8
    fh=f
    fh.seek(lh)
    x=np.fromfile(fh,dtype=k8,count=ltomb)
    x=np.reshape(x,(lxi,let,lze),'F')
    return x

def wrP3dGrid(path,cblock,lxi,let,lze,xyz):
    """Write to Plot3D format"""
    
    fh=open(path+cblock+'.xyz','wb')
    hdr=np.int32(np.array([1,lxi,let,lze]))
    fh.write(hdr)
    for i in range(3):
        fh.write(np.float32(np.transpose(xyz[:,:,:,i]).copy(order='C')))
    fh.close()
    print(cblock+' written!')
    pass

def wrP3dFunction(path,name,lxi,let,lze,f,r):
    """Write to Plot3D format"""
    
    fh=open(path+name+'.f','wb')
    nf=len(r)
    hdr=np.int32(np.array([1,lxi,let,lze,nf]))
    fh.write(hdr)
    for i in r:
        fh.write(np.float32(np.transpose(f[:,:,:,i]).copy(order='C')))
    fh.close()
    print(name+' written!')
    pass

def wrP3dS(path,cblock,size,flowc,f,r):
    """Write to Plot3D format"""
    
    ctime='{:08.4f}'.format(flowc[3])
    fh=open(path+'/solT'+ctime+'b'+cblock+'.q','wb')
    hdr=np.int32(np.array([1,size[0],size[1],size[2]]))
    fh.write(hdr)
    fhdr=np.float32(np.array(flowc))
    fh.write(fhdr)
    for i in r:
        fh.write(np.float32(np.transpose(f[:,:,:,i]).copy(order='C')))
    fh.close()
    print(ctime+' written!')
    pass

#%%
for m in [11]:
    
    cblock='{:03d}'.format(m)
    path='tecplotData/'
    directory=path+'b'+cblock+'/'   
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    cname='output'+cblock+'.plt'
    f=open(path+cname,'rb')
    vbs=False
#%%
    f.seek(0)
    tecVersion=f.read(8);print(tecVersion);
    hSec=rdType(f,'i',vbs)
    fType=rdType(f,'i',vbs)
    title=rdString(f,vbs)
    nvar=rdType(f,'i',vbs)
    cvar=[]
    for i in range(nvar):
        cvar.append(rdString(f,vbs))
    
    znMark0=rdType(f,'f',vbs)
    znName=rdString(f,vbs)
    parentZone=rdType(f,'i',vbs)
    strandID=rdType(f,'i',vbs)
    time=rdType(f,'d',vbs)
    notUsed=rdType(f,'f',vbs)
    znType=rdType(f,'i',vbs)
    varLoc=rdType(f,'i',vbs)
    neighbours=rdType(f,'i',vbs)
    miscNeighb=rdType(f,'i',vbs)
    imax=rdType(f,'i',vbs)
    jmax=rdType(f,'i',vbs)
    kmax=rdType(f,'i',vbs)
    nAux=rdType(f,'i',vbs)
    hMark=rdType(f,'i',vbs)
    znMark1=rdType(f,'f',vbs)
    
    varType=[]
    for i in range(nvar):
        varType.append(rdType(f,'i',vbs))
    
    nPasVar=rdType(f,'i',vbs)
    nShrVar=rdType(f,'i',vbs)
    zBased=rdType(f,'i',vbs)
    minVar=[];maxVar=[]
    for i in range(nvar):
        minVar.append(rdType(f,'d',vbs))
        maxVar.append(rdType(f,'d',vbs))
    lhdr=f.tell()
    
    var=np.zeros((imax,jmax,kmax,nvar))
    for i in range(nvar):
        var[:,:,:,i]=rdVar(f,lhdr,i,imax,jmax,kmax)
    
    #%%
    wrP3dGrid(directory,'gridb'+cblock,imax,jmax,kmax,var)
    wrP3dFunction(directory,'tfunctionb'+cblock,imax,jmax,kmax,var,range(3,12))
    #%%
    nsol0=cvar.index('r000')
    nsol=int((len(cvar)-nsol0)/5);print(str(nsol)+' Snapshots found')
    sizes=[imax,jmax,kmax]
    
    for i in range(nsol):
        flowc=[0.4,5.0,120000.0,i*1.0]
        r0=nsol0+i*5;
        rn=r0+5
        r=range(r0,rn)
        wrP3dS(directory,cblock,sizes,flowc,var,r)
    
    f.close()