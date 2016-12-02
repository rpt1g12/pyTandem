# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 11:04:05 2016

@author: rpt1g12
"""

import numpy as np
import matplotlib.pyplot as plt
from lib.myPlots import *
import struct
import os

def rdType(file,Type,count=1,verbose=False):
    f=file
    if(Type=='i'): l=4
    elif(Type=='f'): l=4
    elif(Type=='c'): l=1
    elif(Type=='d'): l=8
    if count==1:
        v=struct.unpack(Type,f.read(l))[0]
    else:
        v=[]
        for i in range(count):
            v.append(struct.unpack(Type,f.read(l))[0])
#    else:
#        l*=count
#        v=struct.unpack(Type,f.read(l))
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

def rdVar(fh,lh,nvar,size):
    lxi=size[0];let=size[1];lze=size[2]
    ltomb=(lxi)*(let)*(lze)
    k8=np.float32; lk8=np.dtype(k8).itemsize
    lh+=nvar*ltomb*lk8
    fh.seek(lh)
    x=np.fromfile(fh,dtype=k8,count=ltomb)
    x=np.reshape(x,(lxi,let,lze),'F')
    return x

def wrP3dGrid(path,cblock,size,xyz):
    """Write to Plot3D format"""
    lxi=size[0];let=size[1];lze=size[2]
    fh=open(path+cblock+'.xyz','wb')
    hdr=np.int32(np.array([1,lxi,let,lze]))
    fh.write(hdr)
    for i in range(3):
        fh.write(np.float32(np.transpose(xyz[:,:,:,i]).copy(order='C')))
    fh.close()
    print(cblock+' written!')
    pass

def wrP3dFunction(path,name,size,f,r):
    """Write to Plot3D format"""
    lxi=size[0];let=size[1];lze=size[2]
    fh=open(path+name+'.f','wb')
    nf=len(r)
    hdr=np.int32(np.array([1,lxi,let,lze,nf]))
    fh.write(hdr)
    for i in r:
        fh.write(np.float32(np.transpose(f[:,:,:,i]).copy(order='C')))
    fh.close()
    print(name+' written!')
    pass

def wrP3dS(path,size,flowc,f,r,cname='copy'):
    """Write to Plot3D format"""
    
    ctime='{:08.4f}'.format(flowc[3])
    fh=open(path+'/solT'+ctime+'_'+cname+'.q','wb')
    hdr=np.int32(np.array([1,size[0],size[1],size[2]]))
    fh.write(hdr)
    fhdr=np.float32(np.array(flowc))
    fh.write(fhdr)
    for i in r:
        fh.write(np.float32(np.transpose(f[:,:,:,i]).copy(order='C')))
    fh.close()
    print(ctime+' written!')
    pass

def wrP3dHdr(path,nbk,sizes):
    fh=open(path,'wb')
    s=np.int32(np.array([nbk]));fh.write(s)
    for m in range(nbk):
        hdr=np.int32(sizes[:,m]);fh.write(hdr)
    fh.close()

def wrP3dBlk(path,var,r,flowc=None):
    fh=open(path,'ab')
    if (flowc!=None):
            fhdr=np.float32(np.array(flowc))
            fh.write(fhdr)   
    for i in r:
        fh.write(np.float32(np.transpose(var[:,:,:,i]).copy(order='C')))
    fh.close()
