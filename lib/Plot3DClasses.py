# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 11:24:30 2016

@author: rpt1g12
"""
from lib.myPlot3dOperator import *
import numpy as np
import copy as copy

# Classes
class var():
    def __init__(self,size,vid,name='var',val=[]):
        self.id=vid
        self.name=name
        self.size=size
        if len(val)==0:
            self.val=np.zeros(size)
        else:
            self.val=val.copy()
        
    def rdVar(self,fh,lh,nvar,nb,lblk,Type='grid'):
        lxi=self.size[0];let=self.size[1];lze=self.size[2]
        ltomb=(lxi)*(let)*(lze)
        k8=np.float32; lk8=np.dtype(k8).itemsize
        if Type=='grid':
            tvar=3
            fl_hdr=0*lk8
        if Type=='sol':
            tvar=5
            fl_hdr=4*lk8
        
        lh+=fl_hdr*(nb+1)+sum(lblk[:nb])*lk8*tvar+nvar*lblk[nb]*lk8
        #if Type=='sol': print('var',nb,nvar,lh)
        fh.seek(lh)
        self.val=np.fromfile(fh,dtype=k8,count=ltomb)
        self.val=np.reshape(self.val,(lxi,let,lze),order='F')
    
    def getValues(self):
        return self.val.copy()
    def getSize(self):
        return self.size.copy()
    def getName(self):
        return self.name
        
class blk():
    def __init__(self,blk_id,size=(1,1,1)):
        self.id=blk_id
        self.gsize=size
        self.glen=size[0]*size[1]*size[2]
        self.var={}
        self.data=[]
        
    def getSubset(self,xlim=None,ylim=None,zlim=None):
        ndata=len(self.data)
        size=self.gsize
        if xlim==None:
            xlim=range(0,size[0])
        if ylim==None:
            ylim=range(0,size[1])
        if zlim==None:
            zlim=range(0,size[2])
        lsize=(len(xlim),len(ylim),len(zlim))
        sarr=np.zeros(lsize)
        sblk=blk(self.id,lsize)
        print('subset sizes:'+str(lsize))
        for n in range(ndata):
            arr=self.data[n].getValues()
            ii=-1;
            for i in xlim:
                ii+=1;jj=-1
                for j in ylim:
                    jj+=1;kk=-1
                    for k in zlim:
                        kk+=1
                        sarr[ii,jj,kk]=arr[i,j,k]                     
            name=self.data[n].getName()
            sblk.setData(name,sarr,size=lsize)
            sblk.var[name]=sblk.data[n]
        
        return sblk.clone()
        
    def setData(self,vname=None,val=[],vid=None,size=None):
        if vid==None:
            vid=-1
        if vname==None:
            vname='v{:d}'.format(len(self.data))
        if len(val)==0:
            if size==None:
                size=self.gsize
                val=np.zeros(size)
        else:
            size=val.shape
        
        if vid==-1 or len(self.data)<vid+1:
            self.data.append(var(size,len(self.data),vname,val))
            #print(len(self.data),val)
        else:
            myvar=var(size,vid,vname,val)
            
            self.data[vid]=myvar
    
    def clone (self):
        obj=copy.copy(self)
        return obj 
        
class flow():
        def __init__(self,path,gfile,sfile):
            self.path=path
            self.gfile=gfile
            self.sfile=sfile
            self.nbk=1
            self.blk=[]
            self.lblk=[]
        def rdHdr(self):
            fh=open(self.path+self.gfile,'rb')
            #Read number of blocks
            self.nbk=rdType(fh,'i')
            #Read Grid sizes for each block
            for nb in range(self.nbk):
                self.blk.append(blk(blk_id=nb))
                self.blk[nb].gsize=rdType(fh,'i',3)
                #Compute total size of arrays
                for n in range(3):
                    self.blk[nb].glen*=self.blk[nb].gsize[n]
            #Store header length
            self.lhdr=fh.tell()
            
            #Store total sizes for all blocks
            for nb in range(self.nbk):
                self.lblk.append(self.blk[nb].glen)
            fh.close()

        def rdGrid(self):
            fh=open(self.path+self.gfile,'rb')
            nbk=self.nbk
            names=['x','y','z']
            lblk=self.lblk
            for nb in range(nbk):
                size=self.blk[nb].gsize 
                lh=self.lhdr
                for nvar in range(3):
                    vname=names[nvar]
                    self.blk[nb].data.append(var(size,nvar,name=vname))
                    self.blk[nb].data[-1].rdVar(fh,lh,nvar,nb,lblk)
                    self.blk[nb].var[vname]=self.blk[nb].data[-1]
            fh.close()
        
        def rdSol(self,vnames=None):
            fh=open(self.path+self.sfile,'rb')
            nbk=self.nbk
            if vnames==None:
                names=[]
                for n in range(5):
                    names.append('v{:d}'.format(n))
            else:
                names=vnames
            lblk=self.lblk
            for nb in range(nbk):
                size=self.blk[nb].gsize 
                lh=self.lhdr
                for nvar in range(5):
                    vname=names[nvar]
                    self.blk[nb].data.append(var(size,nvar,name=vname))
                    self.blk[nb].data[-1].rdVar(fh,lh,nvar,nb,lblk,Type='sol')
                    self.blk[nb].var[vname]=self.blk[nb].data[-1]
            fh.close()
        
            
        
