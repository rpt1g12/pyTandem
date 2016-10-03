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
        #self.flowc=np.zeros(4)
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
        
        lh+=fl_hdr*(nb)+sum(lblk[:nb])*lk8*tvar
        if Type=='sol':
            fh.seek(lh)
            flowc=np.fromfile(fh,dtype=k8,count=4)
            self.flowc=flowc.copy()
            
            
        lh+=fl_hdr+nvar*lblk[nb]*lk8
        #if Type=='sol': print('var',nb,nvar,lh)
        fh.seek(lh)
        self.val=np.fromfile(fh,dtype=k8,count=ltomb)
        self.val=np.reshape(self.val,(lxi,let,lze),order='F')
    
    def clone(self):
        obj=copy.copy(self)
        return obj
        
    def wrVar(self,fh,lh,nvar,nb,lblk,Type='grid'):
        k8=np.float32; lk8=np.dtype(k8).itemsize
        if Type=='grid':
            tvar=3
            fl_hdr=0*lk8
        if Type=='sol':
            tvar=5
            fl_hdr=4*lk8
        
        lh+=fl_hdr*(nb)+sum(lblk[:nb])*lk8*tvar
        if Type=='sol':
            fh.seek(lh)
            flowc=self.flowc.copy()
            fhdr=np.float32(np.array(flowc))
            fh.write(fhdr)
        
        lh+=fl_hdr+nvar*lblk[nb]*lk8
        #if Type=='sol': print('var',nb,nvar,lh)
        fh.seek(lh)
        v=(np.float32(np.transpose(self.getValues()).copy(order='C')))
        fh.write(v)
    
    def derVar(self,direction):
        direction
        dvar=np.zeros(self.size)
        n=self.size[direction]
        dx=np.zeros(n)
        a=np.zeros(n-2);b=np.zeros(n-1);c=np.zeros(n)
        d=np.zeros(n-1);e=np.zeros(n-2);
        a[-1]=1;e[0]=-1
        b[0:-1]=-1;b[-1]=-4;d[1:]=1;d[0]=4;
        c[0]=-3;c[-1]=3
        m=np.mat(np.diag(a,-2)+np.diag(b,-1)+np.diag(c,0)+np.diag(d,1)+np.diag(e,2))
        if (direction==0):
            for k in range(self.size[2]):
                for j in range(self.size[1]):
                    dx=m*np.mat(self.val[:,j,k]).T
                    dvar[:,j,k]=np.asarray(dx.T)[0]
        if (direction==1):
            for i in range(self.size[0]):
                for k in range(self.size[2]):
                    dx=m*np.mat(self.val[i,:,k]).T
                    dvar[i,:,k]=np.asarray(dx.T)[0]
        if (direction==2):
            for j in range(self.size[1]):
                for i in range(self.size[0]):
                    dx=m*np.mat(self.val[i,j,:]).T
                    dvar[i,j,:]=np.asarray(dx.T)[0]
        return dvar.copy()
            
            

    
    def avgDir(self,direction=2):
        size=self.size        
        if direction==0:
            lsize=size[1:]
            aarr=np.zeros(lsize)
            arr=self.getValues()
            for i in range(size[0]):
                aarr+=arr[i,:,:]
        if direction==1:
            lsize=(size[0],size[2])
            aarr=np.zeros(lsize)
            arr=self.getValues()
            for j in range(size[1]):
                aarr+=arr[:,j,:]
        if direction==2:
            lsize=size[0:2]
            aarr=np.zeros(lsize)
            arr=self.getValues()
            for k in range(size[2]):
                aarr+=arr[:,:,k]
        aarr/=size[direction]
        return aarr.copy()
    
    def getValues(self):
        return self.val.copy()
    def getSize(self):
        return self.size.copy()
    def getName(self):
        return self.name
        
class blk():
    def __init__(self,blk_id,size=(1,1,1)):
        self.id=blk_id
        self.size=size
        self.glen=size[0]*size[1]*size[2]
        self.var={}
        self.data=[]
        self.metFlag=False
        
    def getMetrics(self):
        flag=(not self.metFlag)
        if flag:
            self.mets=[]
            self.imets=[]
            x=self.var['x'].clone()
            y=self.var['y'].clone()
            z=self.var['z'].clone()
            self.imets.append(var(self.size,vid=0,name='dxdxi',val=x.derVar(0)))
            self.imets.append(var(self.size,vid=1,name='dxdet',val=x.derVar(1)))
            self.imets.append(var(self.size,vid=2,name='dxdze',val=x.derVar(2)))  
            self.imets.append(var(self.size,vid=3,name='dydxi',val=y.derVar(0)))
            self.imets.append(var(self.size,vid=4,name='dydet',val=y.derVar(1)))
            self.imets.append(var(self.size,vid=5,name='dydze',val=y.derVar(2)))
            self.imets.append(var(self.size,vid=6,name='dzdxi',val=z.derVar(0)))
            self.imets.append(var(self.size,vid=7,name='dzdet',val=z.derVar(1)))
            self.imets.append(var(self.size,vid=8,name='dzdze',val=z.derVar(2)))
            
            rdum=self.imets[4].getValues()*self.imets[8].getValues()
            rdum-=self.imets[5].getValues()*self.imets[7].getValues()
            self.mets.append(var(self.size,vid=0,name='dxidx',val=rdum))
            rdum*=self.imets[0].getValues()
            self.J=var(self.size,vid=0,name='J',val=rdum)
            
            
            rdum=self.imets[2].getValues()*self.imets[7].getValues()
            rdum-=self.imets[1].getValues()*self.imets[8].getValues()
            self.J.val[:,:,:]+=rdum.copy()*self.imets[3].getValues()           
            self.mets.append(var(self.size,vid=1,name='dxidy',val=rdum))
            
            rdum=self.imets[1].getValues()*self.imets[5].getValues()
            rdum-=self.imets[2].getValues()*self.imets[4].getValues()
            self.J.val[:,:,:]+=rdum.copy()*self.imets[6].getValues()            
            self.mets.append(var(self.size,vid=2,name='dxidz',val=rdum))
            
            rdum=self.imets[5].getValues()*self.imets[6].getValues()
            rdum-=self.imets[3].getValues()*self.imets[8].getValues()
            self.J.val[:,:,:]+=rdum.copy()*self.imets[1].getValues()          
            self.mets.append(var(self.size,vid=3,name='detdx',val=rdum))
            
            rdum=self.imets[0].getValues()*self.imets[8].getValues()
            rdum-=self.imets[2].getValues()*self.imets[6].getValues()
            self.J.val[:,:,:]+=rdum.copy()*self.imets[4].getValues()            
            self.mets.append(var(self.size,vid=4,name='detdy',val=rdum))
            
            rdum=self.imets[2].getValues()*self.imets[3].getValues()
            rdum-=self.imets[0].getValues()*self.imets[5].getValues()
            self.J.val[:,:,:]+=rdum.copy()*self.imets[7].getValues()            
            self.mets.append(var(self.size,vid=5,name='detdz',val=rdum))
            
            rdum=self.imets[3].getValues()*self.imets[7].getValues()
            rdum-=self.imets[4].getValues()*self.imets[6].getValues() 
            self.J.val[:,:,:]+=rdum.copy()*self.imets[2].getValues()
            self.mets.append(var(self.size,vid=6,name='dzedx',val=rdum))
            
            rdum=self.imets[1].getValues()*self.imets[6].getValues()
            rdum-=self.imets[0].getValues()*self.imets[7].getValues()
            self.J.val[:,:,:]+=rdum.copy()*self.imets[5].getValues()           
            self.mets.append(var(self.size,vid=7,name='dzedy',val=rdum))
            
            rdum=self.imets[0].getValues()*self.imets[4].getValues()
            rdum-=self.imets[1].getValues()*self.imets[3].getValues()
            self.J.val[:,:,:]+=rdum.copy()*self.imets[8].getValues()            
            self.mets.append(var(self.size,vid=8,name='dzedz',val=rdum))
            
            self.J.val[:,:,:]=3/(self.J.getValues())
            rdum=1/self.J.getValues()
            self.V=var(self.size,vid=0,name='V',val=rdum)
            
            for n in range(9):
                self.mets[n].val[:,:,:]*=self.J.val[:,:,:]
            
            self.metFlag=True
            
    def derive(self,varName1,varName2):
        if varName2 =='x'
            
        else:
            print('Needs to be derived with respect to: x, y or z')
            
            
    def getSubset(self,xlim=None,ylim=None,zlim=None):
        ndata=len(self.data)
        size=self.size
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
                size=self.size
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
                self.blk[nb].size=rdType(fh,'i',3)
                #Compute total size of arrays
                for n in range(3):
                    self.blk[nb].glen*=self.blk[nb].size[n]
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
                size=self.blk[nb].size 
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
                size=self.blk[nb].size 
                lh=self.lhdr
                for nvar in range(5):
                    vname=names[nvar]
                    self.blk[nb].data.append(var(size,nvar,name=vname))
                    self.blk[nb].data[-1].rdVar(fh,lh,nvar,nb,lblk,Type='sol')
                    self.blk[nb].var[vname]=self.blk[nb].data[-1]
                    
            fh.close()

        def wrSol(self,sfile=None):
            if sfile==None:
                sfile=self.sfile
            fh=open(self.path+sfile,'wb')
            nbk=self.nbk
            lblk=self.lblk
            ghdr=[nbk]
            for nb in range(nbk): 
                lh=self.lhdr
                ghdr.append(self.blk[nb].size[0])
                ghdr.append(self.blk[nb].size[1])
                ghdr.append(self.blk[nb].size[2])
                for nvar in range(5):
                    self.blk[nb].data[nvar+3].wrVar(fh,lh,nvar,nb,lblk,Type='sol')
            fh.seek(0)
            fhdr=np.int32(np.array(ghdr))
            fh.write(fhdr)
            fh.close()
        
            
        
