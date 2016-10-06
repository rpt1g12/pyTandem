# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 11:24:30 2016

@author: rpt1g12
"""
from lib.myPlot3dOperator import *
import numpy as np
import copy as copy
from scipy.interpolate import griddata  

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
    
    def getValues(self,oned=False):
        if oned:
           values=np.reshape(self.val.copy(),self.val.size) 
        else:
           values=self.val.copy()
        return values
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
        if varName1 in self.var:
            var0=self.var[varName1]
            dvdxi=var0.derVar(0)
            dvdet=var0.derVar(1)
            dvdze=var0.derVar(2)         
            if varName2 =='x':           
                dv=(dvdxi*self.mets[0].getValues()+
                   dvdet*self.mets[3].getValues()+
                   dvdze*self.mets[6].getValues())
            elif varName2 =='y':           
                dv=(dvdxi*self.mets[1].getValues()+
                   dvdet*self.mets[4].getValues()+
                   dvdze*self.mets[7].getValues())
            elif varName2 =='z':           
                dv=(dvdxi*self.mets[2].getValues()+
                   dvdet*self.mets[5].getValues()+
                   dvdze*self.mets[8].getValues())
            else:
                print('Needs to be derived with respect to: x, y or z')
            self.setData(vname='d{}d{}'.format(varName1,varName2),val=dv)
        else:
            print('{} is not a valid variable!'.format(varName1))

    def interpolate(self,vname,xi,yi,zi,method='nearest'):
        """Interpolate data form variable (vname) into points (ipts)
           Available methods: 'nearest' (default) or 'linear'
           Returns a block object"""
        x=self.var['x'].getValues(oned=True)
        y=self.var['y'].getValues(oned=True)
        z=self.var['z'].getValues(oned=True)
        values=self.var[vname].getValues(oned=True)
        points=np.array(np.transpose([x,y,z]))
        ipoints=np.array(np.transpose([np.reshape(xi,xi.size),np.reshape(yi,yi.size),np.reshape(zi,zi.size)]))
        isize=list(xi.shape)
        ival=griddata(points,values,ipoints,method=method)
        iblock=blk(self.id,size=isize)
        iblock.setData(vname='x',val=xi)
        iblock.setData(vname='y',val=yi)
        iblock.setData(vname='z',val=zi)
        iblock.setData(vname,val=np.reshape(ival,isize))
        return iblock.clone()

    def interpolate2dk(self,vname,xi,yi,ki,method='cubic'):
        """Interpolate data form variable (vname) into points (ipts)
           at ki span-plane.
           Available methods: 'cubic' (default) ,'linear', 'nearest'
           Returns a block object"""
        size=self.size[0]*self.size[1]
        x=self.var['x'].val[:,:,ki];x=np.reshape(x,size)
        y=self.var['y'].val[:,:,ki];y=np.reshape(y,size)
        z=self.var['z'].val[:,:,ki];z=np.reshape(z,size)
        points=np.array(np.transpose([x,y]))
        values=self.var[vname].val[:,:,ki];values=np.reshape(values,size)
        ipoints=np.array(np.transpose([np.reshape(xi,xi.size),np.reshape(yi,yi.size)]))
        isize=[xi.shape[0],xi.shape[1],1]
        print('Interpolating...')
        ival=griddata(points,values,ipoints,method=method)
        iblock=blk(self.id,size=isize)
        iblock.setData(vname='x',val=xi)
        iblock.setData(vname='y',val=yi)
        zi=xi.copy();zi[:,:]=z[0]
        iblock.setData(vname='z',val=zi)
        iblock.setData(vname,val=np.reshape(ival,isize))
        return iblock.clone()

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
            sblk.setData(name,sarr)
            sblk.var[name]=sblk.data[n]
        
        return sblk.clone()
        
    def setData(self,vname=None,val=[],vid=None):
        if vname==None:
            vname='v{:d}'.format(len(self.data))
        if len(val)==0:
            size=self.size
            val=np.zeros(size)
        else:
            size=val.shape
           
        if (vname in self.var):
            print('Changing {} values'.format(vname))
            self.var[vname].val[:,:,:]=val.copy()
        else:
            print('Adding {} to variable list'.format(vname))
            self.data.append(var(size,len(self.data),vname,val))
            self.var[vname]=self.data[-1]
            
      
    def clone (self):
        obj=copy.copy(self)
        return obj 
        
class flow():
        """This class provides a flow object. 
        A flow object contains blocks from a multiblock structured grid.
        Data is read/writen from/to files stored in Plot3D format
        Parameters:
         * path: (char)file path where the grid and solution files are stored
         * gfile: (char) filename for the grid file
         * sfile: (char) filename for the solution file
         * nbk: (int) number of blocks
         * blk: list(block) containing the block objects
         * lblk: list(int) containing the number of points in each block
         * lhdr: (int) length of the header in files counted in bytes
         """
        def __init__(self,path,gfile,sfile):
            #Path where grid and solution files are stored
            self.path=path
            #Grid file name
            self.gfile=gfile
            #Solution file name
            self.sfile=sfile
            #Number of blocks
            self.nbk=1
            #Blocks list
            self.blk=[]
            #Total number of points per block list
            self.lblk=[]
        def rdHdr(self):
            """Reads header from the grid file and sets up basic class parameters such as:
             * nbk
             * lhdr

             It also initialises the blk list containing block objects

             The header is in the following format:

                int4 #Number of blocks (nbk)
                nbk*3*int4 #Number of points in each direction for each block block.size(0:2)
             """
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
            #Open file for reading in binary
            fh=open(self.path+self.gfile,'rb')
            #Number of blocks
            nbk=self.nbk
            #Variable names
            names=['x','y','z']
            #Number of points per block
            lblk=self.lblk
            #Loop trough the blocks
            for nb in range(nbk):
                #Number of points in each direction
                size=self.blk[nb].size 
                #Find out length of header
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
                    
            self.vnames=names
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

        def shiftK(self,k):
            """Shifts the flow object in the ze direction so the kth point is set to be 0 in the returning flow object
               Returns a flow object"""
            # Clone the current flow object
            sflow=self.clone()
            for nb in range(sflow.nbk):
                for n in sflow.vnames:
                    vr2=sflow.blk[nb].var[n].getValues()
                    lze=sflow.blk[nb].size[2]
                    for sk in range(lze):
                        kk=sk+(k-lze)
                        sflow.blk[nb].var[n].val[:,:,sk]=vr2[:,:,kk].copy()
                        
            return sflow
        
        def mergeBlocks(self,bkx,bky):
            """Merges all blocks into a single one
               Returns a flow object"""
            lxi=0;let=0;lze=self.blk[0].size[2]
            i1=[];i0=[];j0=[];j1=[]
            for i in range(bkx):
                i0.append(lxi)
                blxi=self.blk[i].size[0]
                lxi+=(blxi-1)
                i1.append(lxi)
            for j in range(0,self.nbk,bkx):
                j0.append(let)
                blet=self.blk[j].size[1]
                let+=(blet-1)
                j1.append(let)
            lxi+=1;let+=1
            i1[-1]+=1;j1[-1]+=1
            print(i0,i1,lxi)
            print(j0,j1,let)
            v=np.zeros((lxi,let,lze))
            size=[lxi,let,lze]
            mfl=flow(self.path,'merged_{}'.format(self.gfile),'merged_{}'.format(self.sfile))
            mfl.nbk=1
            mfl.lhdr=16
            mfl.lblk.append(lxi*let*lze)
            mfl.vnames=self.vnames.copy()
            mfl.blk.append(blk(blk_id=0))
            krange=range(0,lze)
            mfl.vnames.append('x')
            mfl.vnames.append('y')
            mfl.vnames.append('z')
            for j in range(bky):
                jrange=range(j0[j],j1[j])
                for i in range(bkx):
                    nb=j*(bkx-1)+i
                    nvar=0
                    irange=range(i0[i],i1[i])
                    for vname in mfl.vnames:
                        vv=self.blk[nb].var[vname].getValues()[:-1,:-1,:]
                        print(vname,nb)
                        print(vv.shape)
                        print(v[i0[i]:i1[i],j0[j]:j1[j],:].shape)
                        #v[irange,jrange,krange]=vv.copy()
                        #mfl.blk[0].data.append(var(size,nvar,name=vname,val=v.copy()))
                        #mfl.blk[nb].var[vname]=mfl.blk[0].data[-1]
                        #nvar+=1
            return mfl.clone()

        def clone(self):
            """Returns a clone of the flow object"""
            obj=copy.copy(self)
            return obj
        
            
        
