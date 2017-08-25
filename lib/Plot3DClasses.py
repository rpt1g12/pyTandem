# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 11:24:30 2016

@author: rpt1g12
"""
#from lib.myPlot3dOperator import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import struct
import copy as copy
from lib.myPlots import *
from scipy.interpolate import griddata  

import glob

import os

#Define a float32
float32=np.float32

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
    if (verbose):
        print(v)
    return v

def getFileNames(path=None,pattern='solT???.????.q'):
    """Get solution files with extension='.ext'
    
    Arguments:
        path (string): Path to look for files. If none use the working directory.
        pattern (string): String that resembles the solution file with character wild cards, i.e. solT???.????.q
    
    Return:
        sfiles (list(string)[nfiles]: List containing file names"""
    
    if path==None:
        path=os.getcwd()

    sfiles=glob.glob(path+pattern)
    nfiles=len(sfiles)

    for i in range(nfiles):
        sfiles[i]=sfiles[i].split('/')[-1]

    sfiles.sort()

    return sfiles

# Classes
class var(object):
    """Defines a variable class.

    Args:
        size (list[3*int]): contains grid sizes in each direction
        vid (int): variable id
        vname (char): variable name
        val (np.ndarray(size*float32)): contains variable values
    """
    def __init__(self,size,vid,name='var',val=[]):
        """Class constructor
        
            If val=[], val is set to zeros

        """
        self.id=vid
        self.vname=name
        self.size=size
        self.flowc=list(np.zeros(4))
        if len(val)==0:
            self.val=np.zeros(size,dtype=float32)
        else:
            self.val=val.copy()
           
    def rdVar(self,fh,lh,nvar,nb,lblk,ntvar=0,Type='grid'):
        """Reads variable values from file to self.val
        
        Args:
            fh (file): file to read from
            lh (int): length of header in bytes
            nvar (int): variable number within the file
            nb (int): block where value is to be read from
            lblk (list[nbk*int]): Total number of points per each block within the file
            Type (char): can be either 'grid'(default) or 'sol'. Depending on that the
                is treated as a Plot3D grid or solution file

        """
        lxi=self.size[0];let=self.size[1];lze=self.size[2]
        ltomb=(lxi)*(let)*(lze)
        k8=np.float32; lk8=np.dtype(k8).itemsize
        if Type=='grid':
            tvar=3
            fl_hdr=0*lk8
        if Type=='fun':
            tvar=ntvar
            fl_hdr=0*lk8
        if Type=='sol':
            tvar=5
            fl_hdr=4*lk8
            fh.seek(lh)
            mach,aoa,Re,time=rdType(fh,'f',4)
        
        lh+=fl_hdr*(nb)+sum(lblk[:nb])*lk8*tvar
        if Type=='sol':
            fh.seek(lh)
            flowc=np.fromfile(fh,dtype=k8,count=4)
            self.flowc=flowc.copy()
            
            
        lh+=fl_hdr+nvar*lblk[nb]*lk8
        fh.seek(lh)
        self.val=np.fromfile(fh,dtype=k8,count=ltomb)
        self.val=np.reshape(self.val,(lxi,let,lze),order='F')
        if Type=='sol':
            return mach,aoa,Re,time

    
    def clone(self):
        """Clones the variable.

            Returns:
                obj (var): copy of the variable class to be returned

        """
        obj=copy.copy(self)
        return obj
        
    def wrVar(self,fh,lh,nvar,nb,lblk,ntvar=0,Type='grid',flInfo=None):
        """Writes variable to file from self.val.

            Args:
            fh (file): file to wrie to
            lh (int): length of header in bytes
            nvar (int): variable number within the file
            nb (int): block where value is to be read from
            lblk (list of ints): Total number of points per each block within the file
            Type (char): can be either 'grid'(default) or 'sol'. Depending on that the
                is treated as a Plot3D grid or solution file
                
        """
        k8=np.float32; lk8=np.dtype(k8).itemsize
        if Type=='grid':
            tvar=3
            fl_hdr=0*lk8
        if Type=='fun':
            tvar=ntvar
            fl_hdr=0*lk8
        if Type=='sol':
            tvar=5
            fl_hdr=4*lk8
        
        lh+=fl_hdr*(nb)+sum(lblk[:nb])*lk8*tvar
        if Type=='sol':
            fh.seek(lh)
            if flInfo==None:
                flowc=self.flowc.copy()
            else:
                flowc=flInfo
            fhdr=np.float32(np.array(flowc))
            fh.write(fhdr)
        
        lh+=fl_hdr+nvar*lblk[nb]*lk8
        #if Type=='sol': print('var',nb,nvar,lh)
        fh.seek(lh)
        v=(np.float32(np.transpose(self.getValues()).copy(order='C')))
        fh.write(v)
    
    def derVar(self,direction):
        """Derives variable with respect to any computational grid direction

            Args:
                direction (int): direction of derivative. Can be 0, 1 or 2.
                    Indicates the index with respect to which the derivative
                    is performed.
            Returns:
                dvar (np.ndarray(float32)): variable derivative

        """
        dvar=np.zeros(self.size,dtype=float32)
        n=self.size[direction]
        dirnames=['xi','eta','zeta']
        if (n<3):
            dvar[:,:,:]=1.0
            print('Not enough points to perform derivation in {} direction'.format(dirnames[direction]))
        else:
            dx=np.zeros(n,dtype=float32)
            a=np.zeros(n-2,dtype=float32);b=np.zeros(n-1,dtype=float32)
            c=np.zeros(n,dtype=float32)
            d=np.zeros(n-1,dtype=float32);e=np.zeros(n-2,dtype=float32);
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
            dvar*=0.5
        return dvar.copy()
    
    def avgDir(self,direction=2):
        """Averages variable in any computational grid direction.

            Args:
                direction (int): Direction of averaging

            Returns:
                aarr (np.ndarray(float32)): Averaged variable
                
        """
        size=self.size        
        if direction==0:
            lsize=size[1:]
            aarr=np.zeros(lsize,dtype=float32)
            arr=self.getValues()
            for i in range(size[0]):
                aarr+=arr[i,:,:]
        if direction==1:
            lsize=(size[0],size[2])
            aarr=np.zeros(lsize,dtype=float32)
            arr=self.getValues()
            for j in range(size[1]):
                aarr+=arr[:,j,:]
        if direction==2:
            lsize=size[0:2]
            aarr=np.zeros(lsize,dtype=float32)
            arr=self.getValues()
            for k in range(size[2]):
                aarr+=arr[:,:,k]
        aarr/=size[direction]
        return aarr.copy()
    
    def getValues(self,oned=False):
        """Extract variable values.

            Args:
                oned (int): If True returns a one-dimensional array. If False(default)
                    returns a three-dimensional array.
            Returns:
                self.val (np.ndarray(float32): Variable values

        """
        if oned:
           values=np.zeros(self.val.size)
           ii=0;lxi=self.size[0];let=self.size[1];lze=self.size[2]
           for k in range(lze):
               for j in range(let):
                   for i in range(lxi):
                        values[ii]=self.val[i,j,k]
                        ii+=1
        else:
           values=self.val.copy()
        return values
    def getSize(self):
        """Gets variable size."""
        return self.size.copy()
    def getName(self):
        """Gets variable name."""
        return self.vname
        
class blk(object):
    """Provides block class
        Parameters:
            id (int): Block id
            size (list[3*int]): grid dimensions in each direction
            glen (int): Total number of points 
            var (dict{var}): Dictionary of var objects
            data (list[var]): Array of var objects
            metFlag (logical): Flag stating if grid mettrics have been 
                already computed
                metFlag (logical): Flag stating if grid mettrics have been 
                    already computed. Default is False.

    """
    def __init__(self,blk_id,size=(1,1,1)):
        self.id=blk_id
        self.size=size
        self.glen=size[0]*size[1]*size[2]
        self.var={}
        self.data=[]
        self.metFlag=False
        
    def getMetrics(self):
        """Compute grid mettrics"""
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
            
    def derive(self,varName1,varName2,r=False):
        """Derive one variable with respect to x, y or z
            
            Args:
                varName1 (char): Variable to be derived
                varName2 (char): Variable with respect to be derived. It must be
                    either 'x', 'y', 'z', 'xi', 'et', or 'ze'.
                r (logic): If true returns the derivative array

        """
        if varName1 in self.var:
            var0=self.var[varName1]
            if varName2 not in ['xi','et','ze']:
                self.getMetrics()
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
            elif varName2 =='xi':           
                dv=var0.derVar(0)
            elif varName2 =='et':           
                dv=var0.derVar(1)
            elif varName2 =='ze':           
                dv=var0.derVar(2)
            else:
                print('Needs to be derived with respect to: x, y, z or xi, et, ze')
            if  r:
                self.setData(vname='d{}d{}'.format(varName1,varName2),val=dv)
                return dv
            else:
                self.setData(vname='d{}d{}'.format(varName1,varName2),val=dv)
        else:
            print('{} is not a valid variable!'.format(varName1))


    def getStrain(self,uvars=['u','v','w'],xvars=['x','y','z'],invariants=False,returnA=False):
        """docstring for getStrain"""
        self.getMetrics()    
        nxi,net,nze=self.size
        A=np.zeros((9,nxi,net,nze))
        ii=0
        for uvar in uvars:
            jj=0
            for xvar in xvars:
                print(ii*3+jj)
                A[ii*3+jj]=self.derive(uvar,xvar,r=True)
                jj+=1
            ii+=1
        if invariants:
            P=A[0]+A[4]+A[8]
            print('P done!')
            Q=A[0]*(A[4]+A[8])+A[4]*A[8]-A[1]*A[3]-A[2]*A[6]-A[5]*A[7]
            print('Q done!')
            R=A[0]*(A[4]*A[8]-A[5]*A[7])+A[1]*(A[5]*A[6]-A[3]*A[8])+A[2]*(A[3]*A[7]-A[4]*A[6])
            print('R done!')
            
            self.setData(vname='P',val=P)
            self.setData(vname='R',val=-R)
            self.setData(vname='Q',val=Q)

        if returnA:
            return A

    def getNormals(self,plane=1,flip=False,tangent=False,returnN=False):
        """docstring for getNormals"""
        self.getMetrics()
        if  plane==1:
            if flip:
                 fctr=-1/np.sqrt((self.mets[3].getValues())**2+(self.mets[4].getValues())**2+(self.mets[5].getValues())**2)
            else:
                 fctr=1/np.sqrt((self.mets[3].getValues())**2+(self.mets[4].getValues())**2+(self.mets[5].getValues())**2)
            nx=fctr*self.mets[3].getValues()[:,:,:]
            self.setData(vname='nx',val=nx)
            if tangent:
             self.setData(vname='ty',val=-nx)

            ny=fctr*self.mets[4].getValues()[:,:,:]
            self.setData(vname='ny',val=ny)
            if tangent:
             self.setData(vname='ty',val=ny)

            nz=fctr*self.mets[5].getValues()[:,:,:]
            self.setData(vname='nz',val=nz)
            if tangent:
             self.setData(vname='tz')
        
        if returnN:
            return nx,ny,nz

    def getWSS(self,plane=1,Re=1,flip=False):
        """docstring for getWSS"""
        self.getViscosity(tflag=False,Re=Re)
        nu=self.var['nu'].getValues()
        nx,ny,nz=self.getNormals(plane=plane,flip=flip,returnN=True)
        A=self.getStrain(returnA=True)

        tau=np.zeros_like(A)
        tw=np.zeros_like(A[0])

        dilatation=(-2.0/3)*(A[0]+A[4]+[8])

        for i in range(3):
            for j in range(3):
                ij=i*3+j
                ji=j*3+i
                if i==j:
                    tau[ij]=dilatation
                tau[ij]+=(A[ij]+A[ji])
                tau[ij]*=nu

        
        xvar=['x','y','z']
        for i in range(3):
            j=i*3
            tw=tau[j]*nx+tau[j+1]*ny+tau[j+2]*nz
            self.setData(vname='tw'+xvar[i],val=tw)
        pass

    def interpolate(self,vname,xi,yi,zi,method='nearest',rblock=False):
        """Interpolate variable at some points

                Args:
                    vname (char): Varible's name to be interpolated
                    xi (np.ndarray(size*float32): X coordinate of interpolation points
                    yi (np.ndarray(size*float32): Y coordinate of interpolation points
                    zi (np.ndarray(size*float32): Z coordinate of interpolation points
                    method (char): Available methods: 'nearest' (default) or 'linear'

                Returns:
                    iblock (block): Block containing interpolated variables

        """
        x=self.var['x'].getValues(oned=True)
        y=self.var['y'].getValues(oned=True)
        z=self.var['z'].getValues(oned=True)
        values=self.var[vname].getValues(oned=True)
        points=np.array(np.transpose([x,y,z]))
        ipoints=np.array(np.transpose([np.reshape(xi,xi.size),np.reshape(yi,yi.size),np.reshape(zi,zi.size)]))
        isize=list(xi.shape)
        print('Interpolating...')
        ival=griddata(points,values,ipoints,method=method)
        if  rblock:
            iblock=blk(self.id,size=isize)
            iblock.setData(vname='x',val=xi)
            iblock.setData(vname='y',val=yi)
            iblock.setData(vname='z',val=zi)
            iblock.setData(vname,val=np.reshape(ival,isize))
            return iblock.clone()
        else:
            return np.reshape(ival,isize)


    def interpolate2dk(self,vname,xi,yi,ki,method='cubic',mode='block'):
        """Interpolate data form variable (vname) into points (ipts)
           at ki span-plane.
                Args:
                    vname (char): Varible's name to be interpolated
                    xi (np.ndarray(size*float32): X coordinate of interpolation points
                    yi (np.ndarray(size*float32): Y coordinate of interpolation points
                    ki (int): Zeta plane on to which perform interpolation
                    method (char): Available methods: 'nearest' (default) or 'linear'
                    mode (char): If mode='block' returns a block object. 
                                 else: returns just a numpy array
                Returns:
                    iblock (blk): if mode='block'
                    iblock (np.ndarray[size*float32]): if mode='values'

                    """
        size=self.size[0]*self.size[1]
        x=self.var['x'].val[:,:,ki];x=np.reshape(x,x.size)
        y=self.var['y'].val[:,:,ki];y=np.reshape(y,y.size)
        z=self.var['z'].val[:,:,ki];z=np.reshape(z,z.size)
        points=np.array(np.transpose([x,y]))
        values=self.var[vname].val[:,:,ki];values=np.reshape(values,size)
        ipoints=np.array(np.transpose([np.reshape(xi,xi.size),np.reshape(yi,yi.size)]))
        if len(xi.shape)>1:
            isize=[xi.shape[0],xi.shape[1],1]
        else:
            isize=[xi.shape[0],1,1]
        print('Interpolating...')
        ival=griddata(points,values,ipoints,method=method)
        if  mode=='block':
            iblock=blk(self.id,size=isize)
            iblock.setData(vname='x',val=xi)
            iblock.setData(vname='y',val=yi)
            zi=xi.copy();zi[:,:]=z[0]
            iblock.setData(vname='z',val=zi)
            iblock.setData(vname,val=np.reshape(ival,isize))
        elif mode=='data':
            if len(xi.shape)>1:
                iblock=np.reshape(ival,isize)
            else:
                iblock=ival
        else:
            print("Invalid mode! select mode='block' or mode='data'")
            return
        
        return iblock

    def getSubset(self,xlim=None,ylim=None,zlim=None,link=False):
        """Extracts a subset of the block data

                Args:
                    xlim (range): Range xi planes from which to extract data
                        If none, uses all available
                    ylim (range): Range eta planes from which to extract data
                        If none, uses all available
                    zlim (range): Range zeta planes from which to extract data
                        If none, uses all available

                Returns:
                    sblk (blk): Block object containing subset

        """
        ndata=len(self.data)
        size=self.size
        if xlim==None:
            xlim=range(0,size[0])
        if ylim==None:
            ylim=range(0,size[1])
        if zlim==None:
            zlim=range(0,size[2])
        lsize=(len(xlim),len(ylim),len(zlim))
        sarr=np.zeros(lsize,dtype=float32)
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
        
        if link:
            return sblk
        else:
            return sblk.clone()
        
    def setData(self,vname=None,val=[],vid=None):
        """Sets values for a given variable in the block. If variable exists, overwrites its values.
            If variable does not exist yet it appeds a new one.

            Args:
                vname (char): Variable's name
                val (np.ndarray[size*float32]): Variable values
        """
        if vname==None:
            vname='v{:d}'.format(len(self.data))
        
        if len(val)==0:
            size=self.size
            val=np.zeros(size,dtype=float32)
        else:
            size=val.shape
           
        if (vname in self.var):
            print('Changing {} values'.format(vname))
            self.var[vname].val[:,:,:]=val.copy()
        else:
            print('Adding {} to variable list'.format(vname))
            self.data.append(var(size,len(self.data),vname,val))
            self.var[vname]=self.data[-1]

    def contourf(self,varname,vmin=None,vmax=None,k=0,plane=2,ax=None,nlvl=11,avg=False,bar=True,cmap=None):
        """Produces a contour plot of the variable varname"""
        if (ax==None):
            f,ax=getFig(varname)
        else:
            f=ax.figure
        if avg:
            if plane==2:
               x_char,y_char='x','y'
            elif plane==0:
               x_char,y_char='y','z'
            elif plane==1:
               x_char,y_char='x','z'
            x=self.var[x_char].avgDir(plane)
            y=self.var[y_char].avgDir(plane)
            v=self.var[varname].avgDir(plane)
        else:
            if plane==2:
               x_char,y_char='x','y'
               x=self.var[x_char].val[:,:,k]
               y=self.var[y_char].val[:,:,k]
               v=self.var[varname].val[:,:,k]
            elif plane==0:
               x_char,y_char='y','z'
               x=self.var[x_char].val[k,:,:]
               y=self.var[y_char].val[k,:,:]
               v=self.var[varname].val[k,:,:]
            elif plane==1:
               x_char,y_char='x','z'
               x=self.var[x_char].val[:,k,:]
               y=self.var[y_char].val[:,k,:]
               v=self.var[varname].val[:,k,:]
        if cmap==None:
            cmap=plt.cm.coolwarm

        if vmin==None:
            vmin=np.min(v)
        if vmax==None:
            vmax=np.max(v)

        print(vmin,vmax)

        lvl=np.linspace(vmin,vmax,nlvl)
            
        im=ax.contourf(x,y,v,levels=lvl,cmap=cmap,vmin=vmin,vmax=vmax,extend='both')
        if bar:
            cb=f.colorbar(im,orientation='vertical')
            cb.set_label(varname,fontsize=20)
        ax.set_xlabel(r'$'+x_char+'$',fontsize=20)
        ax.set_ylabel(r'$'+y_char+'$',fontsize=20)
        if bar:
            axShow(ax)
        
        return f,ax,im

    def contour(self,varname,vmin=None,vmax=None,k=0,plane=2,ax=None,nlvl=11,avg=False,colors=None,dashed=False):
        """Produces a contour plot of the variable varname"""
        if (ax==None):
            f,ax=getFig(varname)
        else:
            f=ax.figure
        if avg:
            if plane==2:
               x_char,y_char='x','y'
            elif plane==0:
               x_char,y_char='y','z'
            elif plane==1:
               x_char,y_char='x','z'
            x=self.var[x_char].avgDir(plane)
            y=self.var[y_char].avgDir(plane)
            v=self.var[varname].avgDir(plane)
        else:
            if plane==2:
               x_char,y_char='x','y'
               x=self.var[x_char].val[:,:,k]
               y=self.var[y_char].val[:,:,k]
               v=self.var[varname].val[:,:,k]
            elif plane==0:
               x_char,y_char='y','z'
               x=self.var[x_char].val[k,:,:]
               y=self.var[y_char].val[k,:,:]
               v=self.var[varname].val[k,:,:]
            elif plane==1:
               x_char,y_char='x','z'
               x=self.var[x_char].val[:,k,:]
               y=self.var[y_char].val[:,k,:]
               v=self.var[varname].val[:,k,:]
        if colors==None:
            colors='black'
        if dashed==False:
            #Set negative contours to be solid instead of dashed:
            matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

        if vmin==None:
            vmin=np.min(v)
        if vmax==None:
            vmax=np.max(v)

        print(vmin,vmax)

        lvl=np.linspace(vmin,vmax,nlvl)
            
        im=ax.contour(x,y,v,levels=lvl,vmin=vmin,vmax=vmax,extend='both',colors=colors)
        ax.set_xlabel(r'$'+x_char+'$',fontsize=20)
        ax.set_ylabel(r'$'+y_char+'$',fontsize=20)
        axShow(ax)
        
        return f,ax,im

    def drawMeshPlane(self,direction=2,pln=0,skp=(1,1),ax=None,color='black',lw=1,showBlock=False,hideAx=False):
        """Plots one plane of the grid
            Args:
                direction (int): 0->yz-plane; 1->xz-plane; 2->xy-plane
                pln (int): Plane of constant xi,eta or zeta to take slice from
                skp (int): Skip index. Plot every skp line.
                ax (plt.axis): axis used to plot. Default is None, so it uses the current axis.
                color (string): Color of the mesh.
                lw (int): Width of the lines
                showBlock (boolean): If True show block boundaries in blue

            Returns:
                Nothing
        """
        skp=np.asarray(skp)

        if len(skp)==1:
            a=np.asarray([skp,skp])
            skp=a.copy()
            del a
        if ax==None:
            f,ax=getFig('Mesh')
        else:
            f=ax.figure


        if direction==0:
            size=[int(np.ceil(self.size[1]/skp[0])),int(np.ceil(self.size[2]/skp[1])),2]
            grid=np.zeros(size)
            grid[:,:,0]=self.var['y'].getValues()[pln,0::skp[0],0::skp[1]]
            grid[:,:,1]=self.var['z'].getValues()[pln,0::skp[0],0::skp[1]]
        elif direction==1:
            size=[int(np.ceil(self.size[0]/skp[0])),int(np.ceil(self.size[2]/skp[1])),2]
            grid=np.zeros(size)
            grid[:,:,0]=self.var['x'].getValues()[0::skp[0],pln,0::skp[1]]
            grid[:,:,1]=self.var['z'].getValues()[0::skp[0],pln,0::skp[1]]
        elif direction==2:
            size=[int(np.ceil(self.size[0]/skp[0])),int(np.ceil(self.size[1]/skp[1])),2]
            grid=np.zeros(size)
            grid[:,:,0]=self.var['x'].getValues()[0::skp[0],0::skp[1],pln]
            grid[:,:,1]=self.var['y'].getValues()[0::skp[0],0::skp[1],pln]
        
        for j in range(size[1]):
            ax.plot(grid[:,j,0],grid[:,j,1],lw=lw,color=color)
        for i in range(size[0]):
            ax.plot(grid[i,:,0],grid[i,:,1],lw=lw,color=color)
        if showBlock:
            ax.plot(grid[:,0,0],grid[:,0,1],lw=lw+1,color='blue')
            ax.plot(grid[0,:,0],grid[0,:,1],lw=lw+1,color='blue')
            ax.plot(grid[:,-1,0],grid[:,-1,1],lw=lw+1,color='blue')
            ax.plot(grid[-1,:,0],grid[-1,:,1],lw=lw+1,color='blue')
        if hideAx:
            ax.axes.get_yaxis().set_visible(False)
            ax.axes.get_xaxis().set_visible(False)
        ax.set_aspect('equal')
        fit(ax,(0.1,0.1))
        if direction==1:
            ax.invert_yaxis()
        f.canvas.manager.window.raise_()
        axShow(ax)
        return f,ax

    def wrVarASCII(self,varnames=None,k=0,fpath=None,direction=2):
        """Writes values of the variables contained in varnames in ASCII format"""
        if varnames==None:
            print('Error! You need to specify varnames')
            return
        if fpath==None:
            print('Error! You need to specify fpath')
            return
        fname=''
        cvariables='VARIABLES='
        for c in varnames:
            fname+=c
            cvariables+='"{}", '.format(c)
        ctitle='TITLE="{} TecPlot file"\n'.format(fname)
        cvariables+='\n'
        fname+='.dat'
        fpath+=fname
        nvar=len(varnames)
        if direction==0:
            nt=self.size[2]*self.size[1]
        elif direction==1:
            nt=self.size[0]*self.size[2]
        elif direction==2:
            nt=self.size[0]*self.size[1]
        csize='ZONE I={:d}, J={:d}, DATAPACKING=POINT\n'.format(self.size[0],self.size[1])
        aout=np.zeros((nvar,nt))


        i=0
        for c in varnames:
            if direction==0:
                aout[i,:]=np.reshape(self.var[c].getValues()[k,:,:],nt,order='F')
            elif direction==1:
                aout[i,:]=np.reshape(self.var[c].getValues()[:,k,:],nt,order='F')
            elif direction==2:
                aout[i,:]=np.reshape(self.var[c].getValues()[:,:,k],nt,order='F')
            i+=1

        fh=open(fpath,'w')
        s=ctitle+cvariables+csize
        for i in range(nt):
            for j in range(nvar):
                s+='{:12.5e}\t'.format(aout[j,i])
            s+='\n'
        print('Writing to {}'.format(fpath))
        fh.write(s)
        fh.close()
        pass

    def wrVarASCII3D(self,varnames=None,fpath=None):
        """Writes values of the variables contained in varnames in ASCII format"""
        if varnames==None:
            print('Error! You need to specify varnames')
            return
        if fpath==None:
            print('Error! You need to specify fpath')
            return
        fname=''
        cvariables='VARIABLES='
        for c in varnames:
            fname+=c
            cvariables+='"{}", '.format(c)
        ctitle='TITLE="{} TecPlot file"\n'.format(fname)
        cvariables+='\n'
        fname+='_3D.dat'
        fpath+=fname
        nvar=len(varnames)
        nt=self.glen
        csize='ZONE I={:d}, J={:d}, K={:d}, DATAPACKING=POINT\n'.format(self.size[0],self.size[1],self.size[2])
        aout=np.zeros((nvar,nt))


        i=0
        for c in varnames:
            aout[i,:]=self.var[c].getValues(True)
            i+=1

        fh=open(fpath,'w')
        s=ctitle+cvariables+csize
        for i in range(nt):
            for j in range(nvar):
                s+='{:12.5e}\t'.format(aout[j,i])
            s+='\n'
        print('Writing to {}'.format(fpath))
        fh.write(s)
        fh.close()
        pass
            
    def getViscosity(self,tflag=True,Re=1):
        """Computes viscosity based on Sutherlands law and stores in a new variable 'nu'
        It also stores Temperature if tflag==True
        """
        rhoI=1/self.var['r'].getValues()
        p=self.var['p'].getValues()
        T=1.4*p*rhoI #T=\gamma*p/\rho (non-dimensional)
        S=111/273 #S=S0/T0 (non-dimensional) 
        C1=(1+S) # Sutherlands constant C1=(T0/T0+S0)*(T0/T0)^{-1.5} (non-dimensional)
        if Re==1:
            nu=(C1*T**(1.5)/(T+S))
        else:
            nu=(C1*T**(1.5)/(T+S))/Re

        if  tflag:
            self.setData(vname='T',val=T)

        self.setData(vname='nu',val=nu)
        pass

    def setTouch(self,x0=0,y0=0,R=0.001,H=3e-8,acc=1e-6):
        """docstring for setTouch"""
        k0=np.log(acc)/(R**2)
        xo=np.sqrt(-1/(2*k0))
        f0=xo*np.exp(-0.5)
        Gamma=H/f0
        x=self.var['x'].getValues()-x0
        y=self.var['y'].getValues()-y0
        r=x**2+y**2
        A=Gamma*np.exp(k0*r)
        fu=-A*y
        fv=A*x
        fU=(fu**2+fv**2)**0.5
        self.setData(vname='fu',val=fu)
        self.setData(vname='fv',val=fv)
        self.setData(vname='fU',val=fU)
        self.setData(vname='FA',val=A)

    def clone (self):
        obj=copy.copy(self)
        return obj 
        
class flow(object):
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
         * vnames (list): Variable names
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
            #Variable names
            self.vnames=[]
            self.gnames=[]
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
            print('Reading header...')
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

        def setHdr(self,nbks,lxib,letb,lzeb):
            """Sets up the flow object and the header file"""
            #Compute total number of blocks
            self.nbk=nbks[0]*nbks[1]*nbks[2]
            #Read Grid sizes for each block
            for k in range(nbks[2]):
                for j in range(nbks[1]):
                    for i in range(nbks[0]):
                        nb=i+j*nbks[0]+k*nbks[1]*nbks[2]
                        self.blk.append(blk(blk_id=nb))
                        self.blk[nb].size=[lxib[i],letb[j],lzeb[k]]
                        #Compute total size of arrays
                        for n in range(3):
                            self.blk[nb].glen*=self.blk[nb].size[n]
            #Store header length
            self.lhdr=4*(1+self.nbk*3)
            #Store total sizes for all blocks
            for nb in range(self.nbk):
                self.lblk.append(self.blk[nb].glen)
            pass

        def wrGrid(self,gfile=None):
            """Writes grid file"""
            if gfile==None:
                gfile=self.gfile
            if not os.path.exists(self.path):
                os.makedirs(self.path)
            fh=open(self.path+gfile,'wb')
            nbk=self.nbk
            lblk=self.lblk
            ghdr=[nbk]
            gnames=['x','y','z']
            self.gnames=gnames
            print('Writing grid: {}'.format(gfile))
            for nb in range(nbk): 
                lh=self.lhdr
                ghdr.append(self.blk[nb].size[0])
                ghdr.append(self.blk[nb].size[1])
                ghdr.append(self.blk[nb].size[2])
                for nvar in range(3):
                    gname=gnames[nvar]
                    self.blk[nb].var[gname].wrVar(fh,lh,nvar,nb,lblk,Type='grid')
            fh.seek(0)
            fhdr=np.int32(np.array(ghdr))
            fh.write(fhdr)
            fh.close()
            pass

        def rdGrid(self):
            #Open file for reading in binary
            fh=open(self.path+self.gfile,'rb')
            #Number of blocks
            nbk=self.nbk
            #Variable names
            names=['x','y','z']
            self.gnames=names
            #Number of points per block
            lblk=self.lblk
            print('Reading grid: {}'.format(self.gfile))
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
        
        def rdFlowInfo(self,sfile=None):
            """Read Flow Operating Conditions"""
            if sfile!=None:
                self.sfile=sfile
            fh=open(self.path+self.sfile,'rb')
            lh=self.lhdr
            fh.seek(lh)
            mach,aoa,Re,time=rdType(fh,'f',4)
            self.mach=mach
            self.aoa=aoa
            self.Re=Re
            self.time=time
            fh.close()
            return mach,aoa,Re,time

        def rdSol(self,vnames=None,sfile=None):
            if sfile!=None:
                self.sfile=sfile
            fh=open(self.path+self.sfile,'rb')
            nbk=self.nbk
            if vnames==None:
                names=[]
                for n in range(5):
                    names.append('v{:d}'.format(n))
            else:
                names=vnames
            lblk=self.lblk
            print('Reading solution: {}'.format(self.sfile))
            for nb in range(nbk):
                size=self.blk[nb].size 
                lh=self.lhdr

                for nvar in range(5):
                    vname=names[nvar]
                    if len(self.blk[nb].data)<8:
                        self.blk[nb].data.append(var(size,nvar,name=vname))
                        nn=-1
                    else:
                        nn=nvar+3
                    self.blk[nb].data[nn].rdVar(fh,lh,nvar,nb,lblk,Type='sol')
                    self.blk[nb].var[vname]=self.blk[nb].data[nn]
                    
            self.vnames=names
            fh.close()

        def rdFun(self,vnames=None,ffile=None):
            """Reads a Plot3D function file"""
            if ffile==None:
                print('Must provide a function file name!!')
                pass
            lblk=self.lblk
            nbk=self.nbk
            k8=np.float32; lk8=np.dtype(k8).itemsize
            lh=self.lhdr+nbk*lk8
            fh=open(self.path+ffile,'rb')
            fh.seek(lh-lk8)
            ntvar=rdType(fh,'i',1)
            print(ntvar)
            if vnames==None:
                names=[]
                for n in range(ntvar):
                    names.append('v{:d}'.format(n))
            else:
                names=vnames
            print('Reading function file: {}'.format(ffile))
            for nb in range(nbk):
                size=self.blk[nb].size 
                lh=self.lhdr
                for nvar in range(ntvar):
                    vname=names[nvar]
                    self.blk[nb].data.append(var(size,nvar,name=vname));nn=-1
                    self.blk[nb].data[nn].rdVar(fh,lh,nvar,nb,lblk,ntvar,Type='fun')
                    self.blk[nb].var[vname]=self.blk[nb].data[nn]
                    
            #self.vnames=names
            fh.close()

            pass

        def wrFun(self,vnames=None,ffile=None,path=None):
            """Writes Plot3D function files"""
            k8=np.float32; lk8=np.dtype(k8).itemsize
            if not os.path.exists(self.path):
                os.makedirs(self.path)
            fh=open(self.path+ffile,'wb')
            nbk=self.nbk
            lblk=self.lblk
            ghdr=[nbk]
            ntvar=len(vnames)
            for nb in range(nbk): 
                lh=self.lhdr+nbk*lk8
                ghdr.append(self.blk[nb].size[0])
                ghdr.append(self.blk[nb].size[1])
                ghdr.append(self.blk[nb].size[2])
                ghdr.append(ntvar)
                for nvar in range(ntvar):
                    name=vnames[nvar]
                    self.blk[nb].var[name].wrVar(fh,lh,nvar,nb,lblk,ntvar,Type='fun')
            fh.seek(0)
            fhdr=np.int32(np.array(ghdr))
            fh.write(fhdr)
            fh.close()
            pass

        def wrSol(self,vnames=None,sfile=None,path=None,flInfo=np.zeros(4)):
            if path==None:
                path=self.path
            if not os.path.exists(path):
                os.makedirs(path)
            if sfile==None:
                sfile=self.sfile
                
            self.vnames=list(self.blk[0].var.keys())
            if vnames==None:
                vnames=[name for name in self.vnames[:6]]
                vnames.sort()
            else:
                if (len(vnames)>5):
                    print('WARNING: Maximum number of variables exceed!')
                    print('Only the first 5 variables will be writen.')
                    vnames=[name for name in self.vnames[:6]]
                elif (len(vnames)<5):
                    print('WARNING: Minimum number of variables not reached!')
                    print('Aborting...')
                    return
            
            
            fh=open(path+sfile,'wb')
            nbk=self.nbk
            lblk=self.lblk
            ghdr=[nbk]
            print('Writing solution: {}'.format(sfile))
            print('variables: {}, {}, {}, {}, {}'.format(vnames[0],
                            vnames[1],vnames[2],vnames[3],vnames[4]))
            for nb in range(nbk): 
                lh=self.lhdr
                ghdr.append(self.blk[nb].size[0])
                ghdr.append(self.blk[nb].size[1])
                ghdr.append(self.blk[nb].size[2])
                for nvar in range(5):
                    name=vnames[nvar]
                    self.blk[nb].var[name].wrVar(fh,lh,nvar,nb,lblk,Type='sol',flInfo=flInfo)
            fh.seek(0)
            fhdr=np.int32(np.array(ghdr))
            fh.write(fhdr)
            fh.close()

        def shiftK(self,k,vnames=None):
            """Shifts the flow object in the ze direction so the kth point is set to be 0 in the returning flow object
               Returns a flow object"""
            # Clone the current flow object
            sflow=self.clone()
            print('Shifting data by {:d}'.format(k))
            if vnames==None:
                names=sflow.vnames
            else:
                names=vnames
            for nb in range(sflow.nbk):
                for n in names:
                    vr2=self.blk[nb].var[n].getValues()
                    lze=sflow.blk[nb].size[2]
                    lze_1=lze-1
                    sflow.blk[nb].var[n].val[:,:,0]=vr2[:,:,k].copy()
                    sflow.blk[nb].var[n].val[:,:,lze_1]=vr2[:,:,k].copy()
                    for sk in range(1,lze_1):
                        kk=sk+(k)
                        if (kk>(lze_1)):
                            kk-=lze_1
                        sflow.blk[nb].var[n].val[:,:,sk]=vr2[:,:,kk].copy()
                        
            return sflow
        
        def mergeBlocks(self,bkx,bky):
            """Merges all blocks into a single one
               Returns a flow object"""
            lxi=0;let=0;lze=self.blk[0].size[2]
            i1=[];i0=[];j0=[];j1=[]
            for i in range(bkx):
                i0.append(lxi)
                lxi+=self.blk[i].size[0]
                i1.append(lxi)
            for j in range(0,self.nbk,bkx):
                j0.append(let)
                let+=self.blk[j].size[1]
                j1.append(let)
            print(i0,i1,lxi)
            print(j0,j1,let)
            v=np.zeros((lxi,let,lze),dtype=float32)
            size=[lxi,let,lze]
            mfl=flow(self.path,'merged_{}'.format(self.gfile),'merged_{}'.format(self.sfile))
            mfl.nbk=1
            mfl.lhdr=16
            mfl.lblk.append(lxi*let*lze)
            mfl.vnames=self.vnames.copy()
            mfl.blk.append(blk(blk_id=0))
            mfl.blk[0].size=size
            mfl.blk[0].glen=lxi*let*lze
            mfl.vnames.insert(0,'z')
            mfl.vnames.insert(0,'y')
            mfl.vnames.insert(0,'x')
            print('Merging Blocks...')
            for vname in mfl.vnames:
                nvar=0
                for j in range(bky):
                    for i in range(bkx):
                        nb=j*(bkx)+i
                        vv=self.blk[nb].var[vname].getValues()
                        v[i0[i]:i1[i],j0[j]:j1[j],:]=vv.copy()
                mfl.blk[0].data.append(var(size,nvar,name=vname,val=v.copy()))
                mfl.blk[0].var[vname]=mfl.blk[0].data[-1]
                nvar+=1
            return mfl.clone()

        def drawMeshPlane(self,direction=2,pln=0,skp=1,ax=None,color='black',lw=1,showBlock=False,hideAx=False):

            for i in range(self.nbk):
                self.blk[i].drawMeshPlane(direction,pln,skp,ax,color,lw,showBlock,hideAx)

        def contourf(self,varname,vmin=None,vmax=None,k=0,plane=2,ax=None,nlvl=11,avg=False,bar=True,cmap=None):
            """Produces a contour plot of the variable varname"""
            f,a,im=self.blk[0].contourf(varname,vmin,vmax,k,plane,ax,nlvl,avg,bar=False,cmap=cmap)
            for i in range(1,self.nbk-1):
                self.blk[i].contourf(varname,vmin,vmax,k,plane,a,nlvl,avg,bar=False,cmap=cmap)
            self.blk[-1].contourf(varname,vmin,vmax,k,plane,a,nlvl,avg,bar,cmap)
            
            return f,a,im

        def contour(self,varname,vmin=None,vmax=None,k=0,plane=2,ax=None,nlvl=11,avg=False,colors=None,dashed=False):
            """Produces a contour plot of the variable varname"""
            f,a,im=self.blk[0].contour(varname,vmin,vmax,k,plane,ax,nlvl,avg,colors,dashed)
            for i in range(1,self.nbk-1):
                self.blk[i].contour(varname,vmin,vmax,k,plane,a,nlvl,avg,colors,dashed)
            self.blk[-1].contour(varname,vmin,vmax,k,plane,a,nlvl,avg,colors,dashed)
            
            return f,a,im

        def setTouch(self,x0=0,y0=0,R=0.001,H=1,acc=1e-6):
            """docstring for setTouch"""
            for i in range(self.nbk):
                self.blk[i].setTouch(x0,y0,R,H,acc)
        
        def getMetrics(self):
            """Compute grid mettrics"""
            for i in range(self.nbk):
                print('Metrics block {:}'.format(i))
                self.blk[i].getMetrics()

        def getStrain(self,uvars=['u','v','w'],xvars=['x','y','z']):
            """docstring for getStrain"""
            for bk in range(self.nbk):
                print('Strain-Rate block {:}'.format(bk))
                self.blk[bk].getStrain(uvars,xvars)


        def getAvg(self,files,vnames,path=None,gfile=None,sfile=None,out=True):
            """Returns time average of files"""
            if path==None:
                path=self.path
            if gfile==None:
                gfile=self.gfile
            if sfile==None:
                sfile=self.sfile
            if out:
                flavg=flow(path,gfile,sfile)
                flavg.rdHdr()
                flavg.rdGrid()
                #flavg.rdSol(vnames,files[0])
                flavg.sfile='solTA.qa'
            else:
                flavg=self
            nt=len(files)

            for bk in range(self.nbk):
                for cvar in vnames:
                    flavg.blk[bk].setData(vname=cvar)
                for cvar in vnames:
                    flavg.blk[bk].setData(vname=cvar+cvar)
                for cvar in ['uv','uw','vw']:
                    flavg.blk[bk].setData(vname=cvar)

            t=np.zeros(len(files))
            for i in range(len(files)):
                f=files[i]
                t[i]=self.rdFlowInfo(f)[3]
            fctr=0.5/(t[-1]-t[0])
            dt=np.zeros_like(t)
            dt[0]=2*fctr*(t[1]-t[0])
            dt[-1]=2*fctr*(t[-1]-t[-2])
            for i in range(1,len(t)-1):
                dt[i]=fctr*(t[i+1]-t[i-1])

            
            ii=0
            for f in files:
                self.rdSol(vnames,f)
                for bk in range(self.nbk):
                    for cvar in vnames:
                        flavg.blk[bk].var[cvar].val[:,:,:]+=self.blk[bk].var[cvar].getValues()*dt[ii]
                        flavg.blk[bk].var[cvar+cvar].val[:,:,:]+=self.blk[bk].var[cvar].getValues()*self.blk[bk].var[cvar].getValues()*dt[ii]
                    for cvar in ['v','w']:
                        flavg.blk[bk].var['u'+cvar].val[:,:,:]+=self.blk[bk].var['u'].getValues()*self.blk[bk].var[cvar].getValues()*dt[ii]
                    flavg.blk[bk].var['vw'].val[:,:,:]+=self.blk[bk].var['v'].getValues()*self.blk[bk].var['w'].getValues()*dt[ii]

            #for bk in range(flavg.nbk):
            #    for cvar in vnames+['rr','uu','vv','ww','pp','uv','uw','vw']:
            #        flavg.blk[bk].var[cvar].val[:,:,:]=flavg.blk[bk].var[cvar].getValues()/nt

            for bk in range(flavg.nbk):
                    for cvar in vnames:
                        flavg.blk[bk].var[cvar+cvar].val[:,:,:]=flavg.blk[bk].var[cvar+cvar].getValues()-flavg.blk[bk].var[cvar].getValues()*flavg.blk[bk].var[cvar].getValues()
                    for cvar in ['v','w']:
                        flavg.blk[bk].var['u'+cvar].val[:,:,:]+=flavg.blk[bk].var['u'+cvar].getValues()-flavg.blk[bk].var['u'].getValues()*flavg.blk[bk].var[cvar].getValues()
                    flavg.blk[bk].var['vw'].val[:,:,:]+=flavg.blk[bk].var['vw'].getValues()-flavg.blk[bk].var['v'].getValues()*flavg.blk[bk].var['w'].getValues()

            if out:
                return flavg
            else:
                pass

        def getSubsets(self,fromBlocks=[0],bkXYZ=[1,1,1],xRanges=[None],yRanges=[None],zRanges=[None],ssName=None,sspath=None,gfile=None,ssfile=None,ssfl=None):
            """Extract subsets from solution"""

            blocks=fromBlocks

            if ssfl==None:
                if ssName==None:
                    ssName='extractedSS'
                if sspath==None:
                    sspath=self.path+ssName+'/'
                if gfile==None:
                    gfile=self.gfile
                if ssfile==None:
                    ssfile=self.sfile
                ss=flow(sspath,gfile,ssfile)
            else:
                ss=ssfl

            #%%SubSet set-up
            ssBlocks=[]
            nxi=[];neta=[];nzeta=[]
            for kb in range(bkXYZ[2]):    
                for jb in range(bkXYZ[1]):
                    for ib in range(bkXYZ[0]):
                        nb=kb*bkXYZ[1]*bkXYZ[0]+jb*bkXYZ[0]+ib
                        ssBlocks.append(self.blk[blocks[nb]].getSubset(xlim=xRanges[ib],ylim=yRanges[jb],zlim=zRanges[kb],link=True))
                        if kb==0:            
                            neta.append(ssBlocks[-1].size[1])
                        if jb==0:            
                            nxi.append(ssBlocks[-1].size[0])
                nzeta.append(ssBlocks[-1].size[2])

            if ssfl==None:
                ss.setHdr(bkXYZ,nxi,neta,nzeta)
                for nb in range(len(blocks)):
                    ss.blk[nb]=ssBlocks[nb]
            
                return ss
            else:
                return ssBlocks

        def clone(self):
            """Returns a clone of the flow object"""
            obj=copy.copy(self)
            return obj
        
class rake(object):
    """Defines a rake of points"""
    def __init__(self, xo,yo,dx,dy,n,l):
        dmod=np.sqrt(dx**2+dy**2)
        dl=l/(n-1)  
        self.dl=dl
        self.dx,self.dy=dl*dx/dmod,dl*dy/dmod        
        self.tx,self.ty=self.dy/dl,-self.dx/dl
        self.x,self.y,self.var=np.zeros((n,1)),np.zeros((n,1)),np.zeros((n,1))

        for i in range(n):
            self.x[i]=xo+i*self.dx
            self.y[i]=yo+i*self.dy
