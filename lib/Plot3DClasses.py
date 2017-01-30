# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 11:24:30 2016

@author: rpt1g12
"""
#from lib.myPlot3dOperator import *
import numpy as np
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
           
    def rdVar(self,fh,lh,nvar,nb,lblk,Type='grid'):
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
        
    def wrVar(self,fh,lh,nvar,nb,lblk,Type='grid'):
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
        if (n<3):
            dvar[:,:,:]=1.0
            print('Not enough points to perform derivation')
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
                    either 'x', 'y' or 'z'.
                r (logic): If true returns the derivative array

        """
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
            if  r:
                return dv
            else:
                self.setData(vname='d{}d{}'.format(varName1,varName2),val=dv)
        else:
            print('{} is not a valid variable!'.format(varName1))

    def interpolate(self,vname,xi,yi,zi,method='nearest'):
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
        iblock=blk(self.id,size=isize)
        iblock.setData(vname='x',val=xi)
        iblock.setData(vname='y',val=yi)
        iblock.setData(vname='z',val=zi)
        iblock.setData(vname,val=np.reshape(ival,isize))
        return iblock.clone()

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
                        If mode='values' returns just a numpy array
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
        isize=[xi.shape[0],xi.shape[1],1]
        print('Interpolating...')
        ival=griddata(points,values,ipoints,method=method)
        if  mode=='block':
            iblock=blk(self.id,size=isize)
            iblock.setData(vname='x',val=xi)
            iblock.setData(vname='y',val=yi)
            zi=xi.copy();zi[:,:]=z[0]
            iblock.setData(vname='z',val=zi)
            iblock.setData(vname,val=np.reshape(ival,isize))
        elif mode=='values':
            iblock=np.reshape(ival,isize)
        
        return iblock

    def getSubset(self,xlim=None,ylim=None,zlim=None):
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

    def contourf(self,varname,vmin=-1,vmax=1,k=0,ax=None,nlvl=11,avg=False,bar=True,cmap=None):
        """Produces a contour plot of the variable varname"""
        if (ax==None):
            f,ax=getFig(varname)
        else:
            f=ax.figure
        if avg:
            x=self.var['x'].avgDir(2)
            y=self.var['y'].avgDir(2)
            v=self.var[varname].avgDir(2)
        else:
            x=self.var['x'].val[:,:,k]
            y=self.var['y'].val[:,:,k]
            v=self.var[varname].val[:,:,k]
        if cmap==None:
            cmap=plt.cm.coolwarm

        lvl=np.linspace(vmin,vmax,nlvl)
            
        im=ax.contourf(x,y,v,levels=lvl,cmap=cmap,vmin=vmin,vmax=vmax,extend='both')
        if bar:
            cb=f.colorbar(im,orientation='vertical')
            cb.set_label(varname,fontsize=20)
        ax.set_xlabel(r'$x$',fontsize=20)
        ax.set_ylabel(r'$y$',fontsize=20)
        if bar:
            axShow(ax)
        return f,ax,im

    def drawMeshPlane(self,direction=2,pln=0,skp=1,ax=None,color='black',lw=1,showBlock=False,hideAx=False):
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
        if ax==None:
            ax=plt.gca()

        if direction==0:
            size=[int(np.ceil(self.size[1]/skp)),int(np.ceil(self.size[2]/skp)),2]
            grid=np.zeros(size)
            grid[:,:,0]=self.var['y'].getValues()[pln,0::skp,0::skp]
            grid[:,:,1]=self.var['z'].getValues()[pln,0::skp,0::skp]
        elif direction==1:
            size=[int(np.ceil(self.size[0]/skp)),int(np.ceil(self.size[2]/skp)),2]
            grid=np.zeros(size)
            grid[:,:,0]=self.var['x'].getValues()[0::skp,pln,0::skp]
            grid[:,:,1]=self.var['z'].getValues()[0::skp,pln,0::skp]
        elif direction==2:
            size=[int(np.ceil(self.size[0]/skp)),int(np.ceil(self.size[1]/skp)),2]
            grid=np.zeros(size)
            grid[:,:,0]=self.var['x'].getValues()[0::skp,0::skp,pln]
            grid[:,:,1]=self.var['y'].getValues()[0::skp,0::skp,pln]
        
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
        axShow(ax)
        pass

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
            
    def getViscosity(self,tflag=True):
        """Computes viscosity based on Sutherlands law and stores in a new variable 'nu'
        It also stores Temperature if tflag==True
        """
        rhoI=1/self.var['r'].getValues()
        p=self.var['p'].getValues()
        T=1.4*p*rhoI #T=\gamma*p/\rho (non-dimensional)
        S=111/273 #S=S0/T0 (non-dimensional) 
        C1=(1+S) # Sutherlands constant C1=(T0/T0+S0)*(T0/T0)^{-1.5} (non-dimensional)
        nu=C1*T**(1.5)/(T+S)

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

        def wrGrid(self):
            """Writes grid file"""
            fh=open(self.path+self.gfile,'wb')
            nbk=self.nbk
            lblk=self.lblk
            ghdr=[nbk]
            gnames=['x','y','z']
            self.gnames=gnames
            print('Writing grid: {}'.format(self.gfile))
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
                    self.blk[nb].data.append(var(size,nvar,name=vname))
                    self.blk[nb].data[-1].rdVar(fh,lh,nvar,nb,lblk,Type='sol')
                    self.blk[nb].var[vname]=self.blk[nb].data[-1]
                    
            self.vnames=names
            fh.close()

        def wrSol(self,sfile=None,vnames=None):
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

            fh=open(self.path+sfile,'wb')
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
                    self.blk[nb].var[name].wrVar(fh,lh,nvar,nb,lblk,Type='sol')
            fh.seek(0)
            fhdr=np.int32(np.array(ghdr))
            fh.write(fhdr)
            fh.close()

        def shiftK(self,k):
            """Shifts the flow object in the ze direction so the kth point is set to be 0 in the returning flow object
               Returns a flow object"""
            # Clone the current flow object
            sflow=self.clone()
            print('Shifting data by {:d}'.format(k))
            for nb in range(sflow.nbk):
                for n in sflow.vnames:
                    vr2=sflow.blk[nb].var[n].getValues()
                    lze=sflow.blk[nb].size[2]
                    lze_1=lze-1
                    for sk in range(lze):
                        kk=sk+(k-1)
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

        def contourf(self,varname,vmin=-1,vmax=1,k=0,ax=None,nlvl=11,avg=False,bar=True,cmap=None):
            """Produces a contour plot of the variable varname"""
            f,a,im=self.blk[0].contourf(varname,vmin,vmax,k,ax,nlvl,avg,bar,cmap)
            for i in range(self.nbk):
                self.blk[i].contourf(varname,vmin,vmax,k,a,nlvl,avg,bar,cmap)
            return f,a,im

        def setTouch(self,x0=0,y0=0,R=0.001,H=1,acc=1e-6):
            """docstring for setTouch"""
            for i in range(self.nbk):
                self.blk[i].setTouch(x0,y0,R,H,acc)

        def clone(self):
            """Returns a clone of the flow object"""
            obj=copy.copy(self)
            return obj
        
            
        
