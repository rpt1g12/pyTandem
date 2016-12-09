# -*- coding: utf-8 -*-
"""
Created on Tue Dec 6 09:00:00 2016

@author: rpt1g12

This module contains the classes needed for multi-block
structured grid generation
"""

# Necessary imports
from lib.myPlots import * 
import numpy as np
import copy as copy
import matplotlib.pyplot as plt
import os
from matplotlib.lines import Line2D
from libGrid.gridFunctions import *

#=========CLASSES================    

class point(object):
    """Defines the point class
        Parameters:
            x (real): x-coordinate
            y (real): y-coordinate
            dxy (real[2]): outward vertical direction [south,north] (dx/dy)
            dyx (real[2]): outward horizontal direction [west,east] (dydx)
            dx (real): horizontal spacing (dx/dxi)
            dy (real): vertical spacing (dy/det)
            dxx (real): horizontal spacing derivative (dx2/dxi2)
            dyy (real): horizontal spacing derivative (dy2/det2)"""
    def __init__(self,x=0.0,y=0.0,
                     dxy=[None,None],dyx=[None,None],
                     dx=None,dy=None,
                     dxx=None,dyy=None):
        """Point constructor"""
        super(point, self).__init__()
        self.x = x
        self.y = y
        self.dxy = dxy
        self.dyx = dyx
        self.dx = dx
        self.dy = dy
        self.dxx = dxx
        self.dyy = dyy
        pass
        

class line(Line2D):
    """
    Defines the line class

    Parameters:
        n (int): Total number of points
        dr (int): Direction. 0 is horizontal (west to east) and 1 is vertical (south to north)
        bc (int): Boundary condition:
            10: Inflow/Outflow
            20: Slip wall
            25: No-Slip wall
            30: Characteristic Interface
            35: Interface
            45: Periodic
        nb (int): Number of segments that define th piece-wise line
        nseg (list(real)[nb]): number of points per segment
        pts (point[1+nb]): Points that the define the piece-wise line
        shape (real): Shape of the line. If shape<=0 normal line, if shape>0 uses NACA00XX series function
        x (np.array(real)[n]): x-coordinates
        dx (np.array(real)[n]): x-coordinate derivatives (dx/dxi)
        dxx (np.array(real)[n]): x-coordinate derivatives (dx2/dxi2)
        dxdy (np.array(real)[n]): x-coordinate derivatives (dx/deta)
        dxdyy (np.array(real)[n]): x-coordinate derivatives (dx/deta2)
        y (np.array(real)[n]): y-coordinates
        dy (np.array(real)[n]): y-coordinate derivatives (dx/deta)
        dyy (np.array(real)[n]): y-coordinate derivatives (dx2/deta2)
        dydx (np.array(real)[n]): y-coordinate derivatives (dx/dxi)
        dydxx (np.array(real)[n]): y-coordinate derivatives (dx/dxi2)

    """
    def __init__(self,n,pts,io,dr,bc=30,shape=0.0):
        super(line,self).__init__(np.zeros(n),np.zeros(n))
        self.n=n
        self.x=np.zeros(n);self.dx=np.zeros(n);self.dxx=np.zeros(n);
        self.dxdy=np.zeros(n);self.dxdyy=np.zeros(n)
        self.y=np.zeros(n);self.dy=np.zeros(n);self.dyy=np.zeros(n);
        self.dydx=np.zeros(n);self.dydxx=np.zeros(n)
        self.dr=dr
        self.bc=bc
        if (((self.bc-20)*(self.bc-25))==0):
            self.shape=0.0
        else:
            self.shape=shape
        self.pts=pts
        self.nb=len(self.pts)-1
        if (len(io)+1)!=self.nb:
            print('Error! Number of segments does not agree with number of points')
            return
        self.io=[0]+io+[self.n-1]
        self.nseg=[]
        for i in range(self.nb):
            self.nseg.append(self.io[i+1]-self.io[i]+1)

        if (self.dr==0):
                if  self.pts[0].dyx[1]==None:
                    self.pts[0].dyx[1]=(self.pts[1].y-self.pts[0].y)/(self.pts[1].x-self.pts[0].x)
                for i in range(1,self.nb):
                    if  self.pts[i].dyx[0]==None:
                        self.pts[i].dyx[0]=(self.pts[i+1].y-self.pts[i-1].y)/(self.pts[i+1].x-self.pts[i-1].x)
                        self.pts[i].dyx[1]=self.pts[i].dyx[0]
                if  self.pts[-1].dyx[0]==None:
                    self.pts[-1].dyx[0]=self.pts[0].dyx[1]
                for i in [0,-1]:
                    if  self.pts[i].dx==None:
                        self.pts[i].dx=-1
                    if  self.pts[i].dxx==None:
                        self.pts[i].dxx=0.0
                for i in range(1,self.nb):
                    if  self.pts[i].dxx==None:
                        self.pts[i].dxx=0.0

        if (self.dr==1):
                if  self.pts[0].dxy[1]==None:
                    self.pts[0].dxy[1]=(self.pts[-1].x-self.pts[0].x)/(self.pts[-1].y-self.pts[0].y)
                if  self.pts[0].dy==None:
                    self.pts[0].dy=-1
                if  self.pts[0].dyy==None:
                    self.pts[0].dyy=0.0
                if  self.pts[-1].dxy[0]==None:
                    self.pts[-1].dxy[0]=self.pts[0].dxy[1]
                if  self.pts[-1].dy==None:
                    self.pts[-1].dy=self.pts[0].dy*1
                if  self.pts[-1].dyy==None:
                    self.pts[-1].dyy=self.pts[0].dyy*1

        self.getPoints()


    def getPoints(self,plot=False):
        if (self.dr==0):
            if (self.shape<=0.0):
                for i in range(self.nb):
                    io=self.io[i];nseg=self.nseg[i];ie=io+nseg
                    dxo=self.pts[i].dx;dxe=self.pts[i+1].dx
                    self.x[io:ie],self.dx[io:ie],self.dxx[io:ie]=gridf(self.pts[i].x,self.pts[i+1].x,
                                            self.pts[i].dx,dxe,
                                            self.nseg[i])

                    self.y[io:ie],tmp,tmp=qhermite(self.x[io:ie],
                                           self.pts[i].y,self.pts[i].dyx[1],
                                           self.pts[i+1].y,self.pts[i+1].dyx[0]) 
                    self.dydx[io:ie]=dxi(self.y[io:ie])
                    self.dydxx[io:ie]=dxi(self.dydx[io:ie])
            else:
               c=self.pts[1].x-self.pts[0].x
               self.x=np.linspace(self.pts[0].x,self.pts[1].x,self.n)
               self.dx=self.x.copy()
               self.dxx=self.x.copy()
               self.y=np.zeros(self.n)

               tol=1.0+1e-4;l=4
               for i in range(l):
                  err=2;
                  dx=self.pts[0].dx
                  while (abs(err)>tol):
                     x=self.x[i]+dx
                     y=naca(x-self.x[0],c,self.shape)
                     dy=y-self.y[i]
                     ds=np.sqrt(dx**2+dy**2)
                     err=ds/self.pts[0].dx
                     dx=dx/err
                  self.x[i+1]=x
                  self.dx[i]=(self.x[i+1]-self.x[i])/2
                  self.y[i+1]=y
               xo=self.x[l];xn=self.pts[1].x
               dxo=np.sum(self.x[l-4:l+1]*np.array([3,-16,36,-48,25]))/12.0
               dxn=self.pts[1].dx
               n=self.n-l
               self.dx[:l+1]=dxi(self.x[:l+1])
               self.dxx[:l+1]=dxi(self.dx[:l+1])
               self.x[l:],self.dx[l:],self.dxx[l:]=gridf(xo,xn,dxo,dxn,n)
               self.y=naca(self.x-self.x[0],c,self.shape)
               self.pts[0].dyx[1]=(self.y[1]-self.y[0])/(self.x[1]-self.x[0])
               self.pts[1].dyx[0]=(self.y[-1]-self.y[-2])/(self.x[-1]-self.x[-2])
            self.dydx=dxi(self.y)
            self.dydxx=dxi(self.dydx)
        if (self.dr==1):
            self.y,self.dy,self.dyy=gridf(self.pts[i].y,self.pts[i+1].y,
                                    self.pts[i].dy,self.pts[i+1].dy,
                                    self.nseg[i])
            self.x,tmp,tmp=qhermite(self.y,
                            self.pts[i].x,self.pts[i+1].dxy[1],
                            self.pts[i+1].x,self.pts[i+1].dxy[0])
            self.dxdy=dxi(self.x)
            self.d2xdyy=dxi(self.dxdy)
        if (plot==True):
            self.plot()
        self.set_xdata(self.x)
        self.set_ydata(self.y)

    def linePlot(self,fig=None):
        self.getPoints()
        if (fig==None): 
            fig,ax=getFig('Line Properties',411)
        else:
            fig.clear()
            fig.canvas.draw()
            fig.canvas.set_window_title('Line Properties')
            fig.set_tight_layout('tight')
            ax=addAx(fig,411)
        ax.plot(self.x,self.y,'b-o')
        fit(ax)
        if (self.dr==0):
            xxi=dxi(self.x)
            xxi2=dxi(self.dx)
            ax2=addAx(fig,412)
            ax2.plot(self.x,'b-o')
            fit(ax2)
            ax3=addAx(fig,413)
            ax3.plot(self.dx,'b-o')
            ax3.plot(xxi,'r-x')
            fit(ax3)
            ax4=addAx(fig,414)
            ax4.plot(self.dxx,'b-o')
            ax4.plot(xxi2,'r-x')
            fit(ax4)
            fig.canvas.draw()
        if (self.dr==1):
            yet=dxi(self.y)
            yet2=dxi(self.dy)
            ax2=addAx(fig,412)
            ax2.plot(self.y,'b-o')
            fit(ax2)
            ax3=addAx(fig,413)
            ax3.plot(self.dy,'b-o')
            ax3.plot(yet,'r-x')
            fit(ax3)
            ax4=addAx(fig,414)
            ax4.plot(self.dyy,'b-o')
            ax4.plot(yet2,'r-x')
            fit(ax4)
            fig.canvas.draw()
        fig.show()

    def clone (self):
        obj=copy.copy(self)
        return obj

class block():
    def __init__(self,lines):
        self.bs=lines[0]
        self.bn=lines[1]
        self.bw=lines[2]
        self.be=lines[3]
        self.lxi=lines[0].n
        self.let=lines[2].n
        self.lze=1
        #self.update(0)
    
    def updateLines(self):
        self.bs.getPoints()
        self.bn.getPoints()
        self.bw.getPoints()
        self.be.getPoints()

    def update(self,opt=1):
        if (opt==1):
            self.updateLines()

        self.bs.dxdet,tmp,tmp=qhermite(self.bs.x,
                             self.bw.dxdet[0],0.0,
                             self.be.dxdet[0],0.0)
        self.bn.dxdet,tmp,tmp=qhermite(self.bn.x,
                             self.bw.dxdet[-1],0.0,
                             self.be.dxdet[-1],0.0)
        self.bs.dydet,tmp,tmp=qhermite(self.bs.x,
                             self.bw.dydet[0],0.0,
                             self.be.dydet[0],0.0)
        self.bn.dydet,tmp,tmp=qhermite(self.bn.x,
                             self.bw.dydet[-1],0.0,
                             self.be.dydet[-1],0.0)
        self.bw.dydxi,tmp,tmp=qhermite(self.bw.y,
                             self.bs.dydxi[0],0.0,
                             self.bn.dydxi[0],0.0)
        self.be.dydxi,tmp,tmp=qhermite(self.be.y,
                             self.bs.dydxi[-1],0.0,
                             self.bn.dydxi[-1],0.0)
        self.bw.dxdxi,tmp,tmp=qhermite(self.bw.y,
                             self.bs.dxdxi[0],0.0,
                             self.bn.dxdxi[0],0.0)
        self.be.dxdxi,tmp,tmp=qhermite(self.be.y,
                               self.bs.dxdxi[-1],0.0,
                               self.bn.dxdxi[-1],0.0)
        self.x=np.zeros((self.lxi,self.let))
        self.y=np.zeros((self.lxi,self.let))
        self.x[:,0]=self.bs.x
        self.x[:,-1]=self.bn.x
        self.y[:,0]=self.bs.y
        self.y[:,-1]=self.bn.y
        self.x[0,:]=self.bw.x
        self.x[-1,:]=self.be.x
        self.y[0,:]=self.bw.y
        self.y[-1,:]=self.be.y
        self.getPoints()

    def getPoints(self,opt=0):
        nxi=self.lxi;net=self.let
        xis=0;xie=nxi-1;ets=0;ete=net-1
        xv=np.zeros((nxi,2));yv=np.zeros((nxi,2))
        xu=np.zeros((net,2));yu=np.zeros((net,2))
        xuv=np.zeros((2,2)); yuv=np.zeros((2,2))
        xv[:,0]=self.bs.dxdet;xv[:,1]=self.bn.dxdet
        xu[:,0]=self.bw.dxdxi;xu[:,1]=self.be.dxdxi
        yv[:,0]=self.bs.dydet;yv[:,1]=self.bn.dydet
        yu[:,0]=self.bw.dydxi;yu[:,1]=self.be.dydxi
        
        tmp=dxi(xv[:,0])
        xuv[0,0]=tmp[0];xuv[1,0]=tmp[-1]
        tmp=dxi(xv[:,1])
        xuv[0,1]=tmp[0];xuv[1,1]=tmp[-1]
        tmp=dxi(yv[:,0])
        yuv[0,0]=tmp[0];yuv[1,0]=tmp[-1]
        tmp=dxi(yv[:,1])
        yuv[0,1]=tmp[0];yuv[1,1]=tmp[-1]

        
        Mx=np.matrix([
        [self.x[xis,ets],self.x[xis,ete],xv[xis,0],xv [xis,1]],
        [self.x[xie,ets],self.x[xie,ete],xv[xie,0],xv [xie,1]],
        [xu[ets,0 ],xu[ete,0 ],xuv[0 ,0],xuv[0  ,1]],
        [xu[ets,1 ],xu[ete,1 ],xuv[1 ,0],xuv[1  ,1]]])

        My=np.matrix([
        [self.y[xis,ets],self.y[xis,ete],yv[xis,0],yv [xis,1]],
        [self.y[xie,ets],self.y[xie,ete],yv[xie,0],yv [xie,1]],
        [yu[ets,0 ],yu[ete,0 ],yuv[0 ,0],yuv[0  ,1]],
        [yu[ets,1 ],yu[ete,1 ],yuv[1 ,0],yuv[1  ,1]]])

        if  opt==0:
            a0 = lambda u: 2*u**3-3*u**2+1
            a1 = lambda u: 3*u**2-2*u**3
            b0 = lambda u: u**3-2*u**2+u
            b1 = lambda u: u**3-u**2
        elif  opt==1:
            a0 = lambda u: (1-u)
            a1 = lambda u: u
            b0 = lambda u: 0.0
            b1 = lambda u: 0.0

        for j in range(0,net):
            v=(j-ets)/(ete-ets)
            abv=np.matrix([[a0(v),a1(v),b0(v),b1(v)]])
            Xv=np.matrix([[self.x [xis, j ]],
                          [self.x [xie, j ]],
                          [xu[ j , 0 ]],
                          [xu[ j , 1 ]]])
            Yv=np.matrix([[self.y [xis, j ]],
                          [self.y [xie, j ]],
                          [yu[ j , 0 ]],
                          [yu[ j , 1 ]]])
            abvt=abv.T
            for i in range(0,nxi):
                u=(i-xis)/(xie-xis)
                abu=np.matrix([[a0(u),a1(u),b0(u),b1(u)]])
                Xu=np.matrix([[self.x [ i ,ets]],
                              [self.x [ i ,ete]],
                              [xv[ i , 0 ]],
                              [xv[ i , 1 ]]])
                Yu=np.matrix([[self.y [ i ,ets]],
                              [self.y [ i ,ete]],
                              [yv[ i , 0 ]],
                              [yv[ i , 1 ]]])


                self.x[i,j]=abu*Xv+abv*Xu-(abu*Mx*abvt)
                self.y[i,j]=abu*Yv+abv*Yu-(abu*My*abvt)
        self.z=np.zeros((self.lxi,self.let))

    def draw(self,ax=None,drw=1,shw=1,stl='b',update=False):
        """docstring for plotgrid"""
        if (ax==None):
            ax=plt.gca()
        if (update):
            self.update()
        nxi=self.lxi;net=self.let
        xis=0;xie=nxi-1;ets=0;ete=net-1
        ax.plot(self.x[:,ets],self.y[:,ets],'ro',markersize=4)
        ax.plot(self.x[:,ete],self.y[:,ete],'ro',markersize=4)
        ax.plot(self.x[xis,:],self.y[xis,:],'ro',markersize=4)
        ax.plot(self.x[xie,:],self.y[xie,:],'ro',markersize=4)
        for j in range(0,net):
            ax.plot(self.x[:,j],self.y[:,j],stl)
        for i in range(0,nxi):
            ax.plot(self.x[i,:],self.y[i,:],stl)
    
        if drw==1:
           ax.figure.canvas.draw()
           ax.set_aspect('equal')
        if shw==1:
           ax.figure.show()
           fit(ax)
        pass

    def clone (self):
        obj=copy.copy(self)
        return obj

    def write(self,number=0):
        """Write to Plot3D format"""
        path=os.getcwd()
        self.number=number
        fh=open(path+'/block'+str(number)+'.xyz','wb')
        hdr=np.int32(np.array([1,self.lxi,self.let,self.lze]))
        fh.write(hdr)
        fh.write(np.float32(np.transpose(self.x).copy(order='C')))
        fh.write(np.float32(np.transpose(self.y).copy(order='C')))
        fh.write(np.float32(np.transpose(self.z).copy(order='C')))
        fh.close()
        print('Block'+str(number)+' written!')
        pass
