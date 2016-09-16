# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 12:13:46 2016

@author: rpt1g12

This module contains the funtions and classes needed for multi-block
structured grid generation
"""

# Necessary imports
from lib.myPlots import * 
import numpy as np
import copy as copy
import matplotlib.pyplot as plt
import os
from matplotlib.lines import Line2D

#=========FUNCTIONS================
def replot(axes=None,doFit=False):
    if (axes==None):
        axes=plt.gca()
    for l in axes.lines:
        l.getPoints()
    if (doFit):
        fit(axes)
    else:
        axes.figure.canvas.draw()
        
        
def naca(x,c,t0):
    """NACA 00XX series equation
    (IN)
    x=> horizontal coordinate
    c=> chord lenght
    t0=> aerofoil thickness in chord percentage (XX in NACA 00XX)
    (OUT)
    y=> vertical coordinate"""
    t=abs(t0)*0.01
    k=0.991148635
    if (t0<0.0):
       fctr=-1.0
    else:
       fctr=1.0
    xc=x/c
    y= fctr*(t*c*k/0.2) \
    * (0.298222773*np.sqrt(xc)-0.127125232*xc-0.357907806*xc**2 \
    +0.291984971*xc**3-0.105174606*xc**4)
    return y
    
    
def gridf(xo,xn,dxo,dxn,mxin):
    """Distributes points along a line:
    (IN)
    xo=> Initial point
    xn=> Final point
    dxo=> Initial grid size
    dxn=> Final grid size
    mxin=> Total number of points
    (OUT)
    x=> np array with points
    xxi=> grid spacings/derivatives along the line"""
    x=np.zeros(mxin);xxi=x.copy()
    mxin=mxin-1
    dxoo=dxo;dxnn=dxn
    if (dxo<0):
        dxoo=(xn-xo)/(mxin-1)
    if(dxn<0):
        dxnn=(xn-xo)/(mxin-1)
    aa=6.0*(xn-xo)-3*mxin*(dxoo+dxnn)
    bb=15*(xo-xn)+mxin*(8*dxoo+7*dxnn)
    cc=10*(xn-xo)-mxin*(6*dxoo+4*dxnn)
    dd=mxin*dxoo; fctr=1.0/mxin
    
    for i in range(mxin+1):
        xi=i*fctr
        x[i]=aa*xi**5+bb*xi**4+cc*xi**3+dd*xi+xo
        xxi[i]=fctr*(5*aa*xi**4+4*bb*xi**3+3*cc*xi**2+dd)
    
    return x,xxi

def dxi(x,fctr=1):
    """Derivative scheme 2nd Order"""
    n=max(x.shape); dx=np.zeros(n)
    a=np.zeros(n-2);b=np.zeros(n-1);c=np.zeros(n)
    d=np.zeros(n-1);e=np.zeros(n-2);
    a[-1]=1;e[0]=-1
    b[0:-1]=-1;b[-1]=-4;d[1:]=1;d[0]=4;
    c[0]=-3;c[-1]=3
    m=(1/(2*fctr))*np.mat(np.diag(a,-2)+np.diag(b,-1)+np.diag(c,0)+np.diag(d,1)+np.diag(e,2))
    dx=m*np.mat(x).T
    dx=np.asarray(dx.T)[0]
    return dx
    
def qhermite(x,y0,yx0,y1,yx1,yxx0=0.0,yxx1=0.0):
    """Quintic Hermitiean interpolation"""
    l=len(x)-1; s=np.zeros(l+1);f=np.zeros(l+1)
    fx=np.zeros(l+1);fxx=np.zeros(l+1)
    x0=x[0];x1=x[l]; dx=x1-x0
    dy=y1-y0
    s[:]=(x[:]-x0)/dx
    p0=y0;p1=yx0*dx;p2=yxx0*dx**2
    p5=y1;p4=yx1*dx;p3=yxx1*dx**2
    p=np.matrix([[p0,p1,p2,p3,p4,p5]])
    if  (yx0==yx1 and y1==(y0+yx0*dx)):
        c=np.matrix([
        [+00.0,-00.0,+00.0,-00.0,-01.0,+01.0],
        [+00.0,-00.0,+00.0,-00.0,+00.0,+00.0],
        [+00.0,-00.0,+00.0,-00.0,+00.0,+00.0],
        [+00.0,-00.0,+00.0,-00.0,+00.0,+00.0],
        [+00.0,-00.0,+00.0,-00.0,+00.0,+00.0],
        [-00.0,+00.0,-00.0,+00.0,+01.0,+00.0]])

        c1=np.matrix([
        [+00.0,-00.0,+00.0,-00.0,+00.0,-01.0],
        [+00.0,-00.0,+00.0,-00.0,+00.0,+00.0],
        [+00.0,-00.0,+00.0,-00.0,+00.0,+00.0],
        [+00.0,-00.0,+00.0,-00.0,+00.0,+00.0],
        [+00.0,-00.0,+00.0,-00.0,+00.0,+00.0],
        [-00.0,+00.0,-00.0,+00.0,+00.0,+01.0]])

        c2=np.matrix([
        [+00.0,-00.0,+00.0,-00.0,+00.0,-00.0],
        [+00.0,-00.0,+00.0,-00.0,+00.0,+00.0],
        [+00.0,-00.0,+00.0,-00.0,+00.0,+00.0],
        [+00.0,-00.0,+00.0,-00.0,+00.0,+00.0],
        [+00.0,-00.0,+00.0,-00.0,+00.0,+00.0],
        [-00.0,+00.0,-00.0,+00.0,+00.0,+00.0]])
    else:
        c=np.matrix([
        [-06.0,+15.0,-10.0,+00.0,+00.0,+01.0],
        [-03.0,+08.0,-06.0,+00.0,+01.0,+00.0],
        [-00.5,+01.5,-01.5,+00.5,+00.0,+00.0],
        [+00.5,-01.0,+00.5,+00.0,+00.0,+00.0],
        [-03.0,+07.0,-04.0,+00.0,+00.0,+00.0],
        [+06.0,-15.0,+10.0,+00.0,+00.0,+00.0]])
        
        c1=np.matrix([
        [+00.0,-30.0,+60.0,-30.0,+00.0,+00.0],
        [+00.0,-15.0,+32.0,-18.0,+00.0,+01.0],
        [+00.0,-02.5,+06.0,-04.5,+01.0,+00.0],
        [+00.0,+02.5,-04.0,+01.5,+00.0,+00.0],
        [+00.0,-15.0,+28.0,-12.0,+00.0,+00.0],
        [+00.0,+30.0,-60.0,+30.0,+00.0,+00.0]])
        
        c2=np.matrix([
        [+000.0,+000.0,-120.0,+180.0,-060.0,+000.0],
        [+000.0,+000.0,-060.0,+096.0,-036.0,+000.0],
        [+000.0,+000.0,-010.0,+018.0,-009.0,+001.0],
        [+000.0,+000.0,+010.0,-012.0,+003.0,+000.0],
        [+000.0,+000.0,-060.0,+084.0,-024.0,+000.0],
        [+000.0,+000.0,+120.0,-180.0,+060.0,+000.0]])
                                                
    for i in range(l+1):
        t=np.matrix([
                    [s[i]**5],
                    [s[i]**4],
                    [s[i]**3],
                    [s[i]**2],
                    [s[i]**1],
                    [s[i]**0]
                            ])
        f[i]=(p*c*t)
        fx[i]=(p*c1*t)/dx
        fxx[i]=(p*c2*t)/(dx**2)
    return f,fx,fxx


#=========CLASSES================    

class corner():
    def __init__(self,x=0.0,y=0.0,
                 dydx=0.0,dxdxi=0.0,d2xdxi2=0.0,
                 dxdy=0.0,dydet=0.0,d2ydet2=0.0):
        self.x=x
        self.y=y
        self.dydx=dydx
        self.dxdxi=dxdxi
        self.d2xdxi2=d2xdxi2
        self.dxdy=dxdy
        self.dydet=dydet
        self.d2ydet2=d2ydet2
        #self.modified=[False,False]
    def clone (self):
        obj=copy.copy(self)
        return obj

class line(Line2D):
    def __init__(self,n,c0,c1,dr,naca=0.0):
        super(line,self).__init__(np.zeros(n),np.zeros(n))
        self.n=n
        self.c0=c0
        self.c1=c1
        self.dr=dr
        self.naca=naca
        if (self.naca!=0):
           self.dr=0
        if (self.dr==0):
                self.c0.dydx=(self.c1.y-self.c0.y)/(self.c1.x-self.c0.x)
                self.c0.dxdxi=(self.c1.x-self.c0.x)/(n-1)
                self.c0.d2xdxi2=0.0
                self.c1.dydx=1*self.c0.dydx
                self.c1.dxdxi=1*self.c0.dxdxi
                self.c1.d2xdxi2=1*self.c0.d2xdxi2
        if (self.dr==1):
                self.c0.dxdy=(self.c1.x-self.c0.x)/(self.c1.y-self.c0.y)
                self.c0.dydet=(self.c1.y-self.c0.y)/(n-1)
                self.c0.d2ydet2=0.0
                self.c1.dxdy=1*self.c0.dxdy
                self.c1.dydet=1*self.c0.dydet
                self.c1.d2ydet2=1*self.c0.d2ydet2
        self.getPoints()

    def clone (self):
        obj=copy.copy(self)
        return obj

    def getPoints(self,plot=False):
        if (self.dr==0):
            if (self.naca==0.0):
                self.x,self.dxdxi,self.d2xdxi2=qhermite(np.arange(self.n),
                                                       self.c0.x,self.c0.dxdxi,
                                                       self.c1.x,self.c1.dxdxi,
                                                       self.c0.d2xdxi2,self.c1.d2xdxi2)
#                self.x,self.dxdxi=gridf(self.c0.x,self.c1.x,
#                                        self.c0.dxdxi,self.c1.dxdxi,
#                                        self.n)
#                self.dxdxi=dxi(self.x)
#                self.d2xdxi2=dxi(self.dxdxi)
                self.y,tmp,tmp=qhermite(self.x,
                                       self.c0.y,self.c0.dydx,
                                       self.c1.y,self.c1.dydx)
                self.dydxi=dxi(self.y)
                self.d2ydxi2=dxi(self.dydxi)
            else:
               c=self.c1.x-self.c0.x
               self.x=np.linspace(self.c0.x,self.c1.x,self.n)
               self.y=np.zeros(self.n)
               for i in range(8):
                  err=2;tol=1.0+1e-4;l=8
                  dx=self.c0.dxdxi
                  while (abs(err)>tol):
                     x=self.x[i]+dx
                     y=naca(x-self.x[0],c,self.naca)
                     dy=y-self.y[i]
                     ds=np.sqrt(dx**2+dy**2)
                     err=ds/self.c0.dxdxi
                     dx=dx/err
                  self.x[i+1]=x
                  self.y[i+1]=y
               dxdxi=np.asarray(np.mat(self.x[l-5:l])*((1.0/12)*np.matrix([[3,-16,36,-48,25]]).T))[0]
               self.x[l-1:],tmp,tmp=qhermite(np.arange(l-1,self.n),
                               self.x[l-1],dxdxi,
                               self.c1.x,self.c1.dxdxi,
                               0.0,self.c1.d2xdxi2)
               self.dxdxi=dxi(self.x)
               self.d2xdxi2=dxi(self.dxdxi)
               self.y=naca(self.x-self.x[0],c,self.naca)
               self.c0.dydx=(self.y[1]-self.y[0])/(self.x[1]-self.x[0])
               self.c1.dydx=(self.y[-1]-self.y[-2])/(self.x[-1]-self.x[-2])
            self.dydxi=dxi(self.y)
            self.d2ydxi2=dxi(self.dydxi)
        if (self.dr==1):
            self.y,self.dydet,self.d2ydet2=qhermite(np.arange(self.n),
                            self.c0.y,self.c0.dydet,
                            self.c1.y,self.c1.dydet,
                            self.c0.d2ydet2,self.c1.d2ydet2)
#            self.y,self.dydet=gridf(self.c0.y,self.c1.y,
#                                    self.c0.dydet,self.c1.dydet,
#                                    self.n)
#            self.dydet=dxi(self.y)
#            self.d2ydet2=dxi(self.dydet)
            self.x,tmp,tmp=qhermite(self.y,
                            self.c0.x,self.c0.dxdy,
                            self.c1.x,self.c1.dxdy)
            self.dxdet=dxi(self.x)
            self.d2xdet2=dxi(self.dxdet)
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
            xxi2=dxi(self.dxdxi)
            ax2=addAx(fig,412)
            ax2.plot(self.x,'b-o')
            fit(ax2)
            ax3=addAx(fig,413)
            ax3.plot(self.dxdxi,'b-o')
            ax3.plot(xxi,'r-x')
            fit(ax3)
            ax4=addAx(fig,414)
            ax4.plot(self.d2xdxi2,'b-o')
            ax4.plot(xxi2,'r-x')
            fit(ax4)
            fig.canvas.draw()
        if (self.dr==1):
            yet=dxi(self.y)
            yet2=dxi(self.dydet)
            ax2=addAx(fig,412)
            ax2.plot(self.y,'b-o')
            fit(ax2)
            ax3=addAx(fig,413)
            ax3.plot(self.dydet,'b-o')
            ax3.plot(yet,'r-x')
            fit(ax3)
            ax4=addAx(fig,414)
            ax4.plot(self.d2ydet2,'b-o')
            ax4.plot(yet2,'r-x')
            fit(ax4)
            fig.canvas.draw()
        fig.show()
        
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