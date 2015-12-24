import numpy as np
import copy as copy
import matplotlib.pyplot as plt

def naca(x,c,t0):
    """NACA 00XX series equation"""
    import numpy as np
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
    
def qhermite(x,y0,yx0,y1,yx1,yxx0=0.0,yxx1=0.0):
    """Quintic Hermitiean interpolation"""
    l=len(x)-1; s=np.zeros(l+1);f=np.zeros(l+1)
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
    else:
        c=np.matrix([
        [-06.0,+15.0,-10.0,+00.0,+00.0,+01.0],
        [-03.0,+08.0,-06.0,+00.0,+01.0,+00.0],
        [-00.5,+01.5,-01.5,+00.5,+00.0,+00.0],
        [+00.5,-01.0,+00.5,+00.0,+00.0,+00.0],
        [-03.0,+07.0,-04.0,+00.0,+00.0,+00.0],
        [+06.0,-15.0,+10.0,+00.0,+00.0,+00.0]])
                                                
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
    return f
    
def qahermite(x,y0,yx0,y1,yx1,side,lnr=0):
    """Quintic Hermitiean interpolation"""
    l=len(x)-1; s=np.zeros(l+1);f=np.zeros(l+1)
    x0=x[0];x1=x[l]; dx=x1-x0
    s[:]=(x[:]-x0)/dx
    p0=y0;p1=yx0*dx
    p3=y1;p2=yx1*dx;
    p=np.matrix([[p0,p1,p2,p3]])
    if (lnr==1):
        c=np.matrix([
        [-00.0,+00.0,-00.0,-01.0,+01.0],
        [-00.0,+00.0,-00.0,+00.0,+00.0],
        [-00.0,+00.0,-00.0,+00.0,+00.0],
        [+00.0,-00.0,+00.0,+01.0,+00.0]])
    else:
        if (side==1):
            c=np.matrix([
            [-03.0,+08.0,-06.0,+00.0,+01.0],
            [-01.0,+03.0,-03.0,+01.0,+00.0],
            [-02.0,+05.0,-03.0,+00.0,+00.0],
            [+03.0,-08.0,+06.0,+00.0,+00.0]])
        else:
            c=np.matrix([
            [+03.0,-04.0,-00.0,-00.0,+01.0],
            [+02.0,-03.0,+00.0,+01.0,+00.0],
            [+01.0,-01.0,+00.0,+00.0,+00.0],
            [-03.0,+04.0,+00.0,+00.0,+00.0]])
    for i in range(l+1):
        t=np.matrix([                                                                                                                                                       
                    [s[i]**4],
                    [s[i]**3],
                    [s[i]**2],
                    [s[i]**1],
                    [s[i]**0]
                            ])
        f[i]=(p*c*t)
    return f
       
def hermite(k1,k2,k3,k4,x0,x1,x):
    """Hermitiean interpolation"""
    l=x1-x0
    invl=1.0/l
    fx=(x-x0)*invl
    a=(2*fx**3-3*fx**2+1)
    b=l*(fx**3-2*fx**2+fx)
    c=(-2*fx**3+3*fx**2)
    d=l*(fx**3-fx**2)
    y=a*k1+b*k2+c*k3+d*k4
    return y

def dxi(x,fctr=1):
    n=max(x.shape); l=n-1; dx=np.zeros(n)
    a=np.zeros(n-2);b=np.zeros(n-1);c=np.zeros(n)
    d=np.zeros(n-1);e=np.zeros(n-2);
    a[-1]=1;e[0]=-1
    b[0:-1]=-1;b[-1]=-4;d[1:]=1;d[0]=4;
    c[0]=-3;c[-1]=3
    m=(1/(2*fctr))*np.mat(np.diag(a,-2)+np.diag(b,-1)+np.diag(c,0)+np.diag(d,1)+np.diag(e,2))
    dx=m*np.mat(x).T
    dx=np.asarray(dx.T)[0]
    return dx

class corner():
    def __init__(self,x=0.0,y=0.0,
                 dydx=0.0,dxdxi=0.0,d2xdxi2=0.0,
                 dxdy=0.0,dydet=0.0,d2ydet2=0.0):
        self.x=x
        self.y=y
    def clone (self):
        obj=copy.copy(self)
        return obj

class line():
    def __init__(self,n,c0,c1,dr,naca=0.0):
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
            self.c1.dydx=self.c0.dydx
            self.c1.dxdxi=self.c0.dxdxi
            self.c1.d2xdxi2=self.c0.d2xdxi2
        if (self.dr==1):
            self.c0.dxdy=(self.c1.x-self.c0.x)/(self.c1.y-self.c0.y)
            self.c0.dydet=(self.c1.y-self.c0.y)/(n-1)
            self.c0.d2ydet2=0.0
            self.c1.dxdy=self.c0.dxdy
            self.c1.dydet=self.c0.dydet
            self.c1.d2ydet2=self.c0.d2ydet2
        self.getPoints()

    def clone (self):
        obj=copy.copy(self)
        return obj

    def getPoints(self,plot=False):
        if (self.dr==0):
            if (self.naca==0.0):
               self.x=qhermite(np.arange(self.n),
                               self.c0.x,self.c0.dxdxi,
                               self.c1.x,self.c1.dxdxi,
                               self.c0.d2xdxi2,self.c1.d2xdxi2)
               self.dxdxi=dxi(self.x)
               self.d2xdxi2=dxi(self.dxdxi)
               self.y=qhermite(self.x,
                               self.c0.y,self.c0.dydx,
                               self.c1.y,self.c1.dydx)
            else:
               c=self.c1.x-self.c0.x
               self.x=np.linspace(self.c0.x,self.c1.x,self.n)
               self.y=np.zeros(self.n)
               for i in range(8):
                  err=2;tol=1.0+1e-4;l=8
                  dx=self.c0.dxdxi
                  while (abs(err)>tol):
                     x=self.x[i]+dx
                     y=naca(x,c,self.naca)
                     dy=y-self.y[i]
                     ds=np.sqrt(dx**2+dy**2)
                     err=ds/self.c0.dxdxi
                     dx=dx/err
                  self.x[i+1]=x
                  self.y[i+1]=y
               dxdxi=np.asarray(np.mat(self.x[l-5:l])*((1.0/12)*np.matrix([[3,-16,36,-48,25]]).T))[0]
               self.x[l-1:]=qhermite(np.arange(l-1,self.n),
                               self.x[l-1],dxdxi,
                               self.c1.x,self.c1.dxdxi,
                               0.0,self.c1.d2xdxi2)
               self.dxdxi=dxi(self.x)
               self.d2xdxi2=dxi(self.dxdxi)
               self.y=naca(self.x,c,self.naca)
               self.c0.dydx=(self.y[1]-self.y[0])/(self.x[1]-self.x[0])
               self.c1.dydx=(self.y[-1]-self.y[-2])/(self.x[-1]-self.x[-2])
            self.dydxi=dxi(self.y)
            self.d2ydxi2=dxi(self.dydxi)
        if (self.dr==1):
            self.y=qhermite(np.arange(self.n),
                            self.c0.y,self.c0.dydet,
                            self.c1.y,self.c1.dydet,
                            self.c0.d2ydet2,self.c1.d2ydet2)
            self.dydet=dxi(self.y)
            self.d2ydet2=dxi(self.dydet)
            self.x=qhermite(self.y,
                            self.c0.x,self.c0.dxdy,
                            self.c1.x,self.c1.dxdy)
            self.dxdet=dxi(self.x)
            self.d2xdet2=dxi(self.dxdet)
        if (plot==True):
            self.plot()

    def plot(self):
        plt.subplot(411)
        plt.plot(self.x,self.y)
        if (self.dr==0):
            plt.subplot(412)
            plt.plot(self.x)
            plt.subplot(413)
            plt.plot(self.dxdxi)
            plt.subplot(414)
            plt.plot(self.d2xdxi2)
            plt.show()
        if (self.dr==1):
            plt.subplot(412)
            plt.plot(self.y)
            plt.subplot(413)
            plt.plot(self.dydet)
            plt.subplot(414)
            plt.plot(self.d2ydet2)
            plt.show()
    def draw(self,show=False):
        plt.plot(self.x,self.y,'b')
        plt.plot(self.x,self.y,'ro')
        plt.draw()
        if (show==True):
           plt.show()

class block():
    def __init__(self,bs,bn,bw,be):
        self.bs=bs
        self.bn=bn
        self.bw=bw
        self.be=be
        self.lxi=bs.n
        self.let=bw.n

        self.bs.dxdet=qhermite(np.arange(self.lxi),
                               self.bw.dxdet[0],0.0,
                               self.be.dxdet[0],0.0)
        self.bn.dxdet=qhermite(np.arange(self.lxi),
                               self.bw.dxdet[-1],0.0,
                               self.be.dxdet[-1],0.0)
        self.bs.dydet=qhermite(np.arange(self.lxi),
                               self.bw.dydet[0],0.0,
                               self.be.dydet[0],0.0)
        self.bn.dydet=qhermite(np.arange(self.lxi),
                               self.bw.dydet[-1],0.0,
                               self.be.dydet[-1],0.0)
        self.bw.dydxi=qhermite(np.arange(self.let),
                               self.bs.dydxi[0],0.0,
                               self.bn.dydxi[0],0.0)
        self.be.dydxi=qhermite(np.arange(self.let),
                               self.bs.dydxi[-1],0.0,
                               self.bn.dydxi[-1],0.0)
        self.bw.dxdxi=qhermite(np.arange(self.let),
                               self.bs.dxdxi[0],0.0,
                               self.bn.dxdxi[0],0.0)
        self.be.dxdxi=qhermite(np.arange(self.let),
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

    def draw(self,drw=1,shw=1,stl='b'):
        """docstring for plotgrid"""
        nxi=self.lxi;net=self.let
        xis=0;xie=nxi-1;ets=0;ete=net-1
        plt.plot(self.x[:,ets],self.y[:,ets],'ro',markersize=4)
        plt.plot(self.x[:,ete],self.y[:,ete],'ro',markersize=4)
        plt.plot(self.x[xis,:],self.y[xis,:],'ro',markersize=4)
        plt.plot(self.x[xie,:],self.y[xie,:],'ro',markersize=4)
        for j in range(0,net):
            plt.plot(self.x[:,j],self.y[:,j],stl)
        for i in range(0,nxi):
            plt.plot(self.x[i,:],self.y[i,:],stl)
    
        plt.axis('equal')
        if drw==1:
           plt.draw()
        if shw==1:
           plt.show()
        pass

    def clone (self):
        obj=copy.copy(self)
        return obj


        



















