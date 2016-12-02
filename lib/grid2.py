from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np

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

def dxi(x,fctr=1):
    """Derivative scheme 2nd Order"""
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
    def __init__(self,x,y):
        self.x=x
        self.y=y

#class line(Line2D):
#    def __init__(self,c0,c1,n=2):
#        self.n=max(n,2)
#        self.c0=c0
#        self.c1=c1
#        self.x=np.zeros(n)
#        self.y=np.zeros(n)
#        self.getPoints()
#        super(line,self).__init__(self.x,self.y)
#
#    def getPoints(self):
#        self.v=((self.c1.x-self.c0.x),(self.c1.y-self.c0.y))
#        self.l=np.sqrt(self.v[0]**2+self.v[1]**2)
#        self.duv=self.v/self.l
#        self.dnv=self.duv*(self.l/(self.n-1))
#        for i in range(self.n):
#            self.x[i]=self.c0.x+i*self.dnv[0]
#            self.y[i]=self.c0.y+i*self.dnv[1]

class line(Line2D):
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
        super(line,self).__init__(self.x,self.y)

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
                     y=naca(x-self.x[0],c,self.naca)
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
               self.y=naca(self.x-self.x[0],c,self.naca)
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

    def linePlot(self):
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
    def lineDraw(self,show=False):
        plt.plot(self.x,self.y,'b')
        plt.plot(self.x,self.y,'ro')
        plt.draw()
        if (show==True):
           plt.show()
