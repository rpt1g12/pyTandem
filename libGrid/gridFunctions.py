# -*- coding: utf-8 -*-
"""
Created on Tue Dec 6 09:00:00 2016

@author: rpt1g12

This module contains the funtions needed for multi-block
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

    Arguments:
        xo (real): Initial point
        xn (real): Final point
        dxo (real): Initial grid size. If dx0<0 free edge
        dxn (real): Final grid size. If dx0<0 free edge
        mxin (int): Total number of points

    Returns:
        x (np.array(real)[mxin]): np array with points
        xxi(np.array(real)[mxin]): grid spacings/derivatives along the line
    """

    x=np.zeros(mxin);dx=x.copy();dxx=dx.copy()
    mxin=mxin-1
    dxoo=dxo;dxnn=dxn
    if (dxo<0):
        dxoo=(xn-xo)/(mxin)
    if(dxn<0):
        dxnn=(xn-xo)/(mxin)
    aa=6.0*(xn-xo)-3*mxin*(dxoo+dxnn)
    bb=15*(xo-xn)+mxin*(8*dxoo+7*dxnn)
    cc=10*(xn-xo)-mxin*(6*dxoo+4*dxnn)
    dd=mxin*dxoo; fctr=1.0/mxin
    
    for i in range(mxin+1):
        xi=i*fctr
        x[i]=aa*xi**5+bb*xi**4+cc*xi**3+dd*xi+xo
        dx[i]=fctr*(5*aa*xi**4+4*bb*xi**3+3*cc*xi**2+dd)
        dxx[i]=fctr*fctr*(20*aa*xi**3+12*bb*xi**2+6*cc*xi)
    
    return x,dx,dxx

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
