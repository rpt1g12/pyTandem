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
#Define a float32
float32=np.float32
import lib.Plot3DClasses as p3d

class circle(p3d.flow):
    """Circle class
    Inheritates from:
        flow class"""
    def __init__(self,path,gfile,sfile,x0,y0,z0,r,theta0,theta1,thetas):
        super(circle, self).__init__(path,gfile,sfile)
        self.x0,y0,z0,r,theta0,theta1,thetas = x0,y0,z0,r,theta0,theta1,thetas
        x,y,z,theta=self.getCirc(x0,y0,z0,r,theta0,theta1,thetas)
        nbks=[1,1,1]
        lxib=[x.shape[0]]
        letb=[x.shape[1]]
        lzeb=[x.shape[2]]
        self.setHdr(nbks,lxib,letb,lzeb)
        self.blk[0].setData('x',x)
        self.blk[0].setData('y',y)
        self.blk[0].setData('z',z)
        self.blk[0].setData('th',theta)

    def getCirc(self,x0,y0,z0,r,theta0,theta1,thetas):
        """Computes the X and Y coorditates of a circle or radius r and centre at (x0,y0). 
        Only coordinates from angle theta0 to theta1 with spacing thetas are computed"""
        theta=np.linspace(theta0,theta1,thetas)
        ntheta=len(theta)
        x=np.zeros((ntheta,1,1))
        y=np.zeros((ntheta,1,1))
        z=np.zeros((ntheta,1,1))
        for n in range(ntheta):
            x[n,0,0]=(r*np.cos(theta[n]))+x0
            y[n,0,0]=(r*np.sin(theta[n]))+y0
        z[:,:,:]=z0
        return x,y,z,theta

class line(p3d.flow):
    """Line class"""
    def __init__(self, path,gfile,sfile,p0,p1,n):
        super(line, self).__init__(path,gfile,sfile)
        p0=np.asarray(p0)
        p1=np.asarray(p1)
        self.p0,p1,n = p0,p1,n
        x,y,z,l=self.getLine(p0,p1,n)
        nbks=[1,1,1]
        lxib=[x.shape[0]]
        letb=[x.shape[1]]
        lzeb=[x.shape[2]]
        self.setHdr(nbks,lxib,letb,lzeb)
        self.blk[0].setData('x',x)
        self.blk[0].setData('y',y)
        self.blk[0].setData('z',z)
        self.blk[0].setData('l',l)
    
    def getLine(self,p0,p1,n):
        """Computes the X and Y coordinates of a line in between points p0 and p1. The number of
        points is n
        
        Args:
            p0 (list[3]): x,y and z coordinates of p0
            p1 (list[3]): x,y and z coordinates of p1
            n (int): Number of points
        
        Returns:
            x (np.array[n,1,1]): x coordinates
            y (np.array[n,1,1]): y coordinates
            z (np.array[n,1,1]): z coordinates
            """
        r=p1-p0
        dr=r/(n-1)
        dr_mag=np.linalg.norm(dr)
        x=np.zeros((n,1,1));y=np.zeros((n,1,1));z=np.zeros((n,1,1));l=np.zeros((n,1,1))
        for i in range(n):
            x[i,0,0]=p0[0]+i*dr[0]
            y[i,0,0]=p0[1]+i*dr[1]
            z[i,0,0]=p0[2]+i*dr[2]
            l[i,0,0]=i*dr_mag
        return x,y,z,l
        
        

