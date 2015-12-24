import matplotlib.pyplot as plt
import numpy as np
import lib.indx as indx

class plot2d():
	"""docstring for plot2d"""
	def __init__(self,fname):
		"""Initialise class with file name and axis titles.
		It also gets the data"""
		self.fname=fname
		self.setAxis('x','y')
		self.getData()

	def setFile(self,fname):
		"""Set the file name where the data is get and get it"""
		self.fname=fname
		self.getData()

	def setAxis(self,xname,yname):
		"""Set axis titles"""
		self.xname=xname
		self.yname=yname

	def getData(self):
		"""Get data from the file, just the 2 first columns"""
		self.x,self.y = np.loadtxt(self.fname,skiprows=1,unpack=True)

	def plot(self,k=0, lxi=150,grid=True,legend=True):
            """Plot data at a given plane"""
            self.lxi=lxi
            ls=indx.wall(0,k,lxi)
            le=indx.wall(lxi,k,lxi)
            self.f=plt.figure()
            self.a=self.f.add_subplot(111)
            self.a.plot(self.x[ls:le],self.y[ls:le],label=str(k))
            if legend is True:
            	plt.legend()
            self.a.grid(grid)
            self.a.set_xlabel(self.xname)
            self.a.set_ylabel(self.yname)

	def addPlot(self,k=0,lxi=150,grid=True,legend=True):
            """Add plot lines at given planes"""
            self.lxi=lxi
            ls=indx.wall(0,k,lxi)
            le=indx.wall(lxi,k,lxi)
            plt.plot(self.x[ls:le],self.y[ls:le],label=str(k))
            if legend is True:
            	plt.legend()
            plt.grid(grid)
            plt.xlabel(self.xname)
            plt.ylabel(self.yname)

	def showPlot(self,grid=True):
            """Draw the plot"""
            plt.xlabel(self.xname)
            plt.ylabel(self.yname)
            plt.grid(grid)
            plt.show()

class plotwall(plot2d):
	"""Class for plotting values of wall distances"""
	def __init__(self,fname='wplus.dat'):
		"""Initialise class from the first 2 columns"""
		self.fname=fname
		self.getData(0,1)

	def setFile(self,fname):
		"""Get data from file from first 2 columns"""
		self.fname=fname
		self.getData(0,1)

	def setAxis(self,xname,yname):
		"""Set axis titles customized for plotwall class"""
		self.xname=xname+'/c'
		self.yname=yname+'+'

	def getData(self,colx,coly):
		"""Get data from file at files colx and coly.
		x=>0, y=>1, z=>2"""
		if colx==0:
			self.xname='x/c'
		if colx==1:
			self.xname='y/c'
		if colx==2:
			self.xname='z/c'
		if coly==0:
		        self.yname='x+'
		if coly==1:
		        self.yname='y+'
		if coly==2:
			self.yname='z+'
		coly=coly+3
		self.x,self.y = np.loadtxt(self.fname,skiprows=1,unpack=True,usecols=(colx,coly))


class plotAvg(plot2d):
	"""Class for plotting average data. There are only two columns 
	at one fixed k plane"""
	def __init__(self, fname='avg.dat'):
		self.fname=fname
		self.getData()
		self.xname='x/c'
		self.yname='f(x)'

	def setAxis(self,xname,yname):
		"""Set axis titles customized for plotAvg class"""
		self.xname=xname+'/c'
		self.yname=yname

	def plot(self,grid=True):
            """Plot data"""
            plt.plot(self.x,self.y)
            plt.xlabel(self.xname)
            plt.ylabel(self.yname)
            plt.grid(grid)
            plt.show()
