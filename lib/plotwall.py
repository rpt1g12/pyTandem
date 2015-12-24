import matplotlib.pyplot as plt
import numpy as np
import indx 

class plotwall:
	"""docstring for plotwall"""
	def __init__(self,fname='wplus.dat'):
		"""docstring for __init__"""
		self.fname=fname
		self.getData(0,1)

	def setFile(self,fname):
		"""docstring for setFile"""
		self.fname=fname
		self.getData(0,1)

	def setAxis(self,xname,yname):
		"""docstring for setAxis"""
		self.xname=xname+'/c'
		self.yname=yname+'+'

	def getData(self,colx,coly):
		"""docstring for getPlot"""
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

	def plot(self,k):
		"""docstring for plot"""
        	ls=indx.wall(0,k,150)
        	le=indx.wall(150,k,150)
        	plt.plot(self.x[ls:le],self.y[ls:le])
        	plt.xlabel(self.xname)
        	plt.ylabel(self.yname)
		plt.grid(True)
        	plt.show()

	def addPlot(self,k):
		"""docstring for plot"""
        	ls=indx.wall(0,k,150)
        	le=indx.wall(150,k,150)
        	plt.plot(self.x[ls:le],self.y[ls:le])

	def showPlot(self):
		"""docstring for plot"""
        	plt.xlabel(self.xname)
        	plt.ylabel(self.yname)
		plt.grid(True)
        	plt.show()

class plotAvg():
	"""docstring for plotAvg"""
	def __init__(self, fname='avg.dat'):
		self.fname=fname
		self.getData
		sefl.xname='x/c'
		sefl.yname='f(x)'
		
	def setFile(self,fname):
		"""docstring for setFile"""
		self.fname=fname
		self.getData()

	def setAxis(self,xname,yname):
		"""docstring for setAxis"""
		self.xname=xname+'/c'
		self.yname=yname

	def getData(self,colx=0,coly=1):
		"""docstring for getPlot"""
		self.x,self.y = np.loadtxt(self.fname,skiprows=1,unpack=True,usecols=(colx,coly))

	def plot(self,k):
		"""docstring for plot"""
        	ls=indx.wall(0,k,150)
        	le=indx.wall(150,k,150)
        	plt.plot(self.x[ls:le],self.y[ls:le])
        	plt.xlabel(self.xname)
        	plt.ylabel(self.yname)
		plt.grid(True)
        	plt.show()

	def addPlot(self,k):
		"""docstring for plot"""
        	ls=indx.wall(0,k,150)
        	le=indx.wall(150,k,150)
        	plt.plot(self.x[ls:le],self.y[ls:le])

	def showPlot(self):
		"""docstring for plot"""
        	plt.xlabel(self.xname)
        	plt.ylabel(self.yname)
		plt.grid(True)
        	plt.show()
