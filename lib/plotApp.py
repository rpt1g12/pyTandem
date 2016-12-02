import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from plots2d import plot2d
import Tkinter as tk
from Tkinter import *
import matplotlib.pyplot as plt
import numpy as np
import indx

LARGE_FONT=('Verdana',12)
NORM_FONT=('Verdana',12)
SMALL_FONT=('Verdana',12)

def popupmsg(msg):
	"""This creates a popup message window"""
	popup=tk.Tk()

	popup.wm_title('Message')
	popup.geometry('200x100')
	label=tk.Label(popup, text=msg, font=NORM_FONT)
	label.pack(pady=10)
	b1=tk.Button(popup, text='OK', command=popup.destroy)
	b1.pack()
	popup.mainloop()

class plotApp(tk.Tk):
	"""docstring for plotApp"""
	def __init__(self, *args, **kwargs):
		tk.Tk.__init__(self,*args,**kwargs)

		tk.Tk.wm_title(self,'Plotting Application')

		container=tk.Frame(self)
		container.pack(side="top",fill="both",expand=True)
		container.grid_rowconfigure(0, weight=1)
		container.grid_columnconfigure(0, weight=1)

		menubar=tk.Menu(container)
		filemenu=tk.Menu(menubar, tearoff=0)
		filemenu.add_command(label='test',
		command=lambda:popupmsg('Sorry,\nNot Implemented Yet'))
		filemenu.add_command(label='Quit', command=self.quit)
		menubar.add_cascade(label='File',menu=filemenu)
		tk.Tk.config(self,menu=menubar)

		self.frames={}

		for F in (startPage,plotPage):
			frame = F(container,self)
			self.frames[F]= frame
			frame.grid(row=0, column=0, sticky='nsew')

		self.show_frame(startPage)

	def show_frame(self, cont):
		"""docstring for show_frame"""
		frame = self.frames[cont]
		frame.tkraise()


class startPage(tk.Frame):
	"""docstring for startPage"""
	def __init__(self, parent, controller):
		tk.Frame.__init__(self,parent)
		lable = tk.Label(self,text='Start Page', font=LARGE_FONT)
		lable.pack(pady=10, padx=10)

		b1 =tk.Button(self,text='Graph Page',
		              command=lambda:controller.show_frame(plotPage))
	 	b1.pack()

		b2 =tk.Button(self,text='Exit',
		              command=parent.quit)
	 	b2.pack()


class plotPage(tk.Frame):
	"""docstring for page1"""
	def __init__(self, parent, controller):
		tk.Frame.__init__(self,parent,)
		lable1 = tk.Label(self,text='Graph Page', font=LARGE_FONT)
		lable1.pack()


		def showPlot(event):
			"""docstring for showPlot"""
		   	self.filename=e1.get()
			self.k=int(e2.get())
			p1=plot2d(self.filename)
			grid=bool(self.pgrid.get())
			legend=bool(self.legend.get())
			p1.plot(self.k,grid=grid,legend=legend)
			self.fig=p1.f
			self.plot=p1.a
			if self.invy.get():
				self.fig.gca().invert_yaxis()

			#Draw plot
			self.canvas.get_tk_widget().pack_forget()
			self.canvas=FigureCanvasTkAgg(self.fig,self.plotFrame)
			self.canvas.show()
			self.canvas.get_tk_widget().pack(side=tk.TOP,
			padx=10,fill=tk.BOTH, expand=True)

			#Add plot toolbar
			self.toolbar.pack_forget()
			self.toolbar=NavigationToolbar2TkAgg(self.canvas,self.plotFrame)
			self.toolbar.update()
			self.canvas._tkcanvas.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)



		#Create Frame for settings
		self.setFrame=tk.Frame(self)
		self.setFrame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True,padx=2)
		#Select file for data extracting
		lable2 = tk.Label(self.setFrame,text='Data File:', font=LARGE_FONT)
		lable2.grid(row=0,column=0,sticky='e')
		self.filename='Cf.dat'
		e1=tk.Entry(self.setFrame,width=12)
		e1.insert(0,self.filename)
		e1.bind('<Return>',showPlot)
		e1.grid(row=0,column=1,sticky='w')
		#Selent the span location to plot
		lable3 = tk.Label(self.setFrame,text='Span Location:', font=LARGE_FONT)
		lable3.grid(row=1,column=0,sticky='e')
		self.k=0
		e2=tk.Entry(self.setFrame,width=5)
		e2.insert(0,'0')
		e2.bind('<Return>',showPlot)
		e2.grid(row=1,column=1,sticky='w')
		#Use grid
		lable4 = tk.Label(self.setFrame,text='Grid:', font=LARGE_FONT)
		lable4.grid(row=2,column=0,sticky='e')
		self.pgrid=BooleanVar()
		c1=tk.Checkbutton(self.setFrame,variable=self.pgrid,
		onvalue=True,offvalue=False)
		c1.grid(row=2,column=1,sticky='w')
		#Use legend
		lable5 = tk.Label(self.setFrame,text='Legend:', font=LARGE_FONT)
		lable5.grid(row=3,column=0,sticky='e')
		self.legend=True
		self.legend=BooleanVar()
		c2=tk.Checkbutton(self.setFrame,variable=self.legend,
		onvalue=True,offvalue=False)
		c2.grid(row=3,column=1,sticky='w')
		#Invert Y-Axis
		lable6=tk.Label(self.setFrame,text='Invert Y-Axis:',font=LARGE_FONT)
		lable6.grid(row=4,column=0,sticky='e')
		self.invy=BooleanVar()
		c3=tk.Checkbutton(self.setFrame,variable=self.invy,onvalue=True,
		offvalue=False)
		c3.grid(row=4,column=1,sticky='w')
		#Create draw button
		b3=tk.Button(self.setFrame,text='Draw!')
		b3.bind('<Button-1>',showPlot)
		b3.grid(row=5,column=0,columnspan=2,sticky='we')

		#Create frame for plotting
		self.plotFrame=tk.Frame(self)
		self.plotFrame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
		self.fig=Figure()
		self.plot=self.fig.add_subplot()
		self.canvas=FigureCanvasTkAgg(self.fig,self.plotFrame)
		self.canvas.show()
		self.canvas.get_tk_widget().pack(side=tk.TOP,
		padx=10,fill=tk.BOTH, expand=True)
		self.toolbar=NavigationToolbar2TkAgg(self.canvas,self.plotFrame)
		self.canvas._tkcanvas.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

		#Create frame for bottom buttons
		self.buttonFrame=tk.Frame(self)
		self.buttonFrame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
		#Go Home button
		b1 =tk.Button(self.buttonFrame,text='Home',
		              command=lambda:controller.show_frame(startPage))
		b1.pack(side=tk.LEFT,pady=10)
		#Exit button
		b2 =tk.Button(self.buttonFrame,text='Exit',
		              command=parent.quit)
		b2.pack(side=tk.LEFT,pady=10)


		



		
