from PyQt4.uic import loadUiType

from lib.grid import *

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import(
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

Ui_MainWindow, QMainWindow = loadUiType('window.ui')

class Main(QMainWindow, Ui_MainWindow):
    def __init__(self, ):
        super(Main, self).__init__()
        self.setupUi(self)

    def addmpl(self, fig):
        self.canvas = FigureCanvas(fig)
        self.mplvl.addWidget(self.canvas)
        self.canvas.draw()
        self.toolbar = NavigationToolbar(self.canvas, 
            self.mplwindow, coordinates=True)
        self.mplvl.addWidget(self.toolbar)


    def rmmpl(self,):
        self.mplvl.removeWidget(self.canvas)
        self.canvas.close()
        self.mplvl.removeWidget(self.toolbar)
        self.toolbar.close()
    

if __name__ == '__main__':
    import sys
    from PyQt4 import QtGui
    import numpy as np

         
    fig1 = Figure()
    ax1f1 = fig1.add_subplot(111)
    l0=line(10,corner(0.0,0.0),corner(1.1,1.1),0)
    ax1f1.add_line(l0)
    ax1f1.figure.canvas.draw()

    fig2 = Figure()
    ax1f2 = fig2.add_subplot(121)
    ax1f2.plot(np.random.rand(5))
    ax2f2 = fig2.add_subplot(122)
    ax2f2.plot(np.random.rand(10))

    app = QtGui.QApplication(sys.argv)
    main = Main()
    main.addmpl(fig1)
    main.show()
    input()
    main.rmmpl()
    main.addmpl(fig2)
    main.show()
    sys.exit(app.exec_())
