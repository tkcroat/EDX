# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 15:34:59 2017

@author: tkc
"""

# from PyQt5.uic import loadUiType
import numpy as np 
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import ( FigureCanvasQTAgg as FigureCanvas, 
                                                NavigationToolbar2QT as NavigationToolbar)
from PyQt5 import QtCore, QtWidgets
from PyQt5 import QtGui # not clear if this will be used 
from PyQt5.uic import loadUI # when loading qtdesigner file
# pyqt5 splits pyqt4 QtGui into QtGui and QtWidgets
# Ui_MainWindow, QMainWindow = loadUiType('mplwindow.ui') # returns GUI application class & custom base class

# below is a rather confusing way to approach embedded mpl since details are hidden in the ui 
# file from qt designer 

class EDXmain(QtWidgets.QWidget):
    ''' Interface to get a list of points from an existing plot '''
    # http://blog.rcnelson.com/building-a-matplotlib-gui-with-qt-designer-part-1/
    def __init__(self, ):
        super(EDXmain, self).__init__() # gets inherited methods of this class
        self.initUI()
        self.fig_dict={} # list of active figures
        # click on item in Qlistwidget triggers change of figures
        self.mplfigs.itemClicked.connect(self.changefig)
        # add select data subset button here 
    
    def initUI(self):
        self.setGeometry(600, 300, 400, 200)
        self.setWindowTitle('EDX main window')     
        self.show()
        
    def addmpl(self, fig):
        ''' adds matplotlib plot to qt container'''
        self.canvas=FigureCanvas(fig) # fig canvas widget
        # mplvl is a QVBoxLayout instance with addWidget method
        self.mplvl.addWidget(self.canvas) # mplvl is layout name within container in .ui window
        # could also be directly added to mplwindow container but then scales differently
        self.canvas.draw()
        self.toolbar = NavigationToolbar(self.canvas, self.mplwindow, coordinates=True)
        self.mplvl.addWidget(self.toolbar) # toolbar added below passed plot
        self.canvas.mpl_connect('key_press_event', self.on_key_press)
       
    def rmmpl(self,):
        ''' remove existing plots  '''
        self.mplvl.removeWidget(self.canvas)
        self.canvas.close()
        self.mplvl.removeWidget(self.toolbar)
        self.toolbar.close()

    def addfig(self, name, fig):
        ''' Adds name, underlying fig to dictionary and displays in window '''
        self.fig_dict[name]=fig
        self.mplfigs.addItem(name) # adds name to QlistWidget
    
    def changefig(self, item):
        text = item.text()
        self.rmmpl() # removes existing
        self.addmpl(self.fig_dict[text])
        
    def on_key_press(self, event):
        ''' get coords on plot '''
        ix, iy = event.xdata, event.ydata
        print ('Coords are: ', str(ix), str(iy))

        
if __name__ == '__main__': # direct run option(rarely used)
    import sys
    from PyQt4 import QtGui
    app = QtGui.QApplication(sys.argv) # Qt GUI event loop
    sys.exit(app.exec_())
    main = Main()
    main.show()
    sys.exit(app.exec_)
    
main = EDXmain() # creates qt4 GUI 
main.show()

fig1=Figure() # matplotlib figure
ax1f1=fig1.add_subplot(111) 
ax1f1.plot(np.random.rand(5))

fig2=Figure() # matplotlib figure
ax1f2=fig2.add_subplot(111) 
ax1f2.plot(np.random.rand(5))

main.addmpl(fig1)
main.addfig('Figure1',fig1)
main.addmpl(fig2)
main.addfig('Figure2',fig2)
main.rmmpl()
main.close() # close method inherited by all QT windows 
