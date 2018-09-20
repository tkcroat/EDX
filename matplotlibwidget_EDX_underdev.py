# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 10:35:25 2014

@author: eegroopm
"""
import numpy as np

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle, Circle, Arrow

from PyQt4 import QtGui
from PyQt4.QtCore import pyqtSlot, pyqtSignal

class MplCanvas(FigureCanvas):

    def __init__(self):
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)
        super(MplCanvas, self).__init__(self.fig)
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)        

class matplotlibWidget(QtGui.QWidget):
    distances = pyqtSignal(str,str,str,str)
    def __init__(self, common, Diffraction, parent = None):
        QtGui.QWidget.__init__(self, parent)
        self.update(EDXdataset) # necessary?
        
        self.canvas = MplCanvas()
        self.vbl = QtGui.QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.setLayout(self.vbl)
        self.Plot_initialize()
        self.canvas.mpl_connect('pick_event', self.on_pick)
        self.x1 = None; self.y1 = None
        
    def update(self, Diffraction):
        self.DSpaces = self.common.DSpaces
        self.Forbidden = self.common.Forbidden
        self.u = self.common.u
        self.v = self.common.v
        self.w = self.common.w
        self.E = self.common.beamenergy
        self.L = self.common.camlength
        self.const = self.common.camconst
        self.lam = self.common.wavelength
        
        self.ZoneAxis = self.common.ZoneAxis
        self.Diffraction = Diffraction
        #self.canvas.mpl_connect('button_press_event', self.on_pick)
#        
#    def setupToolbar(self,canvas,frame):
#        """Setup a custom toolbar"""
#        # Create the navigation toolbar, tied to the canvas
#        self.mpl_toolbar = NavigationToolbar(canvas, frame)
#        #add widgets to toolbar
#        self.comboBox_rotate = QtGui.QComboBox()
#        self.checkBox_labels = QtGui.QCheckBox()
#        self.mpl_toolbar.addWidget(self.comboBox_rotate)
#        self.mpl_toolbar.addWidget(self.checkBox_labels)
#        #add toolbar to tabs
#        self.verticalLayout.addWidget(self.mpl_toolbar)
        
        
    def Plot_initialize(self):
        """Initialize parameters of Matplotlib widget such as axes labels"""
        #label = u'Distance (\u212B\u207B\u00B9)'
        label = r'Distance ($\AA^{-1}$)' #use matplotlib's mathtex rendering: Å⁻¹
        self.canvas.ax.set_xlabel(label,fontsize=14)
        self.canvas.ax.set_ylabel(label,fontsize=14)
        self.canvas.ax.tick_params(axis='both', which='major', labelsize=14, length=6)
        #self.Plot.xaxis.set_units(u'Å⁻¹')
        #self.Plot.yaxis.set_units(u'Å⁻¹')
        self.canvas.fig.tight_layout()
        
    def calc(self,ind1,ind2):
        """Calculates angles from picks"""
        
        p1 = list(self.DSpaces.loc[ind1,['x','y']])
        p2 = list(self.DSpaces.loc[ind2,['x','y']])
        
        recip_d = round(np.sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2),3) #calc distance
        real_d = 1.0/recip_d
        film_d = self.lam*self.L/real_d*self.const
        
        #angle = round(np.degrees(self.Diffraction.AngleAmbiguity(p2[0]-p1[0],p2[1]-p1[1])),1)
        angle = round(np.degrees(np.arctan2((p2[1]-p1[1]),(p2[0]-p1[0]))),2)
        
        return recip_d, real_d,film_d, angle, p1, p2
        
    @pyqtSlot()
    def on_done_pick(self,recip_d, real_d,film_d, angle):
        self.distances.emit(recip_d, real_d,film_d, angle)
        
    
    def on_pick(self, event):
        # The event received here is of the type
        # matplotlib.backend_bases.PickEvent
        #self.DSpaces = self.common.DSpaces

        if isinstance(event.artist, Line2D):
            thisline = event.artist
            
            if self.x1 == None:
                if self.common._x2:
                    #Remove recently done circles
                    self.arr.remove()
                    del self.arr
                    l = self.canvas.ax.lines.pop(-1)
                    del l
                    l = self.canvas.ax.lines.pop(-1)
                    del l
                    self.canvas.draw()
                    self.common._x2=False
                self.x1 = thisline.get_xdata()
                self.y1 = thisline.get_ydata()
                self.ind1 = event.ind[0]
                self.canvas.ax.plot(self.x1[self.ind1],self.y1[self.ind1], linestyle = '', marker='o', markersize = 10,color='r')
                self.canvas.draw()
                
                
            elif self.x1 != None:
                self.update(self.common,self.Diffraction)
                self.common._x2 = True
                self.ind2 = event.ind[0]
                #make names shorter
                #x1 = self.x1[self.ind1]; x2 = self._x2[self.ind2]; y1 = self.y1[self.ind1]; y2 = self.y2[self.ind2]
                recip_d, real_d,film_d, angle, p1, p2 = self.calc(self.ind1,self.ind2)
                
                #reset x1 and y1
                self.x1 = None; self.y1 = None
                #plot colored circle
                self.canvas.ax.plot(p2[0],p2[1], linestyle = '', marker='o', markersize = 10,color='r')
                #plot arrow between selected points
                self.arr = Arrow(p1[0],p1[1],p2[0]-p1[0],p2[1]-p1[1],facecolor='r',width = 1/(5*self.common.a)) #maybe include factor of miller indices. higher miller = larger x,yrange
                self.canvas.ax.add_patch(self.arr)
                self.canvas.draw()
                
                self.on_done_pick(str(recip_d), str(round(real_d,2)),str(round(film_d,2)), str(angle))


import sys
from PyQt4 import QtGui, QtCore

# simple test of spinbox and connections
class spindemo(QtGui.QWidget):
   def __init__(self, parent = None):
      super(spindemo, self).__init__(parent)
      
      layout = QtGui.QVBoxLayout()
      self.l1 = QtGui.QLabel("current value:")
      self.l1.setAlignment(QtCore.Qt.AlignCenter)
      layout.addWidget(self.l1)
      self.sp = QtGui.QSpinBox()
		
      layout.addWidget(self.sp)
      self.sp.valueChanged.connect(self.valuechange)
      self.setLayout(layout)
      self.setWindowTitle("SpinBox demo")
		
   def valuechange(self):
      self.l1.setText("current value:"+str(self.sp.value()))
      print('New value is', str(self.sp.value()))

def main():
   app = QtGui.QApplication(sys.argv)
   ex = spindemo()
   ex.show()
   sys.exit(app.exec_())
	
if __name__ == '__main__':
   main()
   

if QtCore.QCoreApplication.instance() != None:
    app = QtCore.QCoreApplication.instance()
else:
    app = QtGui.QApplication(sys.argv)
mygui = spindemo()
mygui.show()
app.exec_()