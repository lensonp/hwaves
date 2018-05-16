from functools import partial

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import tkinter
from tkinter import Tk, Frame, Canvas, Button, Label, Entry, \
OptionMenu, Scrollbar, Checkbutton, \
StringVar, DoubleVar, BooleanVar, IntVar
import numpy as np

if __name__=='__main__':
    run_hwaves_gui()

def run_hwaves_gui():
    self.gui = Tk()
    self.gui.title('hydrogenic wavefunction plotter')
    
    self.plot_canvas = Canvas(self.gui, width=800, height=1600)
    self.plot_canvas.pack(side=tkinter.LEFT, fill=tkinter.BOTH, expand=tkinter.YES)

    self.control_frame = Frame(self.gui,relief=tkinter.RAISED)

    self.draw_plots()
    # start the tk loop
    self.gui.mainloop()

def draw_plots(self):
    self.ax_rwf.clear()
    self.ax_surf.clear()
    #self.ax_plot.loglog(self.q_I[:,0],self.q_I[:,1],lw=2,color='black')
    I_est = xrsdkit.scattering.compute_intensity(self.q_I[:,0],self.populations,self.src_wl)
    self.ax_plot.semilogy(self.q_I[:,0],I_est,lw=2,color='red')
    #self.ax_plot.loglog(self.q_I[:,0],I_est,lw=2,color='red')
    self.ax_plot.set_xlabel('q (1/Angstrom)')
    self.ax_plot.set_ylabel('Intensity (counts)')
    self.ax_plot.legend(['measured','computed'])
    self.plot_canvas.draw()



