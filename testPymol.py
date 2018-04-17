import pymol
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends import backend_tkagg

def _new_figure_manager(num, *args, **kwargs):
    # import pymol
    #if pymol._ext_gui is None:
    #    return new_figure_manager(num, *args, **kwargs)
    #backend_tkagg.show._needmain = False
    import Tkinter as Tk
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, FigureManagerTkAgg
    FigureClass = kwargs.pop('FigureClass', Figure)
    print kwargs
    figure = FigureClass(*args, **kwargs)
    window = Tk.Toplevel(master=pymol._ext_gui.root)
    canvas = FigureCanvasTkAgg(figure, master=window)
    figManager = FigureManagerTkAgg(canvas, num, window)
    if matplotlib.is_interactive():
        figManager.show()
    return figManager

new_figure_manager = backend_tkagg.new_figure_manager
backend_tkagg.new_figure_manager = _new_figure_manager

matplotlib.interactive(True)
from matplotlib.figure import Figure
from matplotlib.pyplot import *

import numpy
x = numpy.linspace(1, 300, 300)
y = numpy.linspace(1, 300, 300)

fig1 = figure()
plot = fig1.add_subplot(111)
plot.scatter(x, y)

fig2 = figure()
plot = fig2.add_subplot(111)
plot.scatter(x, 2*y)
