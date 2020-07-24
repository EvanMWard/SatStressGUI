#!/usr/bin/python
# -*- coding: utf-8 -*-

# for GUI
#    wx = python extention module that acts as python lang wrapper for wxWidgets
#         (cross platform GUI API written in C++)t

# for manipulating tabluar data
#     wx.grid = allows displaying, editing, customizing tabular data'''
import wx.grid

#     copy/Users/andreismailyan = allows shallow & deep copying operations
#     sys = access syst specific parameters & fcts, variables held by interpreter
#     os = allows more direct interaction with OS
import sys, os

# for trying to print stack traces under program control (i.e. wrapper), probs used for wx?
#     traceback = allows printing/extracting/formatting Python program stack traces
import traceback

# for plotting capabilities
#    matplotlib = 2D
#    scipy.ndimage = multi-D
import matplotlib, scipy.ndimage

# Allows the program to find the time and date.
# Used to save images to a unique file location. - PS 2016
import time

# fig state updated every plot command, but only redraw on explicit calls to draw()
matplotlib.interactive(False)
# sets matplotlib backend to 'WXAgg'
matplotlib.use('WXAgg')

# for plotting
#    basemap allows cartographic
from mpl_toolkits import basemap
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas, \
                                              NavigationToolbar2Wx as NavigationToolbar
from matplotlib.figure import Figure

# for manipulating netCDF files
#import netCDF3
#Does not run in Windows, so we've commented it out here to make for easy copying. -PS 2016

# for math functions
import numpy

# satstress library
from satstress.satstress import *
from satstress.gridcalc import *
from satstress.lineament import plotlinmap, Lineament, lingen_nsr, shp2lins, lins2shp
from satstress.cycloid import Cycloid, SaveCycloidAsShape
from satstress.stressplot import scalar_grid
import satstress.physcon

# Used for mapping and saving to a shapefile.

from SatelliteCalculation import SatelliteCalculation

from CycloidTab import *
from GridTab import *
from PointTab import *
from SatelliteTab import *
from StressesTab import *
# constants set as global variables
seconds_in_year = 31556926.0  # 365.24 days
vector_mult = 4000
display_dpi = 72
scale_left = 0.25
scale_bar_length = 0.38

# ===============================================================================
# Class containing overhead functions for configuring, used in PlotPanel
# ===============================================================================
class Config:
    """
    Class the holds application settings --> specifically?
    """
    #default step for plots <-- what units?
    default_step = 30

    def __init__(self, configfile='config'):
        self.configfile = configfile
        self.conf = {}

    # a is optional arg
    def load(self, *a):
        try:
            c = open(self.configfile)
            self.conf = nvf2dict(c)
            c.close()
            ret = filter(lambda x: x, map(self.conf.get, a))
            if len(a) == 1 and len(ret) == 1:
                return ret[0]
            else:
                return ret
        except:
            self.conf = {}
            return []

    # **kw unpacks the extra dictionary args
    def save(self, **kw):
        for k, v in kw.items():
            self.conf[k] = v   # conf is dictionary
        try:
            c = open(self.configfile, 'w')
            c.writelines([ "%s = %s\n" % (k,v) for k,v in self.conf.items() ])
            c.close()
        except:
            pass

    def load_step(self, step_field='STEP'):
        self.load()
        try:
            return float(self.conf.get(step_field, self.default_step))
        except:
            return self.default_step

    def save_step(self, step, step_field='STEP'):
        self.conf[step_field] = step
        self.save()

#creates a global instance of config
config = Config()


# ===============================================================================
# Global Function for tensile and compressive components, used in ScalarPlotPanel
# ===============================================================================
def vector_points1(stresscalc=None, lons=None, lats=None, time_t=0.0,\
    plot_tens=True, plot_comp=True, plot_greater=True, plot_lesser=True,\
    basemap_ax=None, lonshift=0, w_stress=False,\
    scale=1e8, scale_arr=None, arrow_width=0.008):
    """
    This function and vector_points2 were developed from stressplot.vector_points.

    Display the principal components of the tidal stresses defined by the input
    stresscalc object at the points defined by lons and lats, which are one
    dimensional arrays of equal length, and a time defined by time_t, in
    seconds after periapse.

    The stress vectors are plotted on the map axes defined by basemap_ax.

    By default all the principal components are plotted, but if you wish to see
    only the more or less tensile (less or more compressive) or only those
    principal components which are absolutely compressive or tensile, you may
    exclude some subset of the vectors using the following flags:

       plot_tens  --  if True, plot all tensile stresses.
       plot_comp  --  if True, plot all compressive stresses.
    plot_greater  --  if True, plot the greater (more tensile) principal component
     plot_lesser  --  if True, plot the lesser (more compressive) principal component

    lonshift is a longitudinal displacement added to lons when the stresses are
    calculated, useful in creating plots of lineaments at their current
    location, compared to stresses that they would have experienced at their
    apparent location of formation (i.e. those stresses which best match the
    feature).  For instance, if you wished to show only those stresses which are 
    the more tensile, and which are actually tensile, you would need to set
    the flags: plot_comp=False, plot_lesser=False.

    If w_stress is true, the lengths of the arrows which are used to represent
    the stresses are scaled according to how significant their location within
    the stress field is, i.e. large stresses and anisotropic stresses will be
    more prominent than small stresses and isotropic stresses.

    scale determines the overall size of the arrows representing the stresses.
    A smaller scale means bigger arrows.

    scale_arr is an array of the same length as lons and lats, which is used to
    scale the lengths of the vectors.  Useful in showing the relative
    importance of different segments of a feature having non-uniform lengths.

    arrow_width is passed in to numpy.quiver(), and is the width of the arrow
    shaft, as a proportion of the width of the overall plot.
    """

    # There must be something wrong in this function or in the Polar Wander calculations which makes the vectors point in the wrong direction.
    # There is nothing in this function which makes it dependant on the stress, and the vectors for the other stresses appear correct.
    # We don't know what is causing only the Polar Wander stresses to generate incorrectly.  -PS 2016

    calc_phis   = lons
    calc_thetas = (numpy.pi / 2.0) - lats

    Ttt, Tpt, Tpp = stresscalc.tensor(calc_thetas, calc_phis + lonshift, time_t)

    Tau = numpy.array([[Ttt, Tpt], [Tpt, Tpp]])
    eigensystems = [ numpy.linalg.eig(Tau[:,:,N]) for N in range(len(Tau[0,0,:])) ]
    evals = numpy.array([ e[0] for e in eigensystems ])
    evecs = numpy.array([ e[1] for e in eigensystems ])

    eigval_A = evals[:,0]
    ex_A     = evecs[:,0,1]
    ey_A     = evecs[:,0,0]

    eigval_B = evals[:,1]
    ex_B     = evecs[:,1,1]
    ey_B     = evecs[:,1,0]

    mag1 = numpy.where(eigval_A >  eigval_B, eigval_A, eigval_B)
    ex1  = numpy.where(eigval_A >  eigval_B, ex_A, ex_B)
    ey1  = numpy.where(eigval_A >  eigval_B, ey_A, ey_B)

    mag2 = numpy.where(eigval_A <= eigval_B, eigval_A, eigval_B)
    ex2  = numpy.where(eigval_A <= eigval_B, ex_A, ex_B)
    ey2  = numpy.where(eigval_A <= eigval_B, ey_A, ey_B)

    if numpy.shape(scale_arr) != numpy.shape(mag1):
        scale_arr = numpy.ones(numpy.shape(mag1))
    if numpy.shape(w_stress) == numpy.shape(mag1):
        scale_arr = scale_arr*(mag1 - mag2)/stresscalc.mean_global_stressdiff()

    mag1_comp = numpy.ma.masked_where(mag1 > 0, mag1)
    mag1_tens = numpy.ma.masked_where(mag1 < 0, mag1)

    mag2_comp = numpy.ma.masked_where(mag2 > 0, mag2)
    mag2_tens = numpy.ma.masked_where(mag2 < 0, mag2)

    scaled = {}
    scaled[1, 'comp'] = mag1_comp*scale_arr
    scaled[2, 'comp'] = mag2_comp*scale_arr
    scaled[1, 'tens'] = mag1_tens*scale_arr
    scaled[2, 'tens'] = mag2_tens*scale_arr

    ex = { 1: ex1, 2: ex2 }
    ey = { 1: ey1, 2: ey2 }
    # map coordinates
    dlons, dlats = numpy.degrees(lons), numpy.degrees(lats)
    x,y = basemap_ax(dlons, dlats)
    # new basis
    exx,exy = basemap_ax.rotate_vector(numpy.ones(numpy.shape(lons)), numpy.zeros(numpy.shape(lons)), dlons, dlats)
    eyx,eyy = basemap_ax.rotate_vector(numpy.zeros(numpy.shape(lats)), numpy.ones(numpy.shape(lats)), dlons, dlats)

    rotated = {}
    for i in range(1,3):
        for s in ['comp', 'tens']:
            x1 = scaled[i,s] * ex[i]
            y1 = scaled[i,s] * ey[i]
            rotated[i,s,'x'], rotated[i,s,'y'] = x1*exx + y1*eyx, x1*exy + y1*eyy

    # where is the exclusion done?
    for i in range(1,3):
        for k in range(2):
            basemap_ax.quiver(x, y, (-1)**k*rotated[i,'comp','x'], (-1)**k*rotated[i,'comp','y'],
                    lw=0., width=arrow_width, scale=scale, color='blue', pivot='tip', units='inches')
            basemap_ax.quiver(x, y, (-1)**k*rotated[i,'tens','x'], (-1)**k*rotated[i,'tens','y'],
                    lw=0., width=arrow_width, scale=scale, color='red', pivot='tail', units='inches')

# ===============================================================================
# Global Function for shear and normal components
# ===============================================================================
def vector_points2(stresscalc=None, lons=None, lats=None, time_t=0.0,\
    plot_norm_lon=True, plot_norm_lat=True, plot_shear=True, \
    basemap_ax=None, lonshift=0, \
    scale=1e8, arrow_width=0.008):
    """
    This function and vector_points2 were developed from stressplot.vector_points.

    Display the normal and shear components of the tidal stresses defined by the input
    stresscalc object at the points defined by lons and lats, which are one
    dimensional arrays of equal length, and a time defined by time_t, in
    seconds after periapse.

    The stress vectors are plotted on the map axes defined by basemap_ax.

    lonshift is a longitudinal displacement added to lons when the stresses are
    calculated, useful in creating plots of lineaments at their current
    location, compared to stresses that they would have experienced at their
    apparent location of formation (i.e. those stresses which best match the
    feature).  For instance, if you wished to show only those stresses which are 
    the more tensile, and which are actually tensile, you would need to set
    the flags: plot_comp=False, plot_lesser=False.

    scale determines the overall size of the arrows representing the stresses.
    A smaller scale means bigger arrows.

    arrow_width is passed in to numpy.quiver(), and is the width of the arrow
    shaft, as a proportion of the width of the overall plot.

    """

    calc_phis   = lons
    # because should be co-latitudal (south-positive, [0, pi])
    calc_thetas = (numpy.pi/2.0)-lats
    # plot coordinates
    dlons, dlats = numpy.degrees(lons), numpy.degrees(lats)
    px,py = basemap_ax(dlons, dlats)
    # new basis
    exx,exy = basemap_ax.rotate_vector(numpy.ones(numpy.shape(lons)), numpy.zeros(numpy.shape(lons)), dlons, dlats)
    eyx,eyy = basemap_ax.rotate_vector(numpy.zeros(numpy.shape(lats)), numpy.ones(numpy.shape(lats)), dlons, dlats)

    Ttt, Tpt, Tpp = stresscalc.tensor(calc_thetas, calc_phis+lonshift, time_t)

    def plot_vector(v, x, y, scale, color, pivot):
        vx, vy = v*x*exx + v*y*eyx, v*x*exy + v*y*eyy
        basemap_ax.quiver(px, py,  vx,  vy,\
            lw=0., width=arrow_width, scale=scale, color=color, pivot=pivot, units='inches')
        basemap_ax.quiver(px, py, -vx, -vy,\
            lw=0., width=arrow_width, scale=scale, color=color, pivot=pivot, units='inches')

    def plot_vectors(vs, x, y):
        plot_vector(numpy.ma.masked_where(vs > 0, vs), x, y, scale, 'blue', 'tip')
        plot_vector(numpy.ma.masked_where(vs < 0, vs), x, y, scale, 'red', 'tail')

    if plot_norm_lat:
        plot_vectors(Ttt, 0, 1)
    if plot_norm_lon:
        plot_vectors(Tpp, 1, 0)
    if plot_shear:
        for diag_angle in [numpy.pi/4, 3*numpy.pi/4]:
            plot_vectors(-Tpt, numpy.cos(diag_angle), 1 - numpy.sin(diag_angle))





# ==================================================================================
# MAIN GUI INTERFACE: TABS, PANELS, PLOTS
# ==================================================================================






# ===============================================================================
# PLOT TAB: HELPER CLASSES AND FUNCTIONS
# ===============================================================================
"""
Polar Wander slider is currently disabled.  All of the code has been left in, but it is commented out.
If the calculations for Polar Wander are improved to be viscoelastic, they can be re-implemented.
-Peter Sinclair, 2016
"""
class StepSlider(matplotlib.widgets.Slider):
    """
    Custom designed class for discrete slider control at bottom of plot panel to control 
    satellite's orbital position.

    Used in add_orbit_controls and add_nsr_controls
    """
    def __init__(self, ax, label, valmin, valmax, numsteps, *args, **kw):
        self.steps_n = numsteps
        self.updating = False
        self.prev_val = kw.get('valinit', 0)
        matplotlib.widgets.Slider.__init__(self, ax, label, valmin, valmax, *args, **kw)
        ax.lines.remove(self.vline)

    def on_changed_f(self, val):
        pass

    def on_changed(self, f):
        def f2(val):
            if self.updating:
                return
            self.eventson = False
            self.updating = True
            val += self.valmin
            self.set_stepval(val)
            f(self.val)
            self.updating = False
            self.eventson = True
        self.on_changed_f = f
        matplotlib.widgets.Slider.on_changed(self, f2)

    def set_stepval(self, val):
        if val < self.valmin:
            self.set_val(self.valmin)
        elif val > self.valmax:
            self.set_val(self.valmax)
        elif self.valmax - self.valmin > 0 and self.steps_n > 0:
            step = float(self.valmax - self.valmin)/self.steps_n
            n0 = int((val - self.valmin)/step)
            n1 = n0 + 1
            if abs(val - self.prev_val) > 0.7*step:
                self.prev_val = round((val - self.valmin)/step)*step + self.valmin
                self.set_val(self.prev_val)
            else:
                self.set_val(self.prev_val)

    def reset(self):
        self.updating = True
        matplotlib.widgets.Slider.reset(self)
        self.updating = False

    def first(self):
        self.set_stepval(self.valmin)
        # HERE initial_split()
        # ^ I think that these comments have to do with cycloid generation -PS 2016

    def last(self):
        self.set_stepval(self.valmax)
        # HERE EVERYTHING SHOULD BE GRAPHED

    def next(self):
        step = float(self.valmax - self.valmin)/self.steps_n
        n = int((self.val - self.valmin)/step) + 1
        self.set_stepval(n*step + self.valmin)
        # ONLY GRAPH UP TO THIS POINT

    def prev(self):
        step = float(self.valmax - self.valmin)/self.steps_n
        n = int((self.val - self.valmin)/step) - 1
        self.set_stepval(n*step + self.valmin)
        # ONLY GRAPH UP TO THIS POINT

class CustomPlotToolbar(NavigationToolbar):
    def __init__(self, plotCanvase):
        # create default toolbar
        NavigationToolbar.__init__(self, plotCanvase)

        # remove unwanted button
        # stress plot only exists in rectangular bounds
        # may need to add in later if want movable scale bar
        # or only pan when zoommed
        POSITION_OF_PANNING_BUTTON = 3

        # remove unnecessary button (no subplots)
        POSITION_OF_CONFIGURE_SUBPLOT_BUTTON = 6
        self.DeleteToolByPos(POSITION_OF_CONFIGURE_SUBPLOT_BUTTON)

class MatPlotPanel(wx.Panel):
    """
    GUI object that holds the plot area
    """
    def __init__(self, *args, **kw):
        super(MatPlotPanel, self).__init__(*args, **kw)
        self.figure = Figure(figsize=(6,5),dpi=display_dpi)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.ax = self.figure.add_subplot(111)

        toolbar = CustomPlotToolbar(self.canvas)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, flag=wx.EXPAND|wx.ALL)
        toolbar.Realize()
        sizer.Add(toolbar, flag=wx.EXPAND|wx.ALL)
        self.SetSizer(sizer)
        self.SetMinSize((625, 400))

    def get_axes(self):
        return self.ax

    def draw(self):
        self.canvas.draw()

    def colorbar(self, mappable, *a, **kw):
        #return self.figure.colorbar(mappable, ax=self.ax, *a, **kw)
        return self.figure.colorbar(mappable, *a, **kw)

class StressPlotPanel(MatPlotPanel):
    """
    Contains controls for going through the time frame dictated in "Grid" Tab.
    Specifically, the [< | < | > | >] controls
    """
    scale_y    = 0.15
    orbit_y    = 0.11
    polar_y = 0.01
    nsr_y    = 0.06
    button_l = 0.04
    bbutton_l= 0.12
    slider_h = 0.04
    slider_x = scale_left + scale_bar_length + button_l*2

    def __init__(self, *args, **kw):
        super(StressPlotPanel, self).__init__(*args, **kw)
        self.figure.subplots_adjust(bottom=0.25)
        # creates scale bar for the vectors (arrows) i.e. |-----| 91 kPa
        self.scale_ax = self.figure.add_axes([scale_left, self.scale_y, scale_bar_length, self.slider_h], frame_on=False)
        self.add_orbit()
        #self.add_polar()
        self.add_nsr()

    def get_ax_orbit(self):
        return self.figure.add_axes([scale_left, self.orbit_y, scale_bar_length, self.slider_h])

    def get_ax_nsr(self):
        return self.figure.add_axes([scale_left, self.nsr_y, scale_bar_length, self.slider_h])

    def get_ax_polar(self):
        return self.figure.add_axes([scale_left, self.polar_y, scale_bar_length, self.slider_h])

    def del_orbit(self):
        self.figure.delaxes(self.ax_orbit)
        self.del_orbit_controls()

    def del_nsr(self):
        self.figure.delaxes(self.ax_nsr)
        self.del_nsr_controls()

    def del_polar(self):
        self.figure.delaxes(self.ax_polar)
        self.del_polar_controls()

    def del_orbit_controls(self):
        for a in [self.ax_orbit_first, self.ax_orbit_prev, \
            self.ax_orbit_next, self.ax_orbit_last, self.ax_orbit_save]:
            self.figure.delaxes(a)

    def del_nsr_controls(self):
        for a in [self.ax_nsr_first, self.ax_nsr_prev, \
            self.ax_nsr_next, self.ax_nsr_last, self.ax_nsr_save]:
            self.figure.delaxes(a)

    def del_polar_controls(self):
        for a in [self.ax_polar_first, self.ax_polar_prev, \
            self.ax_polar_next, self.ax_polar_last, self.ax_polar_save]:
            self.figure.delaxes(a)

    def add_orbit(self):
        self.ax_orbit = self.get_ax_orbit()
        self.add_orbit_controls()

    def add_orbit_controls(self):
        x = self.slider_x
        self.ax_orbit_first = self.figure.add_axes([x, self.orbit_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_orbit_prev = self.figure.add_axes([x, self.orbit_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_orbit_next = self.figure.add_axes([x, self.orbit_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_orbit_last = self.figure.add_axes([x, self.orbit_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_orbit_save = self.figure.add_axes([x, self.orbit_y, self.bbutton_l, self.slider_h])

        # Note: StepSlider is custom designed class in/for gui
        self.orbit_slider = StepSlider(self.ax_orbit, 'Orbital position', 0, 1, 10, valinit=0, dragging=False)

        self.orbit_first_button = matplotlib.widgets.Button(self.ax_orbit_first, '[<')
        self.orbit_prev_button = matplotlib.widgets.Button(self.ax_orbit_prev, '<')
        self.orbit_next_button = matplotlib.widgets.Button(self.ax_orbit_next, '>')
        self.orbit_last_button = matplotlib.widgets.Button(self.ax_orbit_last, '>]')
        self.orbit_save_button = matplotlib.widgets.Button(self.ax_orbit_save, 'Save series')

        self.orbit_first_button.on_clicked(lambda e: self.orbit_slider.first())  # lambda functions
        self.orbit_prev_button.on_clicked(lambda e: self.orbit_slider.prev())    # ok, so are empirically necessary, but why?
        self.orbit_next_button.on_clicked(lambda e: self.orbit_slider.next())
        self.orbit_last_button.on_clicked(lambda e: self.orbit_slider.last())
        # hack
        self.orbit_save_button.on_clicked(lambda e: wx.CallLater(125, self.on_save_orbit_series, e))

    def add_polar(self):
        self.ax_polar = self.get_ax_polar()
        self.add_polar_controls()

    def add_polar_controls(self):
        x = self.slider_x
        self.ax_polar_first = self.figure.add_axes([x, self.polar_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_polar_prev = self.figure.add_axes([x, self.polar_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_polar_next = self.figure.add_axes([x, self.polar_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_polar_last = self.figure.add_axes([x, self.polar_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_polar_save = self.figure.add_axes([x, self.polar_y, self.bbutton_l, self.slider_h])

        # Note: StepSlider is custom designed class in/for gui
        self.polar_slider = StepSlider(self.ax_polar, 'Polar position', 0, 1, 10, valinit=0, dragging=False)

        self.polar_first_button = matplotlib.widgets.Button(self.ax_polar_first, '[<')
        self.polar_prev_button = matplotlib.widgets.Button(self.ax_polar_prev, '<')
        self.polar_next_button = matplotlib.widgets.Button(self.ax_polar_next, '>')
        self.polar_last_button = matplotlib.widgets.Button(self.ax_polar_last, '>]')
        self.polar_save_button = matplotlib.widgets.Button(self.ax_polar_save, 'Save series')

        self.polar_first_button.on_clicked(lambda e: self.polar_slider.first())
        self.polar_prev_button.on_clicked(lambda e: self.polar_slider.prev())
        self.polar_next_button.on_clicked(lambda e: self.polar_slider.next())
        self.polar_last_button.on_clicked(lambda e: self.polar_slider.last())
        self.polar_save_button.on_clicked(lambda e: wx.CallLater(125, self.on_save_polar_series, e))

    def add_nsr(self):
        self.ax_nsr = self.get_ax_nsr()
        self.add_nsr_controls()

    def add_nsr_controls(self):
        x = self.slider_x
        self.ax_nsr_first = self.figure.add_axes([x, self.nsr_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_nsr_prev = self.figure.add_axes([x, self.nsr_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_nsr_next = self.figure.add_axes([x, self.nsr_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_nsr_last = self.figure.add_axes([x, self.nsr_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_nsr_save = self.figure.add_axes([x, self.nsr_y, self.bbutton_l, self.slider_h])
        self.nsr_slider = StepSlider(self.ax_nsr, 'NSR position', 0, 1, 10, valinit=0, dragging=False, valfmt="%.1g")
        self.nsr_first_button = matplotlib.widgets.Button(self.ax_nsr_first, '[<')
        self.nsr_prev_button = matplotlib.widgets.Button(self.ax_nsr_prev, '<')
        self.nsr_next_button = matplotlib.widgets.Button(self.ax_nsr_next, '>')
        self.nsr_last_button = matplotlib.widgets.Button(self.ax_nsr_last, '>]')
        self.nsr_first_button.on_clicked(lambda e: self.nsr_slider.first())
        self.nsr_prev_button.on_clicked(lambda e: self.nsr_slider.prev())
        self.nsr_next_button.on_clicked(lambda e: self.nsr_slider.next())
        self.nsr_last_button.on_clicked(lambda e: self.nsr_slider.last())
        self.nsr_save_button = matplotlib.widgets.Button(self.ax_nsr_save, 'Save series')
        self.nsr_save_button.on_clicked(lambda e: wx.CallLater(125, self.on_save_nsr_series, e))

    def change_slider(self, ax, slider, label=None, valmin=None, valmax=None, numsteps=None, valinit=None, valfmt=None):
        if label is None:
            label = slider.label.get_text()
        if valmin is None:
            valmin = slider.valmin
        if valmax is None:
            valmax = slider.valmax
        if numsteps is None:
            numsteps = slider.numsteps
        if valinit is None:
            valinit = slider.valinit
        if valfmt is None:
            valfmt = slider.valfmt
        f = slider.on_changed_f
        slider = StepSlider(ax, label, valmin, valmax, numsteps, valinit=valinit, dragging=False, valfmt=valfmt)
        slider.on_changed(f)
        return slider

    def change_orbit_slider(self, valmin, valmax, numsteps, valinit=None):
        if valinit is None:
            valinit = valmin
        self.figure.delaxes(self.ax_orbit)
        self.ax_orbit = self.get_ax_orbit()
        self.orbit_slider = self.change_slider(
            self.ax_orbit, self.orbit_slider, valmin=valmin, valmax=valmax, numsteps=numsteps, valinit=valinit)

    def change_nsr_slider(self, valmin, valmax, numsteps, valinit=None):
        if valinit is None:
            valinit = valmin
        self.figure.delaxes(self.ax_nsr)
        self.ax_nsr = self.get_ax_nsr()
        self.nsr_slider = self.change_slider(
            self.ax_nsr, self.nsr_slider, valmin=valmin, valmax=valmax, numsteps=numsteps, valinit=valmin, valfmt="%.1g")

    def change_polar_slider(self, valmin, valmax, numsteps, valinit=None):
        if valinit is None:
            valinit = valmin
        self.figure.delaxes(self.ax_polar)
        self.ax_polar = self.get_ax_polar()
        self.polar_slider = self.change_slider(self.ax_polar, self.polar_slider, valmin=valmin, valmax=valmax, numsteps=numsteps, valinit=valmin, valfmt="%.1g")

    def plot_scale(self, scale, valfmt):
        self.scale_ax.clear()
        while self.scale_ax.texts:
            self.scale_ax.texts.pop()
        self.scale_ax.set_xticks([])
        self.scale_ax.set_yticks([])
        self.scale_ax.text(-0.02, 0.5, 'Scale', transform=self.scale_ax.transAxes, va='center', ha='right')
        self.scale_ax.text(0.23, 0.5, valfmt % scale, transform=self.scale_ax.transAxes, va='center', ha='left')
        self.scale_ax.plot([0.00, 0.20], [0.5, 0.5], linestyle='solid', marker='|', color='black', lw=1)
        self.scale_ax.set_xlim(0.0, 1.0)

    def on_save_orbit_series(self, evt):
        try:
            dir_dialog(None,
                message=u"Save calculation series on orbit period",
                style=wx.SAVE,
                action=self.save_orbit_series)
        except LocalError, e:
            error_dialog(self, str(e), e.title)

    def on_save_nsr_series(self, evt):
        try:
            dir_dialog(self,
                message=u"Save calculation series on nsr period",
                style=wx.SAVE,
                action=self.save_nsr_series)
        except LocalError, e:
            error_dialog(self, str(e), e.title)

    def on_save_polar_series(self, evt):
        try:
            dir_dialog(self,
                       message=u"Save calculation series on nsr period",
                       style=wx.SAVE,
                       action=self.save_polar_series)
        except LocalError, e:
            error_dialog(self, str(e), e.title)

class PlotPanel(SatPanel):
    """
    Manages drawing of plots
    """
    step_field = 'STEP'

    def __init__(self, *args, **kw):
        super(PlotPanel, self).__init__(*args, **kw)
        self.sc.parameters['projection'] = 'cyl'
        self.load_step()

    def add_stepspin(self, sz):
        sz.Add(wx.StaticText(self, label=u"Tick mark increment"), flag=wx.ALIGN_CENTER_VERTICAL)
        self.stepspin = wx.SpinCtrl(self, initial=int(self.step), min=0, max=180)
        sz.Add(self.stepspin, flag=wx.ALL|wx.EXPAND)
        self.stepspin.SetValue(self.step)
        self.stepspin.Bind(wx.EVT_SPINCTRL, self.adjust_step)

    def adjust_step(self, evt):
        self.adjust_coord_step(self.stepspin.GetValue())

    def load_step(self):
        self.step = config.load_step(self.step_field)
        return self.step
    def save_step(self):
        config.save_step(self.step, self.step_field)

    def plot(self):
        try:
            self.plot_no_draw()
            self.draw()
        except LocalError, e:
            if self.sc.satellite or self.sc.grid:
                error_dialog(self, str(e), e.title)
        except Exception, e:
            if self.sc.satellite or self.sc.grid:
                if not self.sc.get_stresses():
                    traceback.print_exc()
                    error_dialog(self, 'Stresses are not defined', 'Plot Error')
                else:
                    traceback.print_exc()
                    error_dialog(self, e.__class__.__name__ + ': ' + str(e), "Plot Error")

    def plot_no_draw(self):
        self.grid = self.sc.get_grid()
        self.calc = self.sc.get_calc()
        self.basemap_ax = self.get_basemap_ax()
        self.plot_grid_calc()
        self.draw_coords()
        if self.sc.parameters.get('Polar Wander', False) and self.sc.parameters['to_plot_pw_markers']:
            self.place_pw_coordinates()

    def place_pw_coordinates(self):
        #Places the user-given polar wander coordinates on the map.
        #Also places their antipodes (the anti-jove point and south pole). -PS 2016

        #Plot rotational poles if the coordinates of the initial and final pole differ
        if (self.sc.polarwander_coordinates['thetaRInitial'] != self.sc.polarwander_coordinates['thetaRFinal']
         or self.sc.polarwander_coordinates['phiRInitial'] != self.sc.polarwander_coordinates['phiRFinal']):

            self.basemap_ax.plot(self.sc.polarwander_coordinates['phiRInitial'],
                self.sc.polarwander_coordinates['thetaRInitial'],
                'wo', markersize=10)
            if self.sc.polarwander_coordinates['phiRInitial'] >=0:
                self.basemap_ax.plot(self.sc.polarwander_coordinates['phiRInitial'] - 180,
                    0 - self.sc.polarwander_coordinates['thetaRInitial'],
                    'wo', markersize=10)
            else:
                self.basemap_ax.plot(self.sc.polarwander_coordinates['phiRInitial'] + 180,
                    0- self.sc.polarwander_coordinates['thetaRInitial'],
                    'wo', markersize=10)

            self.basemap_ax.plot(self.sc.polarwander_coordinates['phiRFinal'],
                self.sc.polarwander_coordinates['thetaRFinal'],
                'ko', markersize=10)
            if self.sc.polarwander_coordinates['phiRFinal'] >=0:
                self.basemap_ax.plot(self.sc.polarwander_coordinates['phiRFinal'] - 180,
                    0 - self.sc.polarwander_coordinates['thetaRFinal'],
                    'ko', markersize=10)
            else:
                self.basemap_ax.plot(self.sc.polarwander_coordinates['phiRFinal'] + 180,
                    0- self.sc.polarwander_coordinates['thetaRFinal'],
                    'ko', markersize=10)

        #Plot tidal bulge locations if the coordinates of the initial and final location differ
        if (self.sc.polarwander_coordinates['thetaTInitial'] != self.sc.polarwander_coordinates['thetaTFinal']
            or self.sc.polarwander_coordinates['phiTInitial'] != self.sc.polarwander_coordinates['phiTFinal']
            or self.sc.polarwander_coordinates['Locked']):

            self.basemap_ax.plot(self.sc.polarwander_coordinates['phiTInitial'],
                self.sc.polarwander_coordinates['thetaTInitial'],
                'ws', markersize=10)
            if self.sc.polarwander_coordinates['phiTInitial'] >=0:
                self.basemap_ax.plot(self.sc.polarwander_coordinates['phiTInitial'] - 180,
                    0 - self.sc.polarwander_coordinates['thetaTInitial'],
                    'ws', markersize=10)
            else:
                self.basemap_ax.plot(self.sc.polarwander_coordinates['phiTInitial'] + 180,
                    0- self.sc.polarwander_coordinates['thetaTInitial'],
                    'ws', markersize=10)

            self.basemap_ax.plot(self.sc.polarwander_coordinates['phiTFinal'],
                self.sc.polarwander_coordinates['thetaTFinal'],
                'ks', markersize=10)
            if self.sc.polarwander_coordinates['phiTFinal'] >=0:
                self.basemap_ax.plot(self.sc.polarwander_coordinates['phiTFinal'] - 180,
                    0 - self.sc.polarwander_coordinates['thetaTFinal'],
                    'ks', markersize=10)
            else:
                self.basemap_ax.plot(self.sc.polarwander_coordinates['phiTFinal'] + 180,
                    0 - self.sc.polarwander_coordinates['thetaTFinal'],
                    'ks', markersize=10)


    def basemap_parameters(self, proj):
        p = { 'projection': proj }
        if proj in ['cyl', 'mill', 'merc']:
            if proj == 'merc':
                if self.grid.lat_min <= -89.9:
                    self.grid.lat_min = -89.9
                if self.grid.lat_max >= 89.9:
                    self.grid.lat_max = 89.9
            p.update({
                'llcrnrlon': self.grid.lon_min,
                'llcrnrlat': self.grid.lat_min,
                'urcrnrlon': self.grid.lon_max,
                'urcrnrlat': self.grid.lat_max})
        elif proj == 'ortho':
            p.update({
                'lat_0': int(round((self.grid.lat_min+self.grid.lat_max)/2)),
                'lon_0': int(round((self.grid.lon_min+self.grid.lon_max)/2))})
        else:
            p.update({'boundinglat': 0,
                'lat_0': (self.grid.lat_min+self.grid.lat_max)/2,
                'lon_0': (self.grid.lon_min+self.grid.lon_max)/2})
        return p

    def get_basemap_ax(self):
        ax = self.get_axes()
        ax.clear()
        p = self.basemap_parameters(self.sc.parameters['projection'])
        p.update({'resolution': None, 'ax': ax})
        self.sc.parameters['ax'] = ax
        basemap_ax = basemap.Basemap(**p)
        return basemap_ax

    def draw_coords(self):
        # Draw a grid onto the plot -- independent of actual grid tab
        coord_lons  = numpy.arange(
        numpy.radians(self.grid.lon_min),
        numpy.radians(self.grid.lon_max),
        numpy.radians(self.step))
        coord_lons = numpy.resize(coord_lons, coord_lons.size + 1)
        coord_lons.put(coord_lons.size - 1, numpy.radians(self.grid.lon_max))
        coord_lats  = numpy.arange(
            numpy.radians(self.grid.lat_min),
            numpy.radians(self.grid.lat_max),
            numpy.radians(self.step))
        coord_lats = numpy.resize(coord_lats, coord_lats.size + 1)
        coord_lats.put(coord_lats.size - 1, numpy.radians(self.grid.lat_max))
        parallel_labels = [1,0,0,1]
        parallel_xoffset = 0
        self.meridians = self.basemap_ax.drawmeridians(numpy.around(numpy.degrees(coord_lons)),
            labels=[1,0,0,1], linewidth=0.5, color='gray', yoffset=5)
        self.parallels = self.basemap_ax.drawparallels(numpy.around(numpy.degrees(coord_lats)),
            labels=parallel_labels, linewidth=0.5, color='gray', xoffset=parallel_xoffset)
        self.basemap_ax.drawmapboundary()

    def adjust_coord_step(self, step):
        """
        Change tick step of coordinate axes.
        """
        self.step = step
        self.save_step()
        #ax = self.get_axes()
        def clear(a):
            for ll, tt in a:
                map(self.ax.lines.remove, ll)
                map(self.ax.texts.remove, tt)
        clear(self.meridians.values())
        clear(self.parallels.values())
        self.plot()

class KPaFormatter(matplotlib.ticker.Formatter):
    def __call__(self, x, pos):
        return "%.f kPa" % (x/1000)

class ScalarPlotPanel(PlotPanel):
    """
    Defines the plot panel of the GUI in terms of PlotPanel, which is in term of SatPanel
    """
    step_field = 'SCALAR_PLOT_STEP'

    def __init__(self, *args, **kw):
        super(ScalarPlotPanel, self).__init__(*args, **kw)

        #self.orbit_hidden = self.nsr_hidden = self.polar_hidden = False
        self.orbit_hidden = self.nsr_hidden = False

        main_sz = wx.BoxSizer(orient=wx.VERTICAL)

        main_sz.Add(self.head_text(), flag=wx.EXPAND|wx.ALL)
        main_sz.AddSpacer(5)
        main_sz.Add(self.plot_sizer(), flag=wx.EXPAND|wx.ALL)
        main_sz.AddSpacer(5)

        main_sz.Add(self.lineaments_sizer())
        main_sz.AddSpacer(5)
        main_sz.Add(wx.StaticLine(self), 0, wx.ALL|wx.EXPAND, 5)
        main_sz.AddSpacer(5)
        main_sz.Add(self.cycloids_sizer())

        self.SetSizer(main_sz)
        self.Fit()
        self.update_parameters()
        self.bind_parameters()

    def head_text(self):
        return WrapStaticText(self,
            label=u"Display a rasterized scalar stress field defined by calculation on " +\
            u"satellite and grid parameters at the resolution defined by grid.  " +\
            u"Tension is positive\n " +\
            u"White circles represent initial rotational poles, black circles are final rotational poles." +\
            u"White squares are initial sub- and anti-jove points, black square are final points." +\
            u"Black triangles are cycloids that could not be initiated, white triangles are cycloids that were initiated but unpropagated.")

    def plot_sizer(self):
        self.plot_fields = {}
        self.plot_vectors = {}
        self.n_interp = 10
        self.tick_formatter = KPaFormatter()
        s = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.scp = self.stress_plot_panel()
        self.init_orbit_slider()
        self.init_nsr_slider()
        #self.init_polar_slider()
        p = self.parameters_sizer()
        s.Add(self.scp, flag=wx.ALL|wx.EXPAND)
        s.AddSpacer(10)
        s.Add(p, flag=wx.ALL|wx.EXPAND)
        return s

    def stress_plot_panel(self):
        scp = StressPlotPanel(self)
        scp.canvas.callbacks.connect('motion_notify_event', self.on_move_in_plot)
        scp.orbit_slider.on_changed(self.on_orbit_updated)
        scp.nsr_slider.on_changed(self.on_nsr_updated)
        #scp.polar_slider.on_changed(self.on_polar_updated)
        scp.save_orbit_series = self.save_orbit_series
        scp.save_nsr_series = self.save_nsr_series
        #scp.save_polar_series = self.save_polar_series
        self.orbit_pos = self.sc.get_parameter(int, 'ORBIT_MIN', 0)
        self.nsr_pos = self.sc.get_parameter(float, 'TIME_MIN', 0)
        #self.polar_pos = self.sc.get_parameter(float,'TIME_MIN',0)
        self.updating = False
        scp.Fit()
        return scp

    # sets up the controls and cells to the right of plot in PlotPanel
    def parameters_sizer(self):
        lp = wx.BoxSizer(orient=wx.VERTICAL)

        # layout as two vertical columns (not sure about row parameter)
        spp1 = wx.FlexGridSizer(rows=1, cols=2)

        # Adds widget controlling projection type
        self.add_projection(spp1)
        # Adds tick mar increment widget
        self.add_stepspin(spp1)
        # Adds plot direction widget
        self.add_direction(spp1)
        # Adds blank space
        spp1.AddSpacer(10)
        spp1.AddSpacer(10)

        # Adds stress range (upper/lower bound included) widget
        self.add_scalectrl(spp1)

        spp1.AddSpacer(15)
        spp1.AddSpacer(15)

        # not sure what this does, but is necessary for plotting
        self.add_stress_field(spp1)

        spp1.Add(wx.StaticText(self, label=u'Plot stresses'), flag=wx.ALIGN_TOP)
        spp2 = wx.FlexGridSizer(rows=9, cols=1)
        # Adds set of radiobuttoms
        self.add_to_plot_stresses(spp2)
        spp1.Add(spp2)

        self.scp.plot_scale(self.scale(), "%.f kPa")

        self.ax = self.scp.get_axes()

        spp1.AddSpacer(15)
        spp1.AddSpacer(15)
        # adds widget displaying long, lat, and stress at cursor
        self.add_value_display(spp1)

        lp.Add(spp1)
        lp.AddSpacer(15)

        self.pw_marker_box = wx.CheckBox(self, label='Show Polar Wander coordinates')
        self.pw_marker_box.SetValue(True)
        lp.Add(self.pw_marker_box, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL)
        self.pw_marker_box.Bind(wx.EVT_CHECKBOX, self.show_pw_markers)

        lp.Add(wx.StaticLine(self), 0, wx.ALL|wx.EXPAND, 5)

        return lp


    def show_pw_markers(self, evt):
        if self.pw_marker_box.GetValue():
            self.sc.parameters['to_plot_pw_markers'] = True
            self.plot()
        else:
            self.sc.parameters['to_plot_pw_markers'] = False
            self.plot()

    def update_parameters(self):
        self.show_needed_sliders()
        super(ScalarPlotPanel, self).update_parameters()

    def scale_spin(self, k):
        self.load_scale(k)
        if k < 0 and self.lbound is None:
            self.lbound = -100
        elif k > 0 and self.ubound is None:
            self.ubound = 100
        ctrl = wx.SpinCtrl(self, min=-100000000, max=100000000)
        if k < 0:
            ctrl.SetValue(int(self.lbound))
        else:
            ctrl.SetValue(int(self.ubound))
        ctrl.Bind(wx.EVT_SPINCTRL, self.select_scale)
        return ctrl

    #@into_hbox
    def add_projection(self, sizer):
        self.sc.parameters['projection'] = 'cyl'
        self.parameters.update(add_combobox2_to_sizer(self, sizer, 'projection', u'Display projection',
            [('cyl', u'Cylindrical Equidistant'),
            ('mill', u'Miller Cylindrical'),
            ('merc', u'Mercator'),
            ('ortho', u'Orthographic'),
            ('npaeqd', u'North-Polar'),
            ('spaeqd', u'South-Polar')]))

    def add_direction(self, sizer):
        self.sc.parameters['direction'] = 'east'
        self.parameters.update(add_radiobox2_to_sizer(self, sizer, 'direction', u'Plot direction',
            [('east', u'East Positive'), ('west', u'West Positive')]))

    #@into_hbox
    # function for adding color scalebar/legend of stress plot
    def add_scalectrl(self, sizer):
        sizer.Add(wx.StaticText(self, label=u"Stress range:"), flag=wx.ALIGN_CENTER_VERTICAL)
        sizer.AddSpacer(15)
        sizer.Add(wx.StaticText(self, label=u"Lower bound [kPa]"), flag=wx.ALIGN_CENTER_VERTICAL)
        self.lbound_ctrl = self.scale_spin(-1)
        sizer.Add(self.lbound_ctrl, flag=wx.ALL|wx.EXPAND)
        sizer.Add(wx.StaticText(self, label=u"Upper bound [kPa]"), flag=wx.ALIGN_CENTER_VERTICAL)
        self.ubound_ctrl = self.scale_spin(1)
        sizer.Add(self.ubound_ctrl, flag=wx.ALL|wx.EXPAND)

    def add_stress_field(self, sizer):
        self.sc.parameters['field'] = 'tens'
        self.parameters.update(add_combobox2_to_sizer(self, sizer, 'field', u'Plot gradient',
            [('tens', u'σ1'),
            ('comp', u'σ3'),
            ('mean', u'(σ1 + σ3)/2'),
            ('diff', u'σ1 - σ3'),
            (None, u'None')]))

    def add_to_plot_stresses(self, sizer):
        self.sc.parameters['to_plot_principal_vectors'] = True
        self.sc.parameters['to_plot_shear_vectors'] = \
            self.sc.parameters['to_plot_longitude_vectors'] = \
            self.sc.parameters['to_plot_latitude_vectors'] = False
        self.parameters.update(add_checkboxes_to_sizer(self, sizer,
            [('to_plot_principal_vectors', u'principal'),
            ('to_plot_latitude_vectors', u'latitude'),
            ('to_plot_longitude_vectors', u'longitude'),
            ('to_plot_shear_vectors', u'shear')]))

    def add_value_display(self, sizer):
        self.val_p = add_parameters_to_sizer(self, sizer,
            [ ('LAT', u'Latitude:'), ('LON', u'Longitude:'),('VAL', u'Stress [kPa]:')])
        for p in ['LAT', 'LON', 'VAL']:
            self.val_p[p].SetEditable(False)

    ###########################################################################
    # Plot Tab Load/Save buttons for lineament and helper functions
    def load_save_buttons(self):
        """
        creates and bind the buttons for loading and saving files
        """
        gridSizer = wx.FlexGridSizer(rows=2, cols=2, hgap=15, vgap=5)

        # create and bind buttons
        shapeLoad = wx.Button(self, label=u'Load from shape file')
        shapeLoad.Bind(wx.EVT_BUTTON, self.on_load_shape)

        shapeSave = wx.Button(self, label=u'Save as shape file')
        shapeSave.Bind(wx.EVT_BUTTON, self.on_save_shape)

        netLoad = wx.Button(self, label=u'Load fom NetCDF file')
        netLoad.Bind(wx.EVT_BUTTON, self.on_load_netcdf)

        netSave = wx.Button(self, label=u'Save as NetCDF file')
        netSave.Bind(wx.EVT_BUTTON, self.on_save_netcdf)

        # add widgets to grid
        gridSizer.AddMany([
            (shapeLoad, 0, wx.ALIGN_CENTER|wx.EXPAND),
            (shapeSave, 0, wx.ALIGN_CENTER|wx.EXPAND),
            (netLoad, 0, wx.ALIGN_CENTER|wx.EXPAND),
            (netSave, 0, wx.ALIGN_CENTER| wx.EXPAND)])

        return gridSizer

    def on_load_shape(self, evt):
        try:
            file_dialog(self,
                message = u"Load from shape file",
                style = wx.OPEN,
                wildcard = 'Shape files (*.shp)|*.shp',
                action = self.load_shape)
        except Exception, e:
            error_dialog(self, str(e), u'Shape Load Error')

    def load_shape(self, filename):
        # walk around char const * restriction
        sf = os.path.splitext(str(filename))[0] + '.shp'
        self.loaded['data'] = shp2lins(sf, stresscalc=self.calc)
        self.loaded['lines'] = []
        d = wx.ColourDialog(self, self.loaded['color'])
        if (d.ShowModal() == wx.ID_OK):
            self.loaded['color'] = d.GetColourData()
        self.plot()

    def on_save_shape(self, evt):
        file_dialog(self,
            message = u"Save to shape file",
            style = wx.SAVE | wx.OVERWRITE_PROMPT,
            wildcard = 'Shape files (*.shp)|*.shp',
            defaultFile = 'lineaments.shp',
            action = self.save_shape)

    def save_shape(self, filename):
        lins2shp(self.loaded['data'] + self.generated['data'], filename)

    def on_load_netcdf(self, evt):
        try:
            file_dialog(self,
                message=u"Load from NetCDF file",
                style=wx.OPEN,
                wildcard=u'NetCDF files (*.nc)|*.nc',
                action=self.load_netcdf)
        except LocalError, e:
            error_dialog(self, str(e), e.title)

    def load_netcdf(self, filename):
        self.sc.load_netcdf(filename)
        self.update_parameters()
        self.plot()

    def on_save_netcdf(self, evt):
        try:
            file_dialog(self,
                message=u"Save to NetCDF file",
                style=wx.SAVE | wx.OVERWRITE_PROMPT,
                defaultFile='gridcalc.nc',
                wildcard=u'NetCDF files (*.nc)|*.nc',
                action=self.sc.save_netcdf)
        except LocalError, e:
            error_dialog(self, str(e), e.title)

    ###########################################################################
    # Defining lineament controls and related functions
    def lineaments_sizer(self):
        """
        Defines sizer for controls for lineament plotting
        """
        # define vars
        #self.lin_p = {}  not used anywhere else, could remove
        self.l_count = 2
        self.generated = { 'data': [], 'color': wx.ColourData(), 'lines': [] }
        self.loaded = { 'data': [], 'color': wx.ColourData(), 'lines': [] }
        self.first_run = True
        self.sc.parameters['to_plot_lineaments'] = True

        # define sizers
        lins = wx.BoxSizer(wx.HORIZONTAL)
        lins_ckSizer = wx.BoxSizer(wx.HORIZONTAL)

        # setup widgets
        self.plot_lins = wx.CheckBox(self, label='Show ')
        self.plot_lins.Bind(wx.EVT_CHECKBOX, self.generate_lins)

        self.l_count_tc = wx.TextCtrl(self, style=wx.TE_PROCESS_ENTER)
        self.l_count_tc.SetValue(str(self.l_count))
        self.l_count_tc.Bind(wx.EVT_TEXT, self.generate_lins)

        # construct ckSizer
        lins_ckSizer.AddSpacer(10)
        lins_ckSizer.Add(self.plot_lins, 0, 20)
        lins_ckSizer.Add(self.l_count_tc, wx.ALL|wx.ALIGN_CENTER_VERTICAL)
        lins_ckSizer.Add(wx.StaticText(self, label=u" Lineaments"), wx.ALL|wx.ALIGN_CENTER_VERTICAL)

        # add checkbox
        lins.Add(lins_ckSizer)
        lins.AddSpacer(15)
        # add buttons
        lins.Add(self.load_save_buttons(), wx.ALL|wx.ALIGN_RIGHT)

        return lins

    def generate_lins(self, evt):
        print 'generate_lins'
        try:
            if self.plot_lins.GetValue():     # plot only if box is checked
                self.l_count = int(self.l_count_tc.GetValue())
            else:
                self.l_count = 0
            self.first_run = False
            b = wx.BusyInfo(u"Performing calculations. Please wait.", self)
            wx.SafeYield()
            self.generated['data'] = self.lingen(self.l_count)
            self.generated['lines'] = []
            del b
            self.plot()
        except:
            self.l_count_tc.SetValue(str(self.l_count))
        self.plot_lineaments()
        print 'end generate_lins'

    def lingen(self, number):
        print 'lingen'
        ll = []
        for lat in numpy.linspace(0, numpy.radians(90), number+2)[1:-1]:
            ll.append(lingen_nsr(self.calc, init_lon=0, init_lat=lat))
        ll += [Lineament(lons=l.lons, lats=-l.lats, stresscalc=l.stresscalc) for l in ll]
        ll += [Lineament(lons=l.lons+satstress.physcon.pi, lats=l.lats, stresscalc=l.stresscalc) for l in ll]
        print ll
        return ll

    def plot_lineaments(self):
        for l in [self.generated, self.loaded]:
            if l['data']:
                l['lines'] = plotlinmap(l['data'], map=self.basemap_ax, color=self.mpl_color(l['color'].GetColour()))[0]

    ###########################################################################
    # Plot Tab Load/Save buttons for cycloids and helper functions
    # Replicated from the load_save_buttons function which is used for lineaments.
    # Created by Peter Sinclair (2016)
    def load_save_buttons_cycloids(self):
        """
        creates and bind the buttons for loading and saving files
        """
        gridSizer = wx.FlexGridSizer(rows=2, cols=2, hgap=15, vgap=5)

        # create and bind buttons
        shapeLoad = wx.Button(self, label=u'Load from shape file')
        shapeLoad.Bind(wx.EVT_BUTTON, self.on_load_shape_cycloid)
        shapeSave = wx.Button(self, label=u'Save as shape file')
        shapeSave.Bind(wx.EVT_BUTTON, self.on_save_shape_cycloid)
        netLoad = wx.Button(self, label=u'Load fom NetCDF file')
        netLoad.Bind(wx.EVT_BUTTON, self.on_load_netcdf_cycloid)
        netSave = wx.Button(self, label=u'Save as NetCDF file')
        netSave.Bind(wx.EVT_BUTTON, self.on_save_netcdf_cycloid)

        # add widgets to grid
        gridSizer.AddMany([
            (shapeLoad, 0, wx.ALIGN_CENTER|wx.EXPAND),
            (shapeSave, 0, wx.ALIGN_CENTER|wx.EXPAND),
            (netLoad, 0, wx.ALIGN_CENTER|wx.EXPAND),
            (netSave, 0, wx.ALIGN_CENTER| wx.EXPAND)])

        return gridSizer

    def on_load_shape_cycloid(self, evt):
        try:
            file_dialog(self,
                message = u"Load from shape file",
                style = wx.OPEN,
                wildcard = 'Shape files (*.shp)|*.shp',
                action = self.load_shape_cycloid)
        except Exception, e:
            error_dialog(self, str(e), u'Shape Load Error')

    def load_shape_cycloid(self, filename):
        # walk around char const * restriction
        sf = os.path.splitext(str(filename))[0] + '.shp'
        self.loaded['data'] = shp2lins(sf, stresscalc=self.calc)
        self.loaded['lines'] = []
        d = wx.ColourDialog(self, self.loaded['color'])
        if (d.ShowModal() == wx.ID_OK):
            self.loaded['color'] = d.GetColourData()
        self.plot()

    def on_save_shape_cycloid(self, evt):
        file_dialog(self,
            message = u"Save to shape file",
            style = wx.SAVE | wx.OVERWRITE_PROMPT,
            wildcard = 'Shape files (*.shp)|*.shp',
            defaultFile = 'cycloids.shp',
            action = SaveCycloidAsShape)

    def on_load_netcdf_cycloid(self, evt):
        try:
            file_dialog(self,
                message=u"Load from NetCDF file",
                style=wx.OPEN,
                wildcard=u'NetCDF files (*.nc)|*.nc',
                action=self.load_netcdf_cycloid)
        except LocalError, e:
            error_dialog(self, str(e), e.title)

    def load_netcdf_cycloid(self, filename):
        self.sc.load_netcdf(filename)
        self.update_parameters()
        self.plot()

    def on_save_netcdf_cycloid(self, evt):
        try:
            file_dialog(self,
                message=u"Save to NetCDF file",
                style=wx.SAVE | wx.OVERWRITE_PROMPT,
                defaultFile='gridcalc.nc',
                wildcard=u'NetCDF files (*.nc)|*.nc',
                action=self.sc.save_netcdf)
        except LocalError, e:
            error_dialog(self, str(e), e.title)

    ###########################################################################
    # Defining cycloid controls and related functions
    def cycloids_sizer(self):
        """
        Defines sizer containing controls for cycloid plotting
        """
        self.cycl_generated = { 'cycdata': [], 'color': wx.ColourData(), 'arcs': [] }
        self.cycl_loaded = { 'cycdata': [], 'color': wx.ColourData(), 'arcs': [] }
        self.first_run = True   # for lineaments

        # create sizers
        cycl = wx.BoxSizer(wx.HORIZONTAL)
        ckSizer = wx.BoxSizer(wx.VERTICAL)

        self.plot_cycl = wx.CheckBox(self, label='Show Cycloids')
        # wrap in sizer
        ckSizer.Add(self.plot_cycl, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL)
        # bind to event
        self.plot_cycl.Bind(wx.EVT_CHECKBOX, self.generate_cycl)

        # A checkbox to plot triangles at the cycloid's location if they are unable to start or propagate.  -PS 2016
        self.plot_triangles = wx.CheckBox(self, label='Plot marker if unable to create cycloid')
        ckSizer.Add(self.plot_triangles, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL)
        self.plot_triangles.Bind(wx.EVT_CHECKBOX, self.generate_cycloid_markers)
        self.plot_triangles.SetValue(True) # Starts enabled.

        self.cycl_names_cb = wx.CheckBox(self, label='Show cycloid names')
        ckSizer.Add(self.cycl_names_cb , 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL)
        self.cycl_names_cb.Bind(wx.EVT_CHECKBOX, self.plot_cycl_names)
        self.cycl_names_cb.SetValue(False)

        saveMany = wx.Button(self, label="Save Multiple Cycloids")
        saveMany.Bind(wx.EVT_BUTTON, self.save_many_cycloids)
        ckSizer.AddSpacer(5)
        ckSizer.Add(saveMany)

        # add grid to sizer
        cycl.AddSpacer(10)
        cycl.Add(ckSizer, wx.ALL|wx.ALIGN_LEFT)
        cycl.AddSpacer(5)
        cycl.AddSpacer(15)
        cycl.Add(self.load_save_buttons_cycloids(), wx.ALL|wx.ALIGN_RIGHT)
        return cycl

    def generate_cycl(self, evt):
        if self.plot_cycl.GetValue(): # plot only if box is checked
            self.sc.parameters['to_plot_cycloids'] = True
            self.plot()
        else:
            self.sc.parameters['to_plot_cycloids'] = False
            self.plot()

    def generate_cycloid_markers(self, evt):
        if self.plot_triangles.GetValue():
            self.sc.parameters['to_plot_triangles'] = True
            self.plot()
        else:
            self.sc.parameters['to_plot_triangles'] = False
            self.plot()

    def plot_cycl_names(self, evt):
        s = self.cycl_names_cb.GetValue()
        if s:
            self.sc.parameters['show_cycl_names'] = True
            self.plot()
        else:
            self.sc.parameters['show_cycl_names'] = False
            self.plot()

    def plot_cycloids(self):
        if self.sc.parameters['to_plot_many_cycloids']:
            for i, cycloid_params in enumerate(self.sc.params_for_cycloids.items()):

                if not self.sc.cycloids.has_key(i) or self.sc.many_changed or self.sc.cycloid_changed:
                    self.sc.cycloids[i] = Cycloid(self.calc, **cycloid_params[1])
                self.sc.cycloids[i].plotcoordsonbasemap(self.basemap_ax, self.sc.parameters['ax'],self.orbit_pos, self.sc.parameters['to_plot_triangles'], self.sc.parameters['show_cycl_names'])
            self.sc.many_changed = False
        else:
            if (self.sc.cyc == None or self.sc.cycloid_changed):
                self.sc.cyc = Cycloid(self.calc, self.sc.parameters['cycloid_name'], self.sc.parameters['YIELD'], self.sc.parameters['PROPAGATION_STRENGTH'], self.sc.parameters['PROPAGATION_SPEED'], \
                                      self.sc.parameters['STARTING_LATITUDE'], self.sc.parameters['STARTING_LONGITUDE'], self.sc.parameters['STARTING_DIRECTION'], \
                                      self.sc.parameters['VARY_VELOCITY'],self.sc.parameters['k'],self.sc.get_parameter(float, 'ORBIT_MAX', 360), 0.1)
                self.sc.cycloid_changed = False
            self.sc.cyc.plotcoordsonbasemap(self.basemap_ax, self.sc.parameters['ax'], self.orbit_pos, self.sc.parameters['to_plot_triangles'], self.sc.parameters['show_cycl_names'] )

    def save_many_cycloids(self, evt):
        # if a set of parameters from *.csv hasn't been uploaded, treat it like an error
        # with a popup window
        if not self.sc.parameters["to_plot_many_cycloids"]:
            errorMsg = """Please upload a set of cycloid parameters from *.csv file."""
            msg = wx.MessageDialog(self, errorMsg, "No input file found!", wx.OK | wx.ICON_ERROR)
            msg.ShowModal()
            msg.Destroy()
        # otherwise generate and save plots in designated folder
        else:
            chooseFolder = wx.DirDialog(self, "Choose a directory:", style=wx.DD_DEFAULT_STYLE)
            # so that folderName can accessed outside
            folderName = ""
            if chooseFolder.ShowModal() == wx.ID_OK:
                folderName = chooseFolder.GetPath()
            # Blanks out the entire window, which prevents people from changing tabs
            # or doing anything else, which happens naturally anyways.
            # self.Hide()
            i = 0
            while i < len(self.parameters['YIELD']):
                # create cycloid
                threshold = float(self.parameters['YIELD'][i])
                strength = float(self.parameters['PROPAGATION_STRENGTH'][i])
                speed = float(self.parameters['PROPAGATION_SPEED'][i])
                lon = float(self.parameters['STARTING_LONGITUDE'][i])
                lat = float(self.parameters['STARTING_LATITUDE'][i])
                propdir = self.parameters['STARTING_DIRECTION']
                self.sc.cyc.plotcoordsonbasemap(self.calc, self.basemap_ax,
                                    threshold, strength, speed, lon, lat,
                                    propdir,
                                    self.sc.get_parameter(float, 'ORBIT_MAX', 360), self.sc.parameters['to_plot_triangles'])
                # save cycloid
                plotName = str(threshold) + "_" + str(strength) + "_" +  str(speed) + "_" + str(lat) + "_" + str(lon) + "_" + str(propdir)
                self.scp.figure.savefig(folderName + '/' + plotName + ".png", bbox_inches='tight')
                # To have one cycloid saved per image, clear basemap if cycloid was plotted
                if self.ax.lines != []:
                    # self.ax.lines.pop(0)
                    self.ax.lines = []
                i += 1
    ###########################################################################

    def on_orbit_updated(self, val):
        if self.updating:
            return
        self.orbit_pos = self.scp.orbit_slider.val
        self.updating = True
        self.scp.nsr_slider.first()
        self.nsr_pos = 0
        self.updating = False
        self.plot()

    def on_nsr_updated(self, val):
        if self.updating:
            return
        self.nsr_pos = self.scp.nsr_slider.val
        self.updating = True
        self.scp.orbit_slider.first()
        self.orbit_pos = 0
        self.updating = False
        self.plot()

    def on_polar_updated(self, val):
        if self.updating:
            return
        self.polar_pos = self.scp.polar_slider.val
        self.updating = True
        self.scp.orbit_slider.first()
        self.orbit_pos = 0
        self.updating = False
        self.plot()

    def prepare_plot(self):
        b = wx.BusyInfo(u"Performing calculations. Please wait.", self)
        wx.SafeYield()
        self.prepare_plot_series()
        del b

    def get_grid_time(self):
        if self.orbit_pos > self.nsr_pos:
            s = self.sc.get_satellite()
            return self.orbit_pos/360.0*s.orbit_period()
        else:
            return self.nsr_pos*seconds_in_year

    def mk_change_param(self, k):
        def on_change(evt):
            if k == 'direction':
                #handles change of east-west positivity
                # reverse
                temp_min = -self.sc.get_parameter(float, "LON_MAX")
                temp_max = -self.sc.get_parameter(float, "LON_MIN")
                self.sc.set_parameter("LON_MIN", temp_min)
                self.sc.set_parameter("LON_MAX", temp_max)
                self.plot()
            else:
                self.sc.set_parameter(k, self.parameters[k].GetValue())
                self.plot()
        return on_change

    def load_scale(self, k):
        try:
            if k < 0:
                self.lbound = int(config.load('PLOT_LBOUND'))
            else:
                self.ubound = int(config.load('PLOT_UBOUND'))
        except:
            if k < 0:
                self.lbound = None
            else:
                self.ubound = None
        if k < 0:
            return self.lbound
        else:
            return self.ubound

    def save_scale(self):
        config.save(PLOT_LBOUND=self.lbound, PLOT_UBOUND=self.ubound)

    def select_scale(self, evt):
        l = self.lbound_ctrl.GetValue()
        u = self.ubound_ctrl.GetValue()
        try:
            fl = int(l)
            fu = int(u)
            if self.lbound != l or self.ubound != u:
                self.lbound = l
                self.ubound = u
                self.save_scale()
                self.scp.plot_scale(self.scale(), "%.f kPa")
                self.select_color_range(fl*1000, fu*1000)
        except:
            self.lbound_ctrl.SetValue(self.lbound)
            self.ubound_ctrl.SetValue(self.ubound)

    def select_color_range(self, vmin, vmax):
        self.plot()
        self.cb.set_clim(vmin, vmax)
        self.cb.update_bruteforce(self.im)
        self.cb.draw_all()
        self.draw()

    def get_axes(self):
        return self.ax

    def draw(self):
        self.scp.draw()

    def on_move_in_plot(self, evt):
        if evt.inaxes:
            try:
                x,y = self.basemap_ax(evt.xdata, evt.ydata, inverse=True)
                i = int((x - self.grid.lon_min)/(self.grid.lon_max - self.grid.lon_min + 1e-2)*self.grid.lon_num*self.n_interp)
                j = int((y - self.grid.lat_min)/(self.grid.lat_max - self.grid.lat_min + 1e-2)*self.grid.lat_num*self.n_interp)
                x1, y1, plot_field1 = self.plot_fields[self.get_grid_time()][self.sc.parameters['field']]
                self.val_p['LON'].SetValue("%.2f" % x)
                self.val_p['LAT'].SetValue("%.2f" % y)
                self.val_p['VAL'].SetValue("%.2f" % (plot_field1[i,j]/1000.))
            except:
                pass

    def colorbar(self, replot_colorbar):
        try:
            self.cb
            if replot_colorbar:
                self.adjust_to_tight()
                self.scp.figure.delaxes(self.cb.ax)
                self.cb = self.scp.colorbar(self.im, ax=self.ax, format=self.tick_formatter)
        except Exception, e:
            self.adjust_to_tight()
            self.cb = self.scp.colorbar(self.im, ax=self.ax, format=self.tick_formatter)

    def consider_obliq_lons(self, lx, rx):
        if self.sc.parameters.get('Obliquity'):
            if int(round(lx)) % 90 == 0:
                lx += 1
            if int(round(rx)) % 90 == 0:
                rx -= 1
        return lx, rx

    def consider_obliq_lats(self, ly, hy):
        if self.sc.parameters.get('Obliquity'):
            if int(round(ly)) % 90 == 0:
                ly += 1
            if int(round(hy)) % 90 == 0:
                hy -= 1
        return ly, hy

    def consider_lons(self):
        lx = self.grid.lon_min
        rx = self.grid.lon_max
        lx, rx = self.consider_obliq_lons(lx, rx)
        if self.sc.parameters['projection'] == 'ortho' and rx - lx >= 180:
            cx = int(round((lx + rx)/2))
            lx = cx - 90 + 1
            rx = cx + 90 - 1
        return numpy.linspace(lx, rx, self.grid.lon_num*self.n_interp)

    def consider_lats(self):
        ly = self.grid.lat_min
        hy = self.grid.lat_max
        ly, hy = self.consider_obliq_lats(ly, hy)
        proj = self.sc.parameters['projection']
        if proj == 'spaeqd' and hy > 0 and ly < 0:
            hy = 0
        elif proj == 'npaeqd' and hy > 0 and ly < 0:
            ly = 0
        elif proj == 'ortho' and hy - ly >= 180:
            cy = int(round((hy + ly)/2))
            ly = cy - 90 + 1
            hy = cy + 90 - 1
        return numpy.linspace(ly, hy, self.grid.lat_num*self.n_interp)

    def prepare_plot_series(self):
        self.plot_fields.clear()
        self.plot_vectors.clear()
        sat = self.sc.get_satellite()

        lons = self.consider_lons()
        lats  = self.consider_lats()
        phis, thetas = numpy.meshgrid(lons, lats)
        x,y = self.basemap_ax(phis, thetas)
        i,j = numpy.meshgrid(
            numpy.linspace(0, self.grid.lon_num - 1, self.grid.lon_num*self.n_interp),
            numpy.linspace(0, self.grid.lat_num - 1, self.grid.lat_num*self.n_interp))

        self.vector_mesh_lons, self.vector_mesh_lats = self.vector_meshes()

        # monkey patching not to touch library code
        def imshow(plot_field, cmap=None, **kw):
            plot_field1 = scipy.ndimage.map_coordinates(plot_field, [i,j])
            self.plot_fields[self.plot_time][self.plot_field] = (x, y, plot_field1)

        def quiver(x, y, u, v, **kw):
            self.plot_vectors[self.plot_time][self.plot_vector].append((x, y, u, v, kw))

        _imshow = self.basemap_ax.imshow
        self.basemap_ax.imshow = imshow
        _quiver = self.basemap_ax.quiver
        self.basemap_ax.quiver = quiver

        orbit_period = sat.orbit_period()
        o = self.sc.get_parameter(float, 'ORBIT_MIN', 0)
        om = self.sc.get_parameter(float, 'ORBIT_MAX', 0)
        n = self.sc.get_parameter(float, 'ORBIT_NUM', 0)
        if n > 0:
            s = (om - o)/n
            while o <= om:
                self.plot_time = o/360.0*orbit_period
                self.prepare_plot_for_time()
                o += s
        nm = self.sc.get_parameter(float, 'TIME_MIN', 0)
        s = self.sc.get_parameter(float, 'nsr_time', 0)
        n = self.sc.get_parameter(int, 'TIME_NUM', 0)
        for k in range(0, n+1):
            self.plot_time = (s*k + nm)*seconds_in_year
            self.prepare_plot_for_time()
        self.basemap_ax.imshow = _imshow
        self.basemap_ax.quiver = _quiver

    def prepare_plot_for_time(self):
        # we use self.plot_time instead of passing it as parameter 
        # because it is used in redefined imshow and quiver in function above
        self.plot_fields[self.plot_time] = {}
        lon_min, lon_max = self.consider_obliq_lons(self.grid.lon_min,
                self.grid.lon_max)
        lat_min, lat_max = self.consider_obliq_lats(self.grid.lat_min,
                self.grid.lat_max)
        for self.plot_field in ['tens', 'comp', 'mean', 'diff']:
            scalar_grid(
                stresscalc = self.calc,
                nlons = self.grid.lon_num,
                nlats = self.grid.lat_num,
                min_lon = numpy.radians(lon_min),
                max_lon = numpy.radians(lon_max),
                min_lat = numpy.radians(lat_min),
                max_lat = numpy.radians(lat_max),
                time_t = self.plot_time,
                field = self.plot_field,
                basemap_ax = self.basemap_ax)
        # self.plot_vector for same reasons as self.plot_time
        self.plot_vector = 'principal'
        self.plot_vectors[self.plot_time] = { self.plot_vector: [] }
        # Plots principal stresses
        vector_points1(stresscalc=self.calc,
            lons = self.vector_mesh_lons,
            lats = self.vector_mesh_lats,
            time_t = self.plot_time,
            plot_greater = True,
            plot_lesser = True,
            plot_comp = True,
            plot_tens = True,
            scale = self.scale()*vector_mult,
            basemap_ax = self.basemap_ax)
        for self.plot_vector in ['latitude', 'longitude', 'shear']:
            # Plots lat, lon, and shear stresses
            self.plot_vectors[self.plot_time][self.plot_vector] = []
            vector_points2(stresscalc=self.calc,
                lons = self.vector_mesh_lons,
                lats = self.vector_mesh_lats,
                time_t = self.plot_time,
                plot_norm_lat = (self.plot_vector == 'latitude'),
                plot_norm_lon = (self.plot_vector == 'longitude'),
                plot_shear =( self.plot_vector == 'shear'),
                scale = self.scale()*vector_mult,
                basemap_ax = self.basemap_ax)

    def scale(self):
        def max_abs(*v):
            ''' finds the maximum of the absolute values of [vectors?] '''
            # how diff from max(map(abs, v))?
            return max(*map(abs, v))
        return max_abs(self.ubound, self.lbound)

    def plot_gradient(self):
        try:
            x, y, plot_field1 = self.plot_fields[self.get_grid_time()][self.sc.parameters['field']]
            l = int(self.lbound) * 1000
            u = int(self.ubound) * 1000
            self.im = self.basemap_ax.pcolormesh(x, y, numpy.transpose(plot_field1), cmap='gist_rainbow_r', vmin=l, vmax=u)
        except Exception, e:
            print "%s: %s" % (e.__class__.__name__, e)

    def plot_grid_calc(self):
        replot_colorbar = False
        if self.sc.changed() or self.sc.calc_changed:
            self.orbit_pos = self.sc.get_parameter(int, 'ORBIT_MIN', 0)
            self.nsr_pos = self.sc.get_parameter(float, 'TIME_MIN', 0)
            self.hide_sliders()
            self.show_needed_sliders()
            self.prepare_plot()
            self.sc.calc_changed = False
            replot_colorbar = True
        elif self.sc.projection_changed:
            self.prepare_plot()
            self.sc.projection_changed = False

        if self.sc.parameters['field']:
            self.plot_gradient()
        if self.sc.parameters['to_plot_principal_vectors']:
            self.plot_principal_vectors()
        if self.sc.parameters['to_plot_latitude_vectors'] \
        or self.sc.parameters['to_plot_longitude_vectors'] \
        or self.sc.parameters['to_plot_shear_vectors']:
            self.plot_stress_vectors()
        if self.sc.parameters['to_plot_lineaments']:
            self.plot_lineaments()
        if self.sc.parameters['to_plot_cycloids']:
            self.plot_cycloids()

        self.colorbar(replot_colorbar)

    def adjust_to_tight(self):
        [lat0, lat1, lon0, lon1] = map(float, [ self.sc.parameters[x] for x in ['LAT_MIN', 'LAT_MAX', 'LON_MIN', 'LON_MAX']])
        l = (lon1 - lon0)/(lat1 - lat0)*scale_bar_length
        s = (l - scale_bar_length)/2
        #self.scp.figure.subplots_adjust(left=scale_left - s, right=scale_left + scale_bar_length + s + 0.3*l)
        self.scp.figure.subplots_adjust(left = scale_left - s,# - 0.03,
            right = scale_left + scale_bar_length + 1.5*s + 0.1)

    def vector_meshes(self):
        lon_min, lon_max = self.consider_obliq_lons(self.grid.lon_min,
                self.grid.lon_max)
        lat_min, lat_max = self.consider_obliq_lats(self.grid.lat_min,
                self.grid.lat_max)
        vector_grid_lons  = numpy.linspace(
            numpy.radians(lon_min),
            numpy.radians(lon_max),
            self.grid.lon_num)
        vector_grid_lats  = numpy.linspace(
            numpy.radians(lat_min),
            numpy.radians(lat_max),
            self.grid.lat_num)
        vector_mesh_lons, vector_mesh_lats = numpy.meshgrid(vector_grid_lons, vector_grid_lats)
        vector_mesh_lons = numpy.ravel(vector_mesh_lons)
        vector_mesh_lats = numpy.ravel(vector_mesh_lats)
        return vector_mesh_lons, vector_mesh_lats

    def plot_stress_vectors(self):
        if self.sc.parameters['to_plot_latitude_vectors']:
            for x, y, u, v, kw in self.plot_vectors[self.get_grid_time()]['latitude']:
                self.basemap_ax.quiver(x, y, u, v, **kw)
        if self.sc.parameters['to_plot_longitude_vectors']:
            for x, y, u, v, kw in self.plot_vectors[self.get_grid_time()]['longitude']:
                self.basemap_ax.quiver(x, y, u, v, **kw)
        if self.sc.parameters['to_plot_shear_vectors']:
            for x, y, u, v, kw in self.plot_vectors[self.get_grid_time()]['shear']:
                self.basemap_ax.quiver(x, y, u, v, **kw)

    def plot_principal_vectors(self):
        for x, y, u, v, kw in self.plot_vectors[self.get_grid_time()]['principal']:
            kw['scale'] = float(self.scale()*vector_mult)
            self.basemap_ax.quiver(x, y, u, v, **kw)

    def mpl_color(self, color):
        return map(lambda c: float(c)/255, color[0:3])

    def show_needed_sliders(self):
        if self.sc.parameters.get('Nonsynchronous Rotation', False) \
        and self.sc.parameters.get('TIME_MIN') and self.sc.parameters.get('nsr_time') and self.sc.parameters.get('TIME_NUM'):
            self.reveal_nsr_slider()
        else:
            self.hide_nsr_slider()
        if (self.sc.parameters.get('Diurnal', False) or self.sc.parameters.get('Obliquity', False)) \
        and self.sc.parameters.get('ORBIT_MIN') and self.sc.parameters.get('ORBIT_MAX') and self.sc.parameters.get('ORBIT_NUM'):
            self.reveal_orbit_slider()
        else:
            self.hide_orbit_slider()
        """
        # Polar slider is not shown because it is not currently needed.  -PS 2016
        if self.sc.parameters.get('Polar Wander', False):
            self.reveal_polar_slider()
        else:
            self.hide_polar_slider()
        """

    def init_orbit_slider(self):
        self.scp.change_orbit_slider(
            self.sc.get_parameter(float, 'ORBIT_MIN', 0),
            self.sc.get_parameter(float, 'ORBIT_MAX', 1),
            self.sc.get_parameter(float, 'ORBIT_NUM', 10),
            self.orbit_pos)

    def init_nsr_slider(self):
        nm = self.sc.get_parameter(float, 'TIME_MIN', 0)
        self.scp.change_nsr_slider(
            nm,
            nm + self.sc.get_parameter(float, 'nsr_time', 0)*self.sc.get_parameter(float, 'TIME_NUM', 0),
            self.sc.get_parameter(int, 'TIME_NUM', 1),
            self.nsr_pos)

    def init_polar_slider(self):
        nm = self.sc.get_parameter(float, 'TIME_MIN', 0)
        self.scp.change_polar_slider(
           nm,
           nm + self.sc.get_parameter(float, 'nsr_time', 0)*self.sc.get_parameter(float, 'TIME_NUM', 0),
           self.sc.get_parameter(int, 'TIME_NUM', 1),
           self.polar_pos)

    def hide_orbit_slider(self):
        if not self.orbit_hidden:
            self.orbit_hidden = True
            self.scp.del_orbit()

    def hide_nsr_slider(self):
        if not self.nsr_hidden:
            self.nsr_hidden = True
            self.scp.del_nsr()

    def hide_polar_slider(self):
        if not self.polar_hidden:
            self.polar_hidden = True
            self.scp.del_polar()

    def hide_sliders(self):
        self.hide_nsr_slider()
        self.hide_orbit_slider()
        #self.hide_polar_slider()

    def reveal_orbit_slider(self):
        if self.orbit_hidden:
            self.orbit_hidden = False
            self.scp.add_orbit()
            self.init_orbit_slider()
            self.scp.orbit_slider.on_changed(self.on_orbit_updated)
            self.scp.save_orbit_series = self.save_orbit_series

    def reveal_nsr_slider(self):
        if self.nsr_hidden:
            self.nsr_hidden = False
            self.scp.add_nsr()
            self.scp.nsr_slider.on_changed(self.on_nsr_updated)
            self.init_nsr_slider()
            self.scp.save_nsr_series = self.save_nsr_series

    def reveal_polar_slider(self):
        if self.polar_hidden:
            self.polar_hidden = False
            self.scp.add_polar()
            self.scp.polar_slider.on_changed(self.on_polar_updated)
            self.init_polar_slider()
            self.scp.save_polar_series = self.save_polar_series

    def hide_orbit_controls(self):
        self.scp.del_orbit_controls()
        self.scp.orbit_slider.on_changed(lambda v: v)

    def hide_nsr_controls(self):
        self.scp.del_nsr_controls()
        self.scp.nsr_slider.on_changed(lambda v: v)

    def reveal_orbit_controls(self):
        self.scp.add_orbit_controls()
        self.scp.save_orbit_series = self.save_orbit_series
        self.scp.orbit_slider.on_changed(self.on_orbit_updated)

    def reveal_nsr_controls(self):
        self.scp.add_nsr_controls()
        self.scp.save_nsr_series = self.save_nsr_series
        self.scp.nsr_slider.on_changed(self.on_nsr_updated)

    def save_orbit_series(self, dir='.'):
        b = wx.BusyInfo(u"Saving images. Please wait.", self)
        wx.SafeYield()
        old_orbit_pos = self.orbit_pos
        sat = self.sc.get_satellite()
        orbit_period = sat.orbit_period()
        o = self.sc.get_parameter(float, 'ORBIT_MIN', 0)
        om = self.sc.get_parameter(float, 'ORBIT_MAX', 0)
        n = self.sc.get_parameter(float, 'ORBIT_NUM', 0)
        s = (om - o)/n
        self.hide_orbit_controls()

        localtime = time.asctime(time.localtime(time.time()))
        location = dir + "/" + self.sc.parameters['SYSTEM_ID']
        directory = location + " " + localtime
        if os.path.isdir(location):
            os.mkdir(directory)
        else:
            os.mkdir(location)
            os.mkdir(directory)

        while o <= om:
            self.orbit_pos = o
            self.plot_no_draw()
            self.scp.orbit_slider.set_val(self.orbit_pos)
            self.scp.figure.savefig("%s/orbit_%03d.%02d.png" %
                (directory, int(self.orbit_pos), round(100.*(self.orbit_pos - int(self.orbit_pos)))),
                bbox_inches='tight', pad_inches=1.5)
            o += s
        self.orbit_pos = old_orbit_pos
        self.reveal_orbit_controls()
        self.init_orbit_slider()
        self.scp.orbit_slider.set_val(self.orbit_pos)
        self.plot()
        del b

    def save_nsr_series(self, dir='.'):
        b = wx.BusyInfo(u"Saving images. Please wait.", self)
        wx.SafeYield()
        old_nsr_pos = self.nsr_pos
        nm = self.sc.get_parameter(float, 'TIME_MIN', 0)
        s = self.sc.get_parameter(float, 'nsr_time', 0)
        n = self.sc.get_parameter(int, 'TIME_NUM', 0)
        self.hide_nsr_controls()

        localtime = time.asctime(time.localtime(time.time()))
        location = dir + "/" + self.sc.parameters['SYSTEM_ID']
        directory = location + "/" + localtime
        if os.path.isdir(location):
            os.mkdir(directory)
        else:
            os.mkdir(location)
            os.mkdir(directory)

        for k in range(0, n+1):
            self.nsr_pos = nm + s*k
            self.scp.nsr_slider.set_val(self.nsr_pos)
            self.plot_no_draw()
            self.scp.figure.savefig("%s/nsr_%03d.png" % (directory, k), bbox_inches='tight', pad_inches=0.5)
        self.nsr_pos = old_nsr_pos
        self.reveal_nsr_controls()
        self.init_nsr_slider()
        self.scp.nsr_slider.set_val(self.nsr_pos)
        self.plot()
        del b

    def save_polar_series(self, dir='.'):
        b = wx.BusyInfo(u"Saving images. Please wait.", self)
        wx.SafeYield()
        old_polar_pos = self.polar_pos
        nm = self.sc.get_parameter(float, 'TIME_MIN', 0)
        s = self.sc.get_parameter(float, 'polar_time', 0)
        n = self.sc.get_parameter(int, 'TIME_NUM', 0)
        self.hide_polar_controls()
        for k in range(0, n+1):
            self.polar_pos = nm + s*k
            self.scp.polar_slider.set_val(self.polar_pos)
            self.plot_no_draw()
            self.scp.figure.savefig("%s/polar_%03d.png" % (dir, k), bbox_inches='tight', pad_inches=0.5)
        self.polar_pos = old_polar_pos
        self.reveal_polar_controls()
        self.init_polar_slider()
        self.scp.polar_slider.set_val(self.polar_pos)
        self.plot()
        del b


# ===============================================================================
# PANEL CONTAINING ALL TABS
# ===============================================================================
class SatStressPanel(wx.Panel):
    """
    Defines the panel that contains all GUI pages
    """
    def __init__(self, *args, **kw):
        wx.Panel.__init__(self, *args, **kw)

        self.SetMinSize((1024, 640))
        sz = wx.BoxSizer(orient=wx.VERTICAL)

        self.nb = wx.Notebook(self)

        self.sc = SatelliteCalculation()
        slp = SatelliteLayersPanel(self.nb, satellite_calculation=self.sc)
        stp = StressListPanel(self.nb, satellite_calculation=self.sc)
        self.tp = PointPanel(self.nb, satellite_calculation=self.sc)
        gp = GridCalcPanel(self.nb, satellite_calculation=self.sc)
        self.spp = ScalarPlotPanel(self.nb, satellite_calculation=self.sc)
        self.cy = CycloidsPanel(self.nb, satellite_calculation=self.sc)

        # Assign each panel to a page and give it a name
        self.nb.AddPage(slp, u"Satellite")
        self.nb.AddPage(stp, u"Stresses")
        self.nb.AddPage(self.tp, u"Point")
        self.nb.AddPage(gp, u"Grid")
        self.nb.AddPage(self.cy, u"Cycloids")
        self.nb.AddPage(self.spp, u"Plot")
        # self.nb.AddPage(dummy, u'Test')

        sz.Add(self.nb, 1, wx.ALL|wx.EXPAND)
        #sz.Add(bp, 0, wx.ALIGN_BOTTOM | wx.EXPAND)

        self.SetSizer(sz)
        self.Fit()
        self.sc.parameters['show_cycl_names'] = False
        wx.EVT_NOTEBOOK_PAGE_CHANGED(self, self.nb.GetId(), self.page_change)

    def page_change(self, evt):
        p = self.nb.GetCurrentPage()
        if isinstance(p, SatPanel):
            p.update_parameters()
        if isinstance(p, PlotPanel):
            p.plot()

class SatStressFrame(wx.Frame):
    """
    Actually holds all the tabs? Wrapper for Panel that holds everythings
    """
    def __init__(self, parent, *args, **kw):
        wx.Frame.__init__(self, parent, *args, **kw)
        self.SetSizer(wx.BoxSizer(wx.VERTICAL))
        self.p = SatStressPanel(self)
        self.GetSizer().Add(self.p, 1, wx.ALL|wx.EXPAND, 10)

        menubar = wx.MenuBar()

        ##### 'File' option of menubar #####
        File = wx.Menu()
        export = File.Append(wx.ID_SAVE, '&Export\tCtrl+S', 'Save all variables')
        self.Bind(wx.EVT_MENU,self.onExport, export)
        load = File.Append(wx.ID_OPEN, '&Load\tCtrl+O', 'Load a set of variables')
        self.Bind(wx.EVT_MENU, self.onLoad, load)
        quit = File.Append(wx.ID_ANY, '&Quit\tCtrl+Q', 'Quit Application')
        self.Bind(wx.EVT_MENU, self.onQuit, quit)

        menubar.Append(File,"File")

        ##### 'Information' option of menubar #####        
        Information = wx.Menu()

        About = wx.Menu()
        rights = About.Append(wx.ID_ANY, '&Copyright')
        self.Bind(wx.EVT_MENU, self.onRights, rights)
        updates = About.Append(wx.ID_ANY, '&Version')
        self.Bind(wx.EVT_MENU, self.onUpdates, updates)
        contact = About.Append(wx.ID_ANY, '&Contact')
        self.Bind(wx.EVT_MENU, self.onContacts, contact)
        develop = About.Append(wx.ID_ANY, '&Development')
        self.Bind(wx.EVT_MENU, self.onDevelopment, develop)

        Information.AppendMenu(wx.ID_ANY, "&About", About)
        Information.AppendSeparator()

        References = wx.Menu()
        Diurnalref = References.Append(wx.ID_ANY, '&Diurnal')
        self.Bind(wx.EVT_MENU, self.onDiurnalref, Diurnalref)
        NSRref = References.Append(wx.ID_ANY, '&Nonsynchronous Rotation')
        self.Bind(wx.EVT_MENU, self.onNSRref, NSRref)
        Obliquityref = References.Append(wx.ID_ANY, '&Obliquity')
        self.Bind(wx.EVT_MENU, self.onObliquityref, Obliquityref)
        ISTref = References.Append(wx.ID_ANY, '&Ice Shell Thickening')
        self.Bind(wx.EVT_MENU, self.onISTref, ISTref)
        PWref = References.Append(wx.ID_ANY, '&Polar Wander')
        self.Bind(wx.EVT_MENU, self.onPWref, PWref)
        Cycloidsref = References.Append(wx.ID_ANY, '&Cycloids')
        self.Bind(wx.EVT_MENU, self.onCycloidsref, Cycloidsref)

        Information.AppendMenu(wx.ID_ANY, "&References", References)

        menubar.Append(Information, "&Information")

        ##### 'Help' option of menubar ######
        Help = wx.Menu()
        Tutorial = Help.Append(wx.ID_ANY, '&Getting Started\tf1')
        self.Bind(wx.EVT_MENU, self.onTutorial, Tutorial)
        HelpSat = Help.Append(wx.ID_ANY, '&Satellite Tab')
        self.Bind(wx.EVT_MENU, self.onHelpSat, HelpSat)
        HelpStress = Help.Append(wx.ID_ANY, '&Stresses Tab')
        self.Bind(wx.EVT_MENU, self.onHelpStresses, HelpStress)
        HelpPoint = Help.Append(wx.ID_ANY, '&Point Tab')
        self.Bind(wx.EVT_MENU, self.onHelpPoint, HelpPoint)
        HelpGrid = Help.Append(wx.ID_ANY, '&Grid Tab')
        self.Bind(wx.EVT_MENU, self.onHelpGrid, HelpGrid)
        HelpCycloids = Help.Append(wx.ID_ANY, '&Cycloids Tab')
        self.Bind(wx.EVT_MENU, self.onHelpCycloids, HelpCycloids)
        HelpPlot = Help.Append(wx.ID_ANY, '&Plot Tab')
        self.Bind(wx.EVT_MENU, self.onHelpPlot, HelpPlot)
        menubar.Append(Help, "&Help")

        self.SetMenuBar(menubar)

        exit_id = wx.NewId()
        wx.EVT_MENU(self, exit_id, self.exit)
        accel = wx.AcceleratorTable([
            (wx.ACCEL_CTRL, ord('W'), exit_id)])
        self.SetAcceleratorTable(accel)

        # Bind our events from the close dialog 'x' on the frame
        self.Bind(wx.EVT_CLOSE, self.OnCloseFrame)

        # SetSizeHints(minW, minH, maxW, maxH)
        # This function effectively enforces a lower bound to SatStressGUI window resizing.
        # To allow for unrestricted window resizing, simply remove this line.
        self.SetSizeHints(1045,710,2000, 2000)

        self.Fit()
        self.Show(True)
        self.CenterOnScreen()
        self.p.SetFocus()

    def onExport(self,evt):
        try:
            file_dialog(self,
                    message=u"Save configuration",
                    style=wx.SAVE | wx.OVERWRITE_PROMPT,
                    wildcard='Satstress files (*.sats)|*.sats',
                    action=self.saveFile)
        except Exception, e:
            error_dialog(self, str(e), u'Save Error')

    def onLoad(self,evt):
        try:
            file_dialog(self,
                        message=u"Load configuration",
                        style=wx.OPEN,
                        wildcard='Satstress files (*.sats)|*.sats',
                        action=self.loadFile)
        except Exception, e:
            error_dialog(self, str(e), u'Load Error')

    def loadFile(self,filename):
        f = open(filename)

        for k,v in nvf2dict(f).items():
            if k == 'point_rows':
                self.p.tp.set_num_rows(float(v))
            if str(v)[0] == '[':  #Load in a list
                l = eval(v)
                for i in range(1, len(l)):
                    self.p.sc.set_parameter(k, l[i], point = i)
            else:
                self.p.sc.set_parameter(k,v)

        self.p.sc.grid_changed = True
        self.p.sc.nsr_period_seconds2years()
        self.p.cy.updateFields() #Update the text fields in cycloids tab.
        self.p.nb.GetCurrentPage().update_parameters() #Update the current page's fields

    def saveFile(self,filename):
        f = open(filename,'w')
        for p,v in self.p.sc.parameters.items():
            if v or v == 'to_plot_many_cycloids': #Don't want to save to_plot_many_cycloids simply because this option shouldn't be loaded since the cycloids from the cycloids file aren't saved
                f.write(p + ' = ' + str(v) + '\n')
        f.close()

    def onQuit(self, evt):
        self.Close()


    def onRights(self, evt):
        # indentation (lack thereof) necessary to prevent tab spaces every newline in source code
        # not sure if the need for such indentation or lack thereof is b/c of python or wx
        # alternative is to use concatentation
        spiel = u"""ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any \
commercial use must be negotiated with the Office of Technology Transfer at the \
California Institute of Technology. \n\n
This software may be subject to U.S. export control laws and regulations. \
By accepting this document, the user agrees to comply with all applicable \
U.S. export laws and regulations. User has the responsibility to obtain export \
licenses, or other export authority as may be required before exporting such \
information to foreign countries or providing access to foreign persons. """

        copyright = "Copyright 2016, by the California Institute of Technology."
        #Update year whenever a new version is released.

        self.makeMsgDialog(spiel, copyright)

    def onDevelopment(self, evt):
        spiel = u"""SatStressGUI v4.0 was developed at the Jet Propulsion Laboratory, \
California Institute of Technology and is based on SatStressGUI. \
SatStressGUI was developed by the Planetary Geology Research group at the University of Idaho \
SatStressGUI is based on SatStress, which was designed by Zane Selvans and is available at \
http://code.google.com/p/satstress and most recently at https://github.com/zaneselvans/satstress \
\n\n SatStressGUI 4.0 has been created upon efforts by \
Alex Patthoff, Robert Pappalardo, Jonathan Kay, Lee Tang, \
Simon Kattenhorn, C.M. Cooper, Emily S. Martin, \
David Dubois, Ben J. Ayton, Jessica B. Li, \
Andre Ismailyan, Peter Sinclair."""

        self.makeMsgDialog(spiel, u'Developers')

    def onUpdates(self, evt):
        updates = u"""This is Version 4.0 of SatStressGUI.  For more information, please visit: \n\n\
https://github.com/SatStressGUI/SatStressGUI\n\n\
In this version, several bugs were fixed, and a new stressing mechanism (Polar Wander) was added.\
To find detailed notes of all the changes, suggest improvements, or report bugs, please visit the GitHub page."""

        self.makeMsgDialog(updates, u'Version 4.0')

    def onContacts(self, evt):
        # Create a message dialog box
        self.makeMsgDialog(u"Alex Patthoff via patthoff@jpl.nasa.gov",
                           u"Primary Contact")

    def onDiurnalref(self, evt):
        Resources = u"""Diurnal tidal stresses arise when a satellite is in an eccentric orbit. \
This is due to two reasons. \
First, the amplitude of the planet's gravitational force is greater at periapse than it is at apoapse. \
Secondly, the planet is rotating slightly faster (compared to its synchronous rotation rate) at periapse \
and slightly slower (again compared to its synchronous rotation rate) at apoapse. \
This results in a 'librational tide', where the planet appears to rock back and forth in the sky.\n\n\
For more information on diurnal tides, please see:\n\
Wahr, J., Z. A. Selvans, M. E. Mullen, A. C. Barr, G. C. Collins, \
M. M. Selvans, and R. T. Pappalardo, Modeling stresses on satellites due to non-synchronous rotation \
and orbital eccentricity using gravitational potential theory, \
Icarus, Volume 200, Issue 1, March 2009, Pages 188-206.
"""
        self.makeMsgDialog(Resources, u'About Diurnal Tides')

    def onNSRref(self, evt):
        Resources = u"""Nonsynchronous rotation (NSR) occurs when a satellite's lithosphere is decoupled from its core. \
When this happens, the tidal bulge of the shell causes it to experience a net torque, and could rotate more quickly than the synchronous rate. \
Thus, the planet appears to move across the sky, and the tidal bulge moves beneath the shell. \
This results in surface stresses. \
The period of this rotation should be > 10,000 years.\n\n\
For more information on NSR, please see:\n\
Wahr, J., Z. A. Selvans, M. E. Mullen, A. C. Barr, G. C. Collins, \
M. M. Selvans, and R. T. Pappalardo, Modeling stresses on satellites due to non-synchronous rotation \
and orbital eccentricity using gravitational potential theory, \
Icarus, Volume 200, Issue 1, March 2009, Pages 188-206.
"""
        self.makeMsgDialog(Resources, u'About Nonsynchronous Rotation')

    def onObliquityref(self, evt):
        Resources = u"""A satellite's obliquity (or axial tilt) is the angle between it rotational axis and its orbital axis. \
A satellite of zero obliquity will have a rotational axis perpendicular to its orbital plane. \
However, when the obliquity is nonzero, it causes the stresses due to diurnal tides and non-synchronous rotation to be asymmetric.\n\n\
For more information on stresses due to oblique orbits, see:\n\
Jara-Orue, H. M., & Vermeersen, B. L. (2011). Effects of low-viscous layers and a non-zero \
obliquity on surface stresses induced by diurnal tides and non-synchronous rotation: The \
case of Europa. Icarus, 215(1), 417-438.
"""
        self.makeMsgDialog(Resources, u'About Olibque Orbits')

    def onISTref(self, evt):
        Resources = u"""As satellites age, they could become cooler. \
This would result in more of the liquid ocean freezing, increasing the thickness of the icy crust. \
This process would force the ice shell to expand, putting extensional stress on the surface.\n\n\
For more information on Ice Shell Thickening as a stressing mechanism, please see:\n\
Nimmo, F. (2004). Stresses generated in cooling viscoelastic ice shells: Application \
to Europa. Journal of Geophysical Research: Planets (1991-2012), 109(E12).
"""
        self.makeMsgDialog(Resources, u'About Ice Shell Thickening')


    def onPWref(self, evt):
        Resources = u"""
Polar Wander is the apparent movement of a satellite's rotational pole due to nonsynchronous reorientation of the satellite's crust. \
If a satellite's crust is not coupled to its core, it may experience nonsynchronous rotation (NSR). \
Sometimes, this also results in a reorientation of the poles. \
The north pole appears to wander over the surface as the crust reorients itself. \
This results in stressing, due to the tidal bulge of the core and ocean moving beneath the crust, \
as well as the parent planet appearing to change its location in the sky. \n\n\
This stressing mechanism is calculated using an elastic model.\n\n\
For more information on Polar Wander as a stressing mechanism, please see:\n\
    Matsuyama, Isamu, and Francis Nimmo. "Tectonic patterns on reoriented and despun planetary bodies." Icarus 195, no. 1 (2008): 459-473.\n\
    Matsuyama, Isamu, Francis Nimmo, and Jerry X. Mitrovica. "Planetary reorientation." Annual Review of Earth and Planetary Sciences 42 (2014): 605-634.
"""
        self.makeMsgDialog(Resources, u'About Polar Wander')

    def onCycloidsref(self, evt):
        Resources = u""" Cycloids are arcuate lineaments found on the surface of Europa.  \
They are thought to be created when a fracture in the ice is propagated because of the stresses. \
In order for a cycloid to be created, the tensile stress at the location must exceed the tensile strength of the ice.\
Once the fracture has started, it will propagate through the ice at a certain velocity.\
This velocity could be constant, or could vary depending on the magnitude of the stress.\
During the cycloid's propagation, the satellite will continue orbiting around its primary.\
This causes the stress field on the satellite to change, making the cycloids curve.\
When the stress is no longer greater than the requisite propagation strength, the cycloid stops moving.\
If the stress reaches the propagation strength again, it will continue.\n\n\
For more information, please see:\n\
    Hoppa, G.V., Tufts, B.R., Greenberg, R., Geissler, P.E., 1999b. Formation of cycloidal \
features on Europa. Science 285, 1899-1902"""
        self.makeMsgDialog(Resources, u'About Cycloids')

    def onTutorial(self, evt):
        Tutorial = u"""Welcome to SatStressGUI!  This program is designed to model stresses icy satellites \
experience as they orbit their primary.  For more information on this program and the mathematics behind it, \
check the "Information" menu. \n\n\
1) Input the satellite's physical parameters on the Satellite tab.\n\
2) Select which stresses to apply in the Stresses tab.\n\
- When using Diurnal and NSR, either input Love numbers and check the box marked "Input Love Numbers", or \
leave them blank to allow the program to calculate Love numbers based on the satellite's physical properties.\n\
- Obliquity must be used with either Diurnal or NSR.\n\
3) In the Grid tab, input a latitude and longitude range to examine.\n\
- The number of grid points must be equal for both latitude and longitude.\n\
4) Also in the Grid tab, input the relevant information for the selected stresses.\n\
5) Change to the Plot tab to see the stress maps.\n\
- For more information on how to use the maps, see "Plot" in the Help Menu.\n\
6) Use the Point tab to calculate the stress at discrete points in space and time.
"""
        self.makeMsgDialog(Tutorial, u'Getting Started')

    def onHelpSat(self, evt):
        Help = u"""The Satellite Tab is used to input the physical properties of the satellite.\n\n\
- Each entry should use the units denoted in the square brackets next to the box.\n\
- The viscoelastic model used assumes that the satellite has two icy layers, a liquid ocean, and a solid core.\n\
- The NSR period is usually on the order of 100,000 years.  If you are not using NSR, you can leave it as 'infinity'.\n\
- The orbital eccentricity must be < 0.25.  Otherwise the program cannot reasonably calculate stresses.\n\
- If you have changed a number, but nothing seems to happen, try hitting 'Enter' in the box you changed.\n\
"""
        self.makeMsgDialog(Help, u'The Satellite Tab')

    def onHelpStresses(self, evt):
        Help = u"""The Stresses Tab is used to select which stresses to use.\n\n\
- For Diurnal and NSR stresses, the h2, k2, and l2 boxes should be left blank, unless the user wants to input their own values. \
Checking the "Input Love Numbers" box will allow you to use custom Love numbers. \
When inputting custom love numbers, you must use the format <Re> + <Im>j.  Do not use scientific notation. \
1.2 + 3e-05j would look like 1.2+0.00003j.\n\
- The Obliquity stress must be used with Diurnal or NSR.\n\
- The Thermal Diffusivity of the Ice Shell Thickening stress does not currently function.\n\
- Polar Wander uses an elastic, time-independent calculation, so it should probably not be used with other stresses.\n\
- By turning on the "Assume tidally locked satellite" option, the program will calculate the tidal axis as always perpendicular to the rotational axis.\n\
- If you turn off the tidal locking option and the plot does not update, press 'Enter' in each of the tidal axis text boxes.\n\
- Activating the "Despinning" box allows the user to change the initial and final rotation rate of the satellite.  \
The rotational period should be input in units of hours.\n\
- All coordinates should be input as latitude and longitude; conversion to colatitude is handled by the program.
"""
        self.makeMsgDialog(Help, u'The Stresses Tab')

    def onHelpPoint(self, evt):
        Help = u"""The Point Tab can be used to calculate the stress at any discrete point in space and time.\n\n\
- Enter a latitude, longitude, year, and orbital position for each point.\n\
- Press the "Calculate Stress" button.\n\
- Use the "Save to File" button to save the results as a .cvs file.\n\n\
- θ: Latitude (-90.00 to 90.00) [°]\n\
- φ: Longitude (-180.00 to180.00 (positive West or East to choose from)) [°]\n\
- t: Time since periapse (Periapse = 0) [yrs], used for secular stress calculations\n\
- orbital pos: Orbital position since periapse (Periapse = 0) [°], used for diurnal stress calculations\n\
- Stt: East-West component of stress field [kPa]\n\
- Spt: Off diagonal component of stress field [kPa]\n\
- Spp: North-South component of stress field [kPa]\n\
- σ1: Maximum tension [kPa]\n\
- σ3: Maximum compression [kPa]\n\
- α: The angle between σ1 and due north (clockwise is positive) [°]
"""
        self.makeMsgDialog(Help, u'The Point Tab')


    def onHelpGrid(self, evt):
        Help = u"""The Grid Tab is used to specify what section of the satellite to look at.\n\n\
- For more information about each stress, see the Information menu.\n\
- NOTE: The number of latitude and longitude grid points must be equal.\n\
- To examine the whole moon, use a latitude range from -90 to 90 and a longitude range of -180 to 180.\n\
- Each row will only activate when the appropriate stress is enabled.\n\
- The "Orbital Position" row is used to track diurnal stress from the satellite's orbit.  The satellite starts at the minimum position, and moves to the maximum position. \
Inputting 0 to 360 degrees will be one full orbit.  Additional orbits can be added by increasing the maximum beyond 360 degrees.\n\
- The map will occasionally not work for certain positions.  If this happens, simply change the number of increments or the end position.\n\
- The "Amount of NSR Buildup" row is used to determine how long the ice shell has been rotating. \
The Start Time is when the plotting starts, and the End Time is when the plotting ends.\n\
"""
        self.makeMsgDialog(Help, u'The Grid Tab')

    def onHelpCycloids(self, evt):
        Help = u"""The Cycloids Tab allows the user to generate a cycloidal feature on the map.\n\n\
- The cycloids are modeled and plotted on the Plot Tab.\n\
- The Yield Threshold is how much stress must be put on the crust to break the ice and initiate a fracture.\n\
- The Propagation Strength is how much stress must be put on the crust to make the split continue, and the split continues at the Propagation Speed.\n\
- The Starting Latitude and Longitude determine where the cycloid begins, and the Direction determines the curvature of the cycloid.\n\
- NOTE: The Vary Velocity option is currently untested.\n\
- For more information on cycloids, see the Information menu.
"""
        self.makeMsgDialog(Help, u'The Cycloids Tab')

    def onHelpPlot(self, evt):
        Help = u"""The Plot Tab shows a map of the stresses on the surface of the satellite.\n\n\
- Tension on the map is shown as positive, and compression as negative.
- You can step through the plots by using the buttons to the bottom right of the graph.\n\
- Each individual plot can be saved by using the save button to the lower left of the graph, and the series can be saved using the "Save Series" \
button to the lower right.\n\
- The panel on the right allows manipulation of the map, changing the scale and type of map, as well as the stresses showed.\n\
- The bottom panel enables and disables cycloids.\n\
- When using Polar Wander, the initial and final locations of the rotational poles and/or sub- and anti-jove points will appear on the graph.\n\
  - The initial North and South poles will be white circles.\n\
  - The final North and South poles will be black circles.\n\
  - The initial sub- and anti-jove points will be white squares.\n\
  - The final sub- and anti-jove points will be black squares.\n\
- The vectors created by Polar Wander do not currently appear to be generating correctly.\n\
- When using cycloids, if the program is unable to initiate a cycloid, it will plot a black triangle at the attempted location.\n\
  - If it creates a split, but cannot propagate it, it will plot a white triangle at the location.\n\
- Cycloids can be saved as Shape files via the appropriate button.  Loading of shape files is currently not supported.\n\
- NOTE: The cycloids cannot be saved as netcdf files currently.\n\
- NOTE: The Lineaments features does not function currently.
"""
        self.makeMsgDialog(Help, u'The Plot Tab')

    def makeMsgDialog(self, msg, title):
        msg = wx.MessageDialog(self, msg, title, wx.OK | wx.ICON_INFORMATION)
        msg.ShowModal()
        msg.Destroy()

    # Makes sure the user was intending to quit the application
    # at some point, make this conditional to if not changes have been made, no popup
    def OnCloseFrame(self, event):
        if self.p.sc.saveable_changed():
            dialog = wx.MessageDialog(self,
                message = "To save your parameters and/or plot, return to the relevant tab and click the appropriate button",
                caption = "Are you sure you want to quit without saving?")
            response = dialog.ShowModal() # show and disallows other input until closed

            if (response == wx.ID_OK):
                self.exit(event)
            else:
                event.StopPropagation()
        # if all saveable parameters have been saved, no need for popup window
        else:
            self.exit(event)

    def exit(self, evt):
        """ seems like a 
        """
        sys.exit(0)

# ===============================================================================
# 
class SatStressApp(wx.App):
    def OnInit(self):
        self.frame = SatStressFrame(None, title=u'SatStressGUI V4.0', size=(800,800))
        self.frame.Show(True)
        self.SetTopWindow(self.frame)
        return True

def main():
    #make Mac OS app be able to run calcLoveWahr4Layer from Resources
    #directory in application bundle
    os.environ['PATH'] += os.path.pathsep+os.path.abspath(os.curdir)
    app = SatStressApp(1) # The 0 aka false parameter means "don't redirect stdout and stderr to a window"    app.MainLoop()
    app.MainLoop()


if __name__ == '__main__':
    main()

#GNU Terry Pratchett