from SatPanel import *
import numpy
from Dialog import error_dialog
from exc import LocalError
from Config import *
import traceback

import matplotlib
matplotlib.interactive(False)
# sets matplotlib backend to 'WXAgg'
matplotlib.use('WXAgg')

# for plotting
#    basemap allows cartographic
from mpl_toolkits import basemap

"""
Polar Wander slider is currently disabled.  All of the code has been left in, but it is commented out.
If the calculations for Polar Wander are improved to be viscoelastic, they can be re-implemented.
-Peter Sinclair, 2016
"""

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
