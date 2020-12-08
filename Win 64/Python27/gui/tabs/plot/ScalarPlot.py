# coding=utf-8
import os
import time

import scipy.ndimage

import satstress
from Adders import add_radiobox2_to_sizer
from ComboBox2 import add_combobox2_to_sizer, add_checkboxes_to_sizer, add_parameters_to_sizer
from Dialog import file_dialog
from WrapStaticText import WrapStaticText
from PlotTab import *
from StressPlot import *
from cycloid import Cycloid
from lineament import *
from stressplot import scalar_grid

seconds_in_year = 31556926.0 # 365.24 days

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
            style = wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT,
            #style = wx.SAVE | wx.OVERWRITE_PROMPT, -- This may be deprecated, check back after shp2lins has been written - EW 2020
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
            action = "SaveCycloidAsShape")

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