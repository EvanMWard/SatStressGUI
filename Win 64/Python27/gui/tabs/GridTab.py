# ===============================================================================
# GRID TAB
# ===============================================================================
from Adders import *
from ComboBox2 import *
from SatPanel import SatPanel
from Dialog import *
from exc import LocalError
class GridCalcPanel(SatPanel):
    """
    Defines the grid panel of the GUI
    """
    def __init__(self, *args, **kw):
        super(GridCalcPanel, self).__init__(*args, **kw)

        sz = wx.BoxSizer(orient=wx.VERTICAL)

        sz.Add(wx.StaticText(self, label = u'This tab is used to define the limits of the plot, both spatially and temporally.'))
        sz.AddSpacer(8)

        grid_id_p = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.parameters = add_parameters_to_sizer(self, grid_id_p, [('GRID_ID', u"Grid ID")])

        sb = wx.Button(self, label=u"Save grid")
        lb = wx.Button(self, label=u"Load grid")

        gmcp = wx.FlexGridSizer(0, len(self.sc.grid_vars_d) + 1)
        # grid points
        add_static_texts(self, gmcp, [('','')] + self.sc.grid_vars_d)
        for p, d in self.sc.grid_parameters_d[:2]:
            gmcp.Add(wx.StaticText(self, label=d))
            self.parameters.update(add_text_ctrls(self, gmcp, [ ("%s_%s" % (p,v), '') for v, dv in self.sc.grid_vars_d ]))

        for i in range(4):
            gmcp.AddSpacer(20)
        # orbital
        self.orbit_labels = add_static_texts(self, gmcp, [('',''), ('',u'Minimum'), ('',u'Maximum'), ('',u'Number of increments')])
        p, d = self.sc.grid_parameters_d[3]
        self.orbit_labels.append(wx.StaticText(self, label=d))
        gmcp.Add(self.orbit_labels[-1])
        self.parameters.update(
            add_text_ctrls(self, gmcp, [('%s_%s' % (p,v), '') for v,d1 in self.sc.grid_vars_d ]))
        # nsr time
        for i in range(4):
            gmcp.AddSpacer(20)
        self.nsr_labels = add_static_texts(self, gmcp,
            [('', ''), ('', u'Start Time [yrs]'), ('', u'End Time [yrs]'), ('', u'Number of increments')])
        self.nsr_labels.append(wx.StaticText(self, label=u'Amount of NSR build up'))
        gmcp.Add(self.nsr_labels[-1])
        self.parameters.update(
            add_text_ctrls(self, gmcp, [ ('TIME_MIN', ''), ('nsr_time', ''), ('TIME_NUM', '') ]))
        self.parameters['nsr_time'].SetMinSize((250, 10))
        top = wx.BoxSizer(orient=wx.HORIZONTAL)
        top.Add(grid_id_p)
        top.AddSpacer(6)
        top.Add(lb)
        top.AddSpacer(6)
        top.Add(sb)
        sz.Add(top)
        sz.AddSpacer(15)
        sz.Add(gmcp)
        sz.AddSpacer(15)
        sz.Add(wx.StaticText(self, label = u'Note: Number of latitude and longitude grid points must be equal'))
        sz.Add(wx.StaticText(self, label=u"Sometimes the map will not generate for certain diurnal orbit values."))
        sz.Add(wx.StaticText(self, label=u"If this happens, just change your number of increments or end value."))

        self.SetSizer(sz)

        self.update_parameters()
        self.bind_parameters()

        self.updating_range = False

        wx.EVT_BUTTON(self, sb.GetId(), self.save)
        wx.EVT_BUTTON(self, lb.GetId(), self.load)

        # Used to set default values on the grid tab.  For some reason it only seems to work for Diurnal and NSR.  -PS 2016
        self.orbital_set = 0
        self.nsr_set = 0
        self.parameters['GRID_ID'].SetValue('default')
        self.parameters['LAT_MIN'].SetValue('-90')
        self.parameters['LAT_MAX'].SetValue('90')
        self.parameters['LAT_NUM'].SetValue('10')
        self.parameters['LON_MIN'].SetValue('-180')
        self.parameters['LON_MAX'].SetValue('180')
        self.parameters['LON_NUM'].SetValue('10')


    def enable_nsr(self):
        for p in ['TIME_MIN', 'nsr_time', 'TIME_NUM']:
            self.parameters[p].Enable()
        for sts in self.nsr_labels:
            sts.Enable()
        if not self.nsr_set: #sets default NSR values -PS 2016
            self.parameters['TIME_MIN'].SetValue('0')
            self.parameters['nsr_time'].SetValue('1000000')
            self.parameters['TIME_NUM'].SetValue('10')
            self.nsr_set = 1

    def disable_nsr(self):
        for p in ['TIME_MIN', 'nsr_time', 'TIME_NUM']:
            self.parameters[p].Disable()
        for sts in self.nsr_labels:
            sts.Disable()

    def enable_orbit(self):
        for p in ['ORBIT_MIN', 'ORBIT_MAX', 'ORBIT_NUM']:
            self.parameters[p].Enable()
        for sts in self.orbit_labels:
            sts.Enable()
        if not self.orbital_set: #sets default Diurnal values -PS 2016
            self.parameters['ORBIT_MIN'].SetValue('0')
            self.parameters['ORBIT_MAX'].SetValue('360')
            self.parameters['ORBIT_NUM'].SetValue('10')
            self.orbital_set = 1

    def disable_orbit(self):
        for p in ['ORBIT_MIN', 'ORBIT_MAX', 'ORBIT_NUM']:
            self.parameters[p].Disable()
        for sts in self.orbit_labels:
            sts.Disable()

    def update_parameters(self):
        super(GridCalcPanel, self).update_parameters()
        if self.sc.parameters.get('Nonsynchronous Rotation', False):
            self.enable_nsr()
        else:
            self.disable_nsr()
        if self.sc.parameters.get('Diurnal', False) or \
                self.sc.parameters.get('Obliquity', False):
            self.enable_orbit()
        else:
            self.disable_orbit()
        for p in [ "%s_%s" % (p, v)
                  for p,pd in self.sc.grid_parameters_d
                  for v, vd in self.sc.grid_vars_d ]:
            if self.sc.parameters.get(p) is None:
                if p in self.parameters:
                    self.parameters[p].SetValue('')

    def load_grid(self, filename):
        self.sc.load_grid(filename)
        self.update_parameters()

    def save(self, evt):
        try:
            file_dialog(self,
                message=u"Save to grid file",
                style=wx.SAVE | wx.OVERWRITE_PROMPT,
                wildcard=u'Grid files (*.grid)|*.grid',
                defaultFile=self.sc.parameters['GRID_ID'] + '.grid',
                action=self.sc.save_grid)
        except KeyError, e:
            error_dialog(self, str(e) + ' not defined', 'Grid Error')
        except LocalError, e:
            error_dialog(self, str(e), e.title)

    def load(self, evt):
        try:
            file_dialog(self,
                message=u"Load from grid file",
                style=wx.OPEN,
                wildcard=u'Grid files (*.grid)|*.grid',
                action=self.load_grid)
        except LocalError, e:
            error_dialog(self, str(e), e.title)
