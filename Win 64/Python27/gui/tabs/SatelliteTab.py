import wx

# ===============================================================================
# SATELLITE TAB
# ===============================================================================
from ComboBox2 import *
from SatPanel import SatPanel
from Dialog import *
from Adders import *


class SatelliteLayersPanel(SatPanel):
    """
    Defines the satellite layers panel of the GUI
    """
    def __init__(self, *args, **kw):
        super(SatelliteLayersPanel, self).__init__(*args, **kw)

        sz = wx.BoxSizer(orient=wx.VERTICAL)

        filler = wx.BoxSizer(wx.HORIZONTAL)
        filler.AddSpacer(15)

        top = wx.BoxSizer(orient=wx.VERTICAL)

        bp = wx.BoxSizer(orient=wx.HORIZONTAL)
        load_b = wx.Button(self, label=u'Load from file')
        save_b = wx.Button(self, label=u'Save to file')
        bp.Add(load_b, 1, wx.ALL|wx.EXPAND, 3)
        bp.Add(save_b, 1, wx.ALL|wx.EXPAND, 3)

        # satellite parameters
        # FlexGridSizer organizes visual elements into grid layout
        sp = wx.FlexGridSizer(1,2)
        self.parameters = add_parameters_to_sizer(self, sp, self.sc.satellite_vars)
        # layers parameters
        lp = wx.FlexGridSizer(1, len(self.sc.layer_vars_d))
        add_static_texts(self, lp, self.sc.layer_vars_d)
        lv = []
        for l, v in self.sc.satlayers_d:
            for p, d in self.sc.layer_vars_d:
                #if p == 'YOUNG':
                #    lv.append(('LAME_MU_%d' % l, ''))
                #if p == 'POISSON':
                #    lv.append(('LAME_LAMBDA_%d' % l, ''))
                lv.append(("%s_%d" % (p, l), ''))
        self.parameters.update(add_text_ctrls(self, lp, lv))
        self.update_parameters()
        self.bind_parameters()
        for l, v in self.sc.satlayers_d:
            self.parameters["LAYER_ID_%d" % l].SetEditable(True)
        #for l, v in self.sc.satlayers_d:
            #self.parameters["TENSILE_STR_%d" % l].SetValue('0')
            #self.parameters["TENSILE_STR_%d" % l].Disable()
        # end

        top.Add(bp, 0, wx.ALL|wx.EXPAND)
        top.Add(filler)
        top.Add(sp)

        sz.Add(top)
        sz.Add(filler)
        sz.Add(lp)

        sz.AddSpacer(10)
        HelpText = wx.StaticText(self, label=u'For help in using this program, select "Getting Started", in the Help menu.')
        HelpFont = wx.Font(18, wx.DEFAULT, wx.NORMAL, wx.BOLD) # Sets the font and size of the text. -PS 2016
        HelpText.SetFont(HelpFont)
        sz.Add(HelpText)
        sz.AddSpacer(10)
        # This text was added to provide new users with important information. -PS 2016
        sz.Add(wx.StaticText(self, label=u'This model makes several assumptions when calculating stresses.'))
        sz.Add(wx.StaticText(self, label=u'-The body is assumed to be composed of four layers, with the third layer being a liquid ocean.'))
        sz.Add(wx.StaticText(self, label=u'-It is assumed to behave in a viscoelastic manner.'))
        sz.Add(wx.StaticText(self, label=u'-Each layer is considered to be homogenous throughout, with no differences in density or thickness based on location, but decreasing in mass out from the core.'))
        sz.Add(wx.StaticText(self, label=u'-The Polar Wander stress assumes that the body is in a circular, zero-inclination, synchronous orbit.'))
        sz.Add(wx.StaticText(self, label=u'-Polar Wander stress is calculated using an elastic model.'))
        sz.Add(wx.StaticText(self, label=u'-The orbit is assumed to have an eccentricity of <0.25, and the primary\'s mass be at least 10 times the satellite\'s mass.'))

        self.SetSizer(sz)
        wx.EVT_BUTTON(self, load_b.GetId(), self.load)
        wx.EVT_BUTTON(self, save_b.GetId(), self.save)

    def load(self, evt):
        try:
            file_dialog(self,
                message=u"Load from satellite file",
                style=wx.OPEN,
                wildcard='Satellite files (*.satellite;*.sat)|*.satellite;*.sat',
                action=self.load_entries)
        except Exception, e:
            error_dialog(self, str(e), u'Satellite Error')

    def save(self, evt):
        file_dialog(self,
            message=u"Save to satellite file",
            style=wx.SAVE | wx.OVERWRITE_PROMPT,
            wildcard='Satellite files (*.satellite;*.sat)|*.satellite;*.sat',
            defaultFile='satellite.satellite',
            action=self.sc.save_satellite)

    def load_entries(self, filename):
        self.sc.load_satellite(filename)
        self.update_parameters()
