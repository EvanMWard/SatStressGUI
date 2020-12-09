# coding=utf-8
from ..share.SatPanel import SatPanel
from ..share.WrapStaticText import *
from ..share.ComboBox2 import *
import re
from ..share.Dialog import *
from ..exc import LocalError
from ..share.Adders import *
# ===============================================================================
# STRESSES TAB
# ===============================================================================
class StressListPanel(SatPanel):
    """
    Defines the stresses panel of the GUI. Contains stress type selections and allows
    user to input their own love numbers for further calculations
    """
    def __init__(self, *args, **kw):
        super(StressListPanel, self).__init__(*args, **kw)

        filler = wx.BoxSizer(wx.HORIZONTAL)
        filler.AddSpacer(15)

        topsizer = wx.BoxSizer(wx.VERTICAL)
        othersz = wx.BoxSizer(wx.HORIZONTAL)

        sz = wx.BoxSizer(orient=wx.VERTICAL)

        topsizer.Add(WrapStaticText(self, label=
            u'Select the stress types to use in further computation, such as Love numbers, stress tensors, plotting of stress trajectories.'),
            0, wx.ALL|wx.EXPAND)
        topsizer.AddSpacer(10)

        sz.AddSpacer(23)

        # for Diurnal
        self.parameters = add_checkboxes_to_sizer(self, sz, [ ('Diurnal', 'Diurnal') ])
        sz.AddSpacer(8)

        # for NSR
        self.parameters.update(add_checkboxes_to_sizer(self, sz,
            [ ('Nonsynchronous Rotation', 'Nonsynchronous Rotation') ]))

        sz.AddSpacer(8)
        # Added this text to provide users with important information. -PS 2016
        sz.Add(wx.StaticText(self, label=u'To input custom Love numbers, use the format <Re> +/- <Im>j.'))
        sz.Add(wx.StaticText(self, label=u'Do not use scientific notation when inputting custom Love numbers.'))
        sz.Add(wx.StaticText(self, label=u'"3.0-1.0e-03j" should be written as "3.0-0.001j".'))
        sz.Add(wx.StaticText(self, label=u"If no Love Numbers are input, the program will calculate them automatically."))

        sz.AddSpacer(8)

        # for Diurnal w/ Obliquity
        self.parameters.update(add_checkboxes_to_sizer(self, sz,
            [ ('Obliquity', 'Obliquity') ]))
        DiObliq_sz = wx.BoxSizer(wx.VERTICAL)
        # include arg of periapsis parameter for Diurnal w/ Obliquity
        peri_sz = wx.BoxSizer(orient=wx.HORIZONTAL)
        peri_sz.AddSpacer(28)
        self.periapsis_label = wx.StaticText(self,
           label=u'Argument of Periapsis [°]  ')
        peri_sz.Add(self.periapsis_label, flag=wx.ALIGN_CENTER_VERTICAL)
        self.parameters.update(add_text_ctrls(self, peri_sz,
           [ ('periapsis_arg', 'periapsis_arg') ]))
        DiObliq_sz.Add(peri_sz)
        DiObliq_sz.AddSpacer(5)
        # include degree of obliquity for Diurnal w/ Obliquity
        obliq_sz = wx.BoxSizer(orient=wx.HORIZONTAL)
        obliq_sz.AddSpacer(30)
        self.obliq_label = wx.StaticText(self, label=u'Degree of Obliquity [°]  ')
        obliq_sz.Add(self.obliq_label, flag=wx.ALIGN_CENTER_HORIZONTAL)
        obliq_sz.Add(filler)
        self.parameters.update(add_text_ctrls(self, obliq_sz,
            [ ('obliquity', 'obliquity') ]))
        obliq_sz.AddSpacer(5)
        DiObliq_sz.Add(obliq_sz)

        sz.Add(DiObliq_sz)

        self.parameters.update(add_checkboxes_to_sizer(self, sz,
            [ ('Ice Shell Thickening', 'Ice Shell Thickening') ]))
        ISTParams_sz = wx.BoxSizer(wx.VERTICAL)
        # include ice thickness parameter for IST aka Ice Shell Thickening
        delta_tc_sz = wx.BoxSizer(orient=wx.HORIZONTAL)
        delta_tc_sz.AddSpacer(28)
        self.delta_label = wx.StaticText(self, label=u'Change in Thickness [km] ')
        delta_tc_sz.Add(self.delta_label, flag=wx.ALIGN_CENTER_VERTICAL)
        self.parameters.update(add_text_ctrls(self, delta_tc_sz, [ ('delta_tc', 'delta_tc') ]))
        ISTParams_sz.Add(delta_tc_sz)
        ISTParams_sz.AddSpacer(5)
        # include thermal diffusivity parameter for IST
        #diff_sz = wx.BoxSizer(orient=wx.HORIZONTAL)
        #diff_sz.AddSpacer(28)
        #self.diffusivity_label = wx.StaticText(self, label=u'Thermal Diffusivity [m\u00b2/s]  ')
        #diff_sz.Add(self.diffusivity_label, flag=wx.ALIGN_CENTER_VERTICAL)
        #self.parameters.update(add_text_ctrls(self, diff_sz, [ ('diffusivity', 'diffusivity') ]))
        #ISTParams_sz.Add(diff_sz)
        sz.Add(ISTParams_sz)

        # Removed the diffusivity option from the GUI because it is not implemented. -PS 2016

        self.parameters.update(add_checkboxes_to_sizer(self, sz, [ ('Polar Wander', 'Polar Wander') ]))

        sz.Add(wx.StaticText(self, label=u"Polar wander is not completely tested, so it may not be accurate"))
        sz.Add(wx.StaticText(self, label=u"to combine it with other stresses."))
        sz.Add(wx.StaticText(self, label=u"The stress map from Polar Wander appears to be correct,"))
        sz.Add(wx.StaticText(self, label=u"but the principal stress vectors are rotated 180° for some reason."))

        Polargrid = wx.FlexGridSizer(rows=9, cols=3, hgap=3, vgap=5) # A GridSizer to hold the polar wander coordinates.  -PS 2016
        self.Latitude_label = wx.StaticText(self, label=u'Latitude [°]')
        self.Longitude_label = wx.StaticText(self, label=u'Longitude [°]')
        self.Blank_label = wx.StaticText(self, label=u' ')
        self.Blank_label2 = wx.StaticText(self, label=u' ')
        self.Blank_label3 = wx.StaticText(self, label=u' ')
        self.Blank_label4 = wx.StaticText(self, label=u' ')
        self.PoleInitial = wx.StaticText(self, label=u'Initial Pole Location')
        self.PoleFinal = wx.StaticText(self, label=u'Final Pole Location')
        self.TidalInitial = wx.StaticText(self, label=u'Initial Tidal Bulge Location')
        self.TidalFinal = wx.StaticText(self, label=u'Final Tidal Bulge Location')
        self.InitialSpin_label = wx.StaticText(self, label=u'Initial Period')
        self.FinalSpin_label = wx.StaticText(self, label=u'Final Period')

        self.PWthetaRi = wx.TextCtrl(self, wx.ID_ANY, '', style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT, self.set_thetaRi, self.PWthetaRi)
        self.PWphiRi = wx.TextCtrl(self, wx.ID_ANY, '', style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT, self.set_phiRi, self.PWphiRi)
        self.PWthetaRf = wx.TextCtrl(self, wx.ID_ANY, '', style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT, self.set_thetaRf, self.PWthetaRf)
        self.PWphiRf = wx.TextCtrl(self, wx.ID_ANY, '', style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT, self.set_phiRf, self.PWphiRf)
        self.PWthetaTi = wx.TextCtrl(self, wx.ID_ANY, '', style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT, self.set_thetaTi, self.PWthetaTi)
        self.PWphiTi = wx.TextCtrl(self, wx.ID_ANY, '', style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT, self.set_phiTi, self.PWphiTi)
        self.PWthetaTf = wx.TextCtrl(self, wx.ID_ANY, '', style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT, self.set_thetaTf, self.PWthetaTf)
        self.PWphiTf = wx.TextCtrl(self, wx.ID_ANY, '', style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT, self.set_phiTf, self.PWphiTf)
        self.TidalLock = wx.CheckBox(self, wx.ID_ANY, style=wx.ALIGN_RIGHT, label=u'Assume tidally locked satellite')
        self.TidalLock.SetValue(True)
        self.Bind(wx.EVT_CHECKBOX, self.Lock_Body, self.TidalLock)
        # Despinning could be implemented as a separate stress, but I did not have time to do this. -PS 2016
        self.DespinningBox = wx.CheckBox(self, wx.ID_ANY, style=wx.ALIGN_RIGHT, label=u'Despinning [hours]')
        self.Bind(wx.EVT_CHECKBOX, self.Despinning, self.DespinningBox)
        self.InitialSpinPeriod = wx.TextCtrl(self, wx.ID_ANY, '', style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT, self.set_InitialSpinPeriod, self.InitialSpinPeriod)
        self.FinalSpinPeriod = wx.TextCtrl(self, wx.ID_ANY, '', style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT, self.set_FinalSpinPeriod, self.FinalSpinPeriod)

        Polargrid.AddMany([
            (self.Blank_label, 0, wx.ALL|wx.EXPAND), (self.Latitude_label, 0, wx.ALL|wx.EXPAND), (self.Longitude_label, 0, wx.ALL|wx.EXPAND),
            (self.PoleInitial, 0, wx.ALL|wx.EXPAND), (self.PWthetaRi, 0, wx.ALL|wx.EXPAND), (self.PWphiRi, 0, wx.ALL|wx.EXPAND),
            (self.PoleFinal, 0, wx.ALL|wx.EXPAND), (self.PWthetaRf, 0, wx.ALL|wx.EXPAND), (self.PWphiRf, 0, wx.ALL|wx.EXPAND),
            (self.TidalInitial, 0, wx.ALL|wx.EXPAND), (self.PWthetaTi, 0, wx.ALL|wx.EXPAND), (self.PWphiTi, 0, wx.ALL|wx.EXPAND),
            (self.TidalFinal, 0, wx.ALL|wx.EXPAND), (self.PWthetaTf, 0, wx.ALL|wx.EXPAND), (self.PWphiTf, 0, wx.ALL|wx.EXPAND),
            (self.TidalLock, 0, wx.ALL|wx.EXPAND), (self.Blank_label2, 0, wx.ALL|wx.EXPAND), (self.Blank_label3, 0, wx.ALL|wx.EXPAND),
            (self.Blank_label4, 0, wx.ALL|wx.EXPAND), (self.InitialSpin_label, 0, wx.ALL|wx.EXPAND), (self.FinalSpin_label, 0, wx.ALL|wx.EXPAND),
            (self.DespinningBox, 0, wx.ALL|wx.EXPAND), (self.InitialSpinPeriod, 0, wx.ALL|wx.EXPAND), (self.FinalSpinPeriod, 0, wx.ALL|wx.EXPAND)])

        sz.Add(Polargrid)

        sz.AddSpacer(15)
        save_love_bt = wx.Button(self, label='Save Love numbers')
        wx.EVT_BUTTON(self, save_love_bt.GetId(), self.on_save_love)
        sz.Add(save_love_bt)

        ##### Create boxes for inputting love number #####
        ## todo: display calcuated love numbers in these boxes also ##
        grid = wx.FlexGridSizer(rows=4, cols=4, hgap=0, vgap=5)

        self.h2 = wx.StaticText(self, label=u'h\u2082')
        self.k2 = wx.StaticText(self, label=u'k\u2082')
        self.l2 = wx.StaticText(self, label=u'l\u2082')

        self.h2Diurn = wx.TextCtrl(self, wx.ID_ANY, '', style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT, self.set_h2Diurn, self.h2Diurn)
        self.k2Diurn = wx.TextCtrl(self, wx.ID_ANY, '', style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT, self.set_k2Diurn, self.k2Diurn)
        self.l2Diurn = wx.TextCtrl(self, wx.ID_ANY, '', style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT, self.set_l2Diurn, self.l2Diurn)
        self.userDiurn = wx.CheckBox(self, wx.ID_ANY, label='Input Love Numbers')
        self.Bind(wx.EVT_CHECKBOX, self.useUserLove_diurn, self.userDiurn)

        self.h2NSR = wx.TextCtrl(self, wx.ID_ANY, '', style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT, self.set_h2NSR, self.h2NSR)
        self.k2NSR = wx.TextCtrl(self, wx.ID_ANY, '', style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT, self.set_k2NSR, self.k2NSR)
        self.l2NSR = wx.TextCtrl(self, wx.ID_ANY, '', style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT, self.set_l2NSR, self.l2NSR)
        self.userNSR = wx.CheckBox(self, wx.ID_ANY, label='Input Love Numbers')
        self.Bind(wx.EVT_CHECKBOX, self.useUserLove_nsr, self.userNSR)

        grid.AddMany([
            (self.h2, 0, wx.ALL|wx.EXPAND), (self.k2, 0, wx.ALL|wx.EXPAND), (self.l2, 0, wx.ALL|wx.EXPAND), (self.Blank_label, 0, wx.ALL|wx.EXPAND),
            (self.h2Diurn, 0, wx.ALL|wx.EXPAND), (self.k2Diurn, 0, wx.ALL|wx.EXPAND), (self.l2Diurn, 0, wx.ALL|wx.EXPAND), (self.userDiurn, 0, wx.ALL|wx.EXPAND),
            (self.h2NSR, 0, wx.ALL|wx.EXPAND), (self.k2NSR, 0, wx.ALL|wx.EXPAND), (self.l2NSR, 0, wx.ALL|wx.EXPAND), (self.userNSR, 0, wx.ALL|wx.EXPAND)
            ])


        othersz.Add(sz, 5, wx.ALL|wx.EXPAND)
        othersz.Add(grid, 5, wx.ALL|wx.EXPAND)

        topsizer.Add(othersz, wx.ALL|wx.EXPAND)
        self.SetSizer(topsizer)

        self.update_parameters()
        self.bind_parameters()

        self.parameters['Diurnal'].Bind(wx.EVT_CHECKBOX, self.on_set_diurn)

        self.parameters['Nonsynchronous Rotation'].Bind(wx.EVT_CHECKBOX, self.on_set_nsr)

        self.parameters['Ice Shell Thickening'].Bind(wx.EVT_CHECKBOX, self.on_set_ist)

        self.parameters['Obliquity'].Bind(wx.EVT_CHECKBOX, self.on_set_obliq)

        self.parameters['Polar Wander'].Bind(wx.EVT_CHECKBOX,self.on_set_polar)

        # Causes the inputs for individual stresses to be disabled until the stresses are selected.  -PS 2016
        self.disable_display_diurnlove()
        self.disable_display_nsrlove()
        self.disable_istparams()
        self.disable_obliq()
        self.disable_polar()


    # These functions disable and enable various stress parameters on the GUI. -PS 2016
    def disable_display_diurnlove(self):
        for widg in [self.h2, self.k2, self.l2,
                     self.h2Diurn, self.k2Diurn, self.l2Diurn,
                     self.userDiurn]:
            widg.Disable()

    def enable_display_diurnlove(self):
        for widg in [self.h2, self.k2, self.l2,
                     self.h2Diurn, self.k2Diurn, self.l2Diurn,
                     self.userDiurn]:
            widg.Enable()

    def disable_display_nsrlove(self):
        for widg in [self.h2, self.k2, self.l2,
                     self.h2NSR, self.k2NSR, self.l2NSR,
                     self.userNSR]:
            widg.Disable()

    def enable_display_nsrlove(self):
        for widg in [self.h2, self.k2, self.l2,
                     self.h2NSR, self.k2NSR, self.l2NSR,
                     self.userNSR]:
            widg.Enable()

    def disable_istparams(self):
        for e in [self.delta_label, self.parameters['delta_tc']]:
            e.Disable()
        """
        for e in [self.delta_label, self.parameters['delta_tc'],
                  self.diffusivity_label, self.parameters['diffusivity'] ]:
            e.Disable()
        """

    def enable_istparams(self):
        """Don't yet enable diffusivity as it is only relevant for the viscoelastic case."""
        for e in [self.delta_label, self.parameters['delta_tc'] ]:
            e.Enable()

    def disable_obliq(self):
        for e in [self.obliq_label, self.parameters['obliquity'],
               self.periapsis_label, self.parameters['periapsis_arg'] ]:
            e.Disable()

    def enable_obliq(self):
        for e in [self.obliq_label, self.parameters['obliquity'],
               self.periapsis_label, self.parameters['periapsis_arg'] ]:
            e.Enable()

    def enable_polar(self):
        for e in [
         self.PWthetaRi, self.PWphiRi,
         self.PWthetaRf, self.PWphiRf,
         self.Longitude_label, self.Latitude_label,
         self.PoleInitial, self.PoleFinal,
         self.TidalLock,
         self.DespinningBox]:
            e.Enable()
        if not self.TidalLock.GetValue():
            for e in [
             self.TidalInitial, self.TidalFinal]:
                e.Enable()
        if self.DespinningBox.GetValue():
            for e in [
            self.InitialSpinPeriod, self.FinalSpinPeriod,
            self.InitialSpin_label, self.FinalSpin_label]:
                e.Enable()

    def disable_polar(self):
        for e in [
         self.PWthetaRi, self.PWphiRi,
         self.PWthetaRf, self.PWphiRf,
         self.PWthetaTi, self.PWphiTi,
         self.PWthetaTf, self.PWphiTf,
         self.Longitude_label, self.Latitude_label,
         self.PoleInitial, self.PoleFinal,
         self.TidalInitial, self.TidalFinal,
         self.TidalLock,
         self.DespinningBox,
         self.InitialSpinPeriod, self.FinalSpinPeriod,
         self.InitialSpin_label, self.FinalSpin_label]:
            e.Disable()

    def on_set_diurn(self, evt):
        state = self.parameters['Diurnal'].GetValue()
        self.sc.set_parameter('Diurnal', state)
        if state:
            self.enable_display_diurnlove()
        else:
            self.disable_display_diurnlove()

    def on_set_nsr(self, evt):
        state = self.parameters['Nonsynchronous Rotation'].GetValue()
        self.sc.set_parameter('Nonsynchronous Rotation', state)
        if state:
            self.enable_display_nsrlove()
        else:
            self.disable_display_nsrlove()

    def on_set_ist(self, evt):
        s = self.parameters['Ice Shell Thickening'].GetValue()
        self.sc.set_parameter('Ice Shell Thickening', s)
        if s:
            self.enable_istparams()
        else:
            self.disable_istparams()

    def on_set_obliq(self, evt):
        name = 'Obliquity'
        s = self.parameters[name].GetValue()
        self.sc.set_parameter(name, s)
        if s:
            self.enable_obliq()
        else:
            self.disable_obliq()

    def on_set_polar(self,evt):
        s = self.parameters['Polar Wander'].GetValue()
        self.sc.set_parameter('Polar Wander', s)
        if s:
            self.enable_polar()
        else:
            self.disable_polar()

    def parse_complex(self, string):
        real, imag = re.split(r'[+-]', string)
        if imag.startswith('i') or imag.startswith('j'):
            return float(real), float(imag[1:])
        elif imag.endswith('i') or imag.endswith('j'):
            return float(real), float(imag[:-1])

    def useUserLove_diurn(self, evt):
        if self.userDiurn.GetValue():
            self.sc.stress_d['Diurnal'].useUser = True
        else:
            self.sc.stress_d['Diurnal'].useUser = False

    def useUserLove_nsr(self, evt):
        if self.userDiurn:
            self.sc.stress_d['Nonsynchronous Rotation'].useUser = True
        else:
            self.sc.stress_d['Nonsynchronous Rotation'].useUser = False

    # These functions are used to pass values from the GUI to satstress.py without using the file-saving method.  -PS 2016

    def set_h2Diurn(self, evt):
        self.sc.stresses_changed = True
        self.sc.stress_d['Diurnal'].loveUser.update_h2(self.parse_complex(evt.GetString()))

    def set_k2Diurn(self, evt):
        self.sc.stresses_changed = True
        self.sc.stress_d['Diurnal'].loveUser.update_k2(self.parse_complex(evt.GetString()))

    def set_l2Diurn(self, evt):
        self.sc.stresses_changed = True
        self.sc.stress_d['Diurnal'].loveUser.update_l2(self.parse_complex(evt.GetString()))

    def set_h2NSR(self, evt):
        self.sc.stresses_changed = True
        self.sc.stress_d['Nonsynchronous Rotation'].loveUser.update_h2(self.parse_complex(evt.GetString()))

    def set_k2NSR(self, evt):
        self.sc.stresses_changed = True
        self.sc.stress_d['Nonsynchronous Rotation'].loveUser.update_k2(self.parse_complex(evt.GetString()))

    def set_l2NSR(self, evt):
        self.sc.stresses_changed = True
        self.sc.stress_d['Nonsynchronous Rotation'].loveUser.update_l2(self.parse_complex(evt.GetString()))

    def set_thetaRi(self, evt):
        self.sc.stresses_changed = True
        self.sc.stress_d['Polar Wander'].UserCoordinates.update_thetaRi(float(evt.GetString()))
        self.sc.polarwander_coordinates['thetaRInitial'] = float(evt.GetString())
        if self.sc.polarwander_coordinates['Locked']:
            if float(evt.GetString()) >= 0:
                self.sc.polarwander_coordinates['thetaTInitial'] = float(evt.GetString()) - 90
            else:
                self.sc.polarwander_coordinates['thetaTInitial'] = float(evt.GetString()) + 90

    def set_phiRi(self, evt):
        self.sc.stresses_changed = True
        self.sc.stress_d['Polar Wander'].UserCoordinates.update_phiRi(float(evt.GetString()))
        self.sc.polarwander_coordinates['phiRInitial'] = float(evt.GetString())
        if self.sc.polarwander_coordinates['Locked']:
            self.sc.polarwander_coordinates['phiTInitial'] = float(evt.GetString())

    def set_thetaRf(self, evt):
        self.sc.stresses_changed = True
        self.sc.stress_d['Polar Wander'].UserCoordinates.update_thetaRf(float(evt.GetString()))
        self.sc.polarwander_coordinates['thetaRFinal'] = float(evt.GetString())
        if self.sc.polarwander_coordinates['Locked']:
            if float(evt.GetString()) >= 0:
                self.sc.polarwander_coordinates['thetaTFinal'] = float(evt.GetString()) - 90
            else:
                self.sc.polarwander_coordinates['thetaTFinal'] = float(evt.GetString()) + 90

    def set_phiRf(self, evt):
        self.sc.stresses_changed = True
        self.sc.stress_d['Polar Wander'].UserCoordinates.update_phiRf(float(evt.GetString()))
        self.sc.polarwander_coordinates['phiRFinal'] = float(evt.GetString())
        if self.sc.polarwander_coordinates['Locked']:
            self.sc.polarwander_coordinates['phiTFinal'] = float(evt.GetString())

    def set_thetaTi(self, evt):
        self.sc.stresses_changed = True
        self.sc.stress_d['Polar Wander'].UserCoordinates.update_thetaTi(float(evt.GetString()))
        self.sc.polarwander_coordinates['thetaTInitial'] = float(evt.GetString())

    def set_phiTi(self, evt):
        self.sc.stresses_changed = True
        self.sc.stress_d['Polar Wander'].UserCoordinates.update_phiTi(float(evt.GetString()))
        self.sc.polarwander_coordinates['phiTInitial'] = float(evt.GetString())

    def set_thetaTf(self, evt):
        self.sc.stresses_changed = True
        self.sc.stress_d['Polar Wander'].UserCoordinates.update_thetaTf(float(evt.GetString()))
        self.sc.polarwander_coordinates['thetaTFinal'] = float(evt.GetString())

    def set_phiTf(self, evt):
        self.sc.stresses_changed = True
        self.sc.stress_d['Polar Wander'].UserCoordinates.update_phiTf(float(evt.GetString()))
        self.sc.polarwander_coordinates['phiTFinal'] = float(evt.GetString())

    def Lock_Body(self, evt):
        """
        If the "Assume tidally locked" option is selected, the program will calculate the tidal axis automatically.
        It will assume that the tidal axis should be 90 degrees apart from the rotational axis.
        -PS 2016
        """
        self.sc.stresses_changed = True
        if self.TidalLock.GetValue():
            for e in [
            self.PWthetaTi, self.PWphiTi,
            self.PWthetaTf, self.PWphiTf,
            self.TidalInitial, self.TidalFinal]:
                e.Disable()
            self.sc.stress_d['Polar Wander'].UserCoordinates.lock_body(True)

            if not self.sc.polarwander_coordinates['Locked']:
                self.sc.stress_d['Polar Wander'].UserCoordinates.update_thetaRi(self.sc.polarwander_coordinates['thetaRInitial'])
                if self.sc.polarwander_coordinates['thetaRInitial'] >= 0:
                    self.sc.polarwander_coordinates['thetaTInitial'] = self.sc.polarwander_coordinates['thetaRInitial'] - 90
                else:
                    self.sc.polarwander_coordinates['thetaTInitial'] = self.sc.polarwander_coordinates['thetaRInitial'] + 90
                self.sc.stress_d['Polar Wander'].UserCoordinates.update_phiRi(self.sc.polarwander_coordinates['phiRInitial'])
                self.sc.polarwander_coordinates['phiTInitial'] = self.sc.polarwander_coordinates['phiRInitial']

                self.sc.stress_d['Polar Wander'].UserCoordinates.update_thetaRf(self.sc.polarwander_coordinates['thetaRFinal'])
                if self.sc.polarwander_coordinates['thetaRFinal'] >= 0:
                    self.sc.polarwander_coordinates['thetaTFinal'] = self.sc.polarwander_coordinates['thetaRFinal'] - 90
                else:
                    self.sc.polarwander_coordinates['thetaTFinal'] = self.sc.polarwander_coordinates['thetaRFinal'] + 90
                self.sc.stress_d['Polar Wander'].UserCoordinates.update_phiRf(self.sc.polarwander_coordinates['phiRFinal'])
                self.sc.polarwander_coordinates['phiTFinal'] = self.sc.polarwander_coordinates['phiRFinal']

            self.sc.polarwander_coordinates['Locked'] = True

        else:
            for e in [
            self.PWthetaTi, self.PWphiTi,
            self.PWthetaTf, self.PWphiTf,
            self.TidalInitial, self.TidalFinal]:
                e.Enable()
            self.sc.stress_d['Polar Wander'].UserCoordinates.lock_body(False)
            self.sc.polarwander_coordinates['Locked'] = False

    def Despinning(self, evt):
        # Despinning is calculated by changing the rotation rate in the flattening coefficient for polar wander. -PS 2016
        self.sc.stresses_changed = True
        if self.DespinningBox.GetValue():
            self.sc.polarwander_coordinates['Despinning'] = True
            for e in [
             self.InitialSpinPeriod, self.FinalSpinPeriod,
             self.InitialSpin_label, self.FinalSpin_label]:
                e.Enable()
            self.sc.stress_d['Polar Wander'].UserCoordinates.spin_change(True)
        else:
            self.sc.polarwander_coordinates['Despinning'] = False
            for e in [
             self.InitialSpinPeriod, self.FinalSpinPeriod,
             self.InitialSpin_label, self.FinalSpin_label]:
                e.Disable()
            self.sc.stress_d['Polar Wander'].UserCoordinates.spin_change(False)

    def set_InitialSpinPeriod(self, evt):
        self.sc.stresses_changed = True
        self.sc.stress_d['Polar Wander'].UserCoordinates.update_InitialSpin(float(evt.GetString()))

    def set_FinalSpinPeriod(self, evt):
        self.sc.stresses_changed = True
        self.sc.stress_d['Polar Wander'].UserCoordinates.update_FinalSpin(float(evt.GetString()))

    def on_save_love(self, evt):
        # Saves the love numbers generated by the program.
        try:
            file_dialog(self,
                message=u"Save Love Numbers",
                style=wx.SAVE | wx.OVERWRITE_PROMPT,
                wildcard='Text files (*.txt)|*.txt',
                defaultFile='love.txt',
                action=self.sc.save_love)
        except LocalError, e:
            error_dialog(self, str(e), e.title)
