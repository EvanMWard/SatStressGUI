# ===============================================================================
# CYCLOID TAB
# ===============================================================================
from ..share.SatPanel import SatPanel
import random
from ..share.WrapStaticText import WrapStaticText
from ..share.Dialog import *
import os
from ..exc import LocalError
import csv

class CycloidsPanel(SatPanel):
    """
    Defines the cycloids panel of the GUI
    NTS: Should restrict decimal places at some pt
    """
    def __init__(self, *args, **kw):
        super(CycloidsPanel, self).__init__(*args, **kw)
        self.cyc = None
        self.textCtrls = {} #keys are the name of the fields (YIELD,etc) and values are the textCtrl objects for those fields.
        # initialize sizers
        sz = wx.BoxSizer(wx.VERTICAL)
        gridSizer    = wx.FlexGridSizer(rows=7, cols=2, hgap=5, vgap=0)

        dirSizer = wx.BoxSizer(wx.HORIZONTAL)
        buttonSizer = wx.BoxSizer(wx.HORIZONTAL)
        filler = wx.BoxSizer(wx.HORIZONTAL)
        varyvSizer = wx.BoxSizer(wx.HORIZONTAL)
        # whould eventually replace with add_combobox2_to_sizer()
        # create combobox that chooses (initial?) direction
        which_dir = wx.StaticText(self, wx.ID_ANY, 'Propagation Direction: ')
        dirSizer.Add(which_dir, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL, 0)
        all_dir = ['East', 'West']
        self.start_dir = wx.ComboBox(self, size=(100, 50) ,choices=all_dir, style=wx.CB_DROPDOWN|wx.CB_READONLY)
        # bind
        self.Bind(wx.EVT_COMBOBOX, self.EvtSetDir, self.start_dir)
        self.parameters['STARTING_DIRECTION'] = self.start_dir

        # create load/save buttons
        save_bt = wx.Button(self, label='Save to file')
        save_bt.Bind(wx.EVT_BUTTON, self.on_save_cyclparams)
        load_bt = wx.Button(self, label='Load from file')
        load_bt.Bind(wx.EVT_BUTTON, self.on_load_cyclparams)
        buttonSizer.Add(load_bt, wx.ALIGN_CENTER, 10)
        buttonSizer.AddSpacer(5)
        buttonSizer.Add(save_bt, wx.ALIGN_CENTER)

        self.vary = wx.CheckBox(self, wx.ID_ANY, 'Vary Velocity   k = ')
        self.Bind(wx.EVT_CHECKBOX, self.EvtSetVary, self.vary)
        self.parameters['VARY_VELOCITY'] = self.vary

        self.constant = wx.TextCtrl(self, wx.ID_ANY, '0', style=wx.TE_PROCESS_ENTER)
        self.Bind(wx.EVT_TEXT, self.EvtSetConstant, self.constant)
        self.constant.Disable()
        self.parameters['k'] = self.constant

        self.use_multiple = wx.CheckBox(self, wx.ID_ANY, 'Use loaded CSV file')
        self.Bind(wx.EVT_CHECKBOX, self.EvtSetUseMultiple, self.use_multiple)
        self.use_multiple.Disable()
        self.parameters['to_plot_many_cycloids'] = self.use_multiple

        varyvSizer.Add(self.vary, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL)
        varyvSizer.Add(self.constant, 0, wx.ALL|wx.ALIGN_CENTER_VERTICAL)

        # add widgets into grid
        # Set the TextCtrl to expand on resize

        fieldsToAdd = [ ('cycloid_name', 'Cycloid Name'), ('YIELD', 'Yield (Threshold) [kPa]: '),('PROPAGATION_STRENGTH','Propagation Strength [kPa]: '),('PROPAGATION_SPEED','Propagation Speed [m/s]: '), ('STARTING_LONGITUDE', 'Starting Longitude: '), ('STARTING_LATITUDE', 'Starting Latitude: ')]

        self.textCtrls.update(self.add_text_control(gridSizer, fieldsToAdd ))
        self.sc.parameters['cycloid_name'] = ""

        gridSizer.Add(dirSizer)
        gridSizer.Add(self.start_dir)
        gridSizer.Add(varyvSizer)
        many_params = wx.Button(self, label='Load Multiple Cycloid Parameters')
        wx.EVT_BUTTON(self, many_params.GetId(), self.load_many)

        # add to overarching sizer sz

        sz.Add(WrapStaticText(self,
            label=u'This tab calculates cycloids through combined diurnal and NSR stresses. Cycloids ' +
            u'are arcuate lineaments found on the surface of Europa. ' +
            u'They can be modeled and plotted on the following ' +
            u'Plot tab. The Yield Strength is the threshold that initiates fracture in the ice. ' +
            u'This fracture will propagate as long as the strength is below this threshold and greater than the ' +
            u'Propagation Strength. The Propagation Speed is usually <10 m/s. ' +
            u'For further information on cycloids see the Help menu.'),
            flag=wx.ALL|wx.EXPAND)
            #sz.Add(filler)

        sz.Add(buttonSizer, 0, wx.ALL, 5)
        #sz.Add(filler2)
        sz.Add(gridSizer, 0, wx.ALL|wx.EXPAND, 5)
        sz.Add(self.use_multiple,0, wx.ALL|wx.EXPAND, 5)
        sz.Add(many_params)
        self.SetSizer(sz)
        sz.Fit(self)

    def add_text_control(self,  sz, parameters_d):
        txtCtrls = {}
        for p, d in parameters_d:
            sz.Add(wx.StaticText(self, label=d), flag=wx.ALIGN_CENTER_VERTICAL)
            txtCtrlObj = wx.TextCtrl(self, -1,name=p)
            #txtCtrlObj.Bind(wx.EVT_CHAR,self.OnChar)
            txtCtrlObj.Bind(wx.EVT_TEXT, self.OnText)
            txtCtrls[p] = txtCtrlObj
            sz.Add(txtCtrlObj, flag=wx.EXPAND|wx.ALL)
            self.parameters[p] = txtCtrlObj
        return txtCtrls

    def OnChar(self,event):
        # Checks to make sure only numbers are entered into a box.
        # Some characters are allowed, for imaginary numbers and scientific notation. -PS 2016
        charEntered= event.GetKeyCode()
        if (charEntered >= 48 and charEntered <= 57) or charEntered == 8 or charEntered == 9 or charEntered == 13 or charEntered == 45 or charEntered ==46:
            event.Skip()

    def OnText(self,event):
        self.sc.cycloid_changed = True
        if not event.GetEventObject().GetValue() == 'None':
            try:
                self.sc.parameters[event.GetEventObject().GetName()] = float(event.GetEventObject().GetValue())
            except:
                self.sc.parameters[event.GetEventObject().GetName()] = event.GetEventObject().GetValue()
        else:
            self.sc.parameters[event.GetEventObject().GetName()] = None

    def load_many(self, evt):
        try:
            file_dialog(self,
                message=u'Load from .csv file',
                style=wx.OPEN,
                wildcard=u'*.csv',
                action=self.load_many_params)
        except LocalError, e:
            error_dialog(self, str(e), e.title)

    def on_save_cyclparams(self, evt):
        try:
            file_dialog(self,
                message = u'Save cycloid parameters to file',
                style = wx.SAVE | wx.OVERWRITE_PROMPT,
                wildcard = 'Cycloid files (*.cyc)|*.cyc',
                defaultFile = 'cycloid_params.cyc',
                action = self.save_cyclparams)
        except Exception, e:
            error_dialog(self, str(e), u'Error saving cycloid parameters')

    def save_cyclparams(self, filename):
        tmp = False
        if filename is None:
            filename = os.tempnam(None, 'grid')
            tmp = True
        f = open(filename, 'w')
        for k,v in self.sc.cycloid_parameters_d.items():
            if k == 'VARY_VELOCITY' and not v:
                f.write(k + " = False" + "\n")
            else:
                if self.sc.parameters.has_key(k):
                    f.write(k + " = " + str(self.sc.parameters[k]) + "\n")
                else:
                    f.write(k + " = None" + "\n")
        f.close()
        if not tmp:
            self.grid_save_changed = False
        return filename, tmp

    def on_load_cyclparams(self, evt):
        try:
            file_dialog(self,
                message = u"Load cycloid parameters from file",
                style = wx.OPEN,
                wildcard = 'Cycloid files (*.cyc)|*.cyc',
                action = self.load_cyclparams)
        except Exception, e:
            error_dialog(self, str(e), u'Error loading cycloid parameters')

    def load_cyclparams(self, filename):
        try:
            f = open(filename)
        except:
            error_dialog(self, 'File error', 'Cannot open file')
        for p, v in nvf2dict(f).items():
            if not p in ('k','VARY_VELOCITY', 'STARTING_DIRECTION'):
                if v == 'None':
                    self.sc.parameters[p] = 'None'
                else:
                    self.sc.parameters[p] = float(v)
            elif p == 'k':
                if v == 'None':
                    self.sc.parameters[p] = 0
                    self.constant.SetValue(0)
                else:
                    self.sc.parameters[p] = float(v)
            elif p == 'VARY_VELOCITY':
                if v == 'True' or v == '1':
                    self.constant.Enable()
                    self.vary.SetValue(1)
                else:
                    self.constant.Disable()
                self.sc.parameters[p] = v
            elif p == 'STARTING_DIRECTION':
                self.sc.parameters[p] = v
                self.start_dir.SetValue(v)
        self.updateFields()
        self.cycloid_saved = True
        f.close()

    #For loading multiple cycloids
    def load_many_params(self, filename):
        #the CSV headers need to be the same name as the parameters for the Cycloids init object (threshold, propgataion_speed, etc)
        self.use_multiple.Enable()
        self.use_multiple.SetValue(True)
        self.EvtSetUseMultiple(None)
        self.sc.parameters['to_plot_many_cycloids'] = True
        self.sc.many_changed = True
        paramFile = open(filename, 'rU')
        rows = list(csv.reader(paramFile))
        params_to_load = rows[0]
        self.sc.params_for_cycloids = {}
        i = 0
        for row in rows[1:]:
            self.sc.params_for_cycloids[i] = {}
            for j, param in enumerate(params_to_load):
                self.sc.params_for_cycloids[i].update({param: str(row[j]) })
                self.sc.params_for_cycloids[i].update({'degree_step':0.3}) # degree step can be higher
            i += 1
        paramFile.close()

    def updateFields(self):
        if self.sc.parameters['VARY_VELOCITY'] == 'True' or self.sc.parameters['VARY_VELOCITY'] == '1':
            self.vary.SetValue(True)
            if self.sc.parameters.has_key('k'):
                self.constant.Enable()
                self.constant.SetValue(str(self.sc.parameters['k']))
        if self.sc.parameters.has_key('STARTING_DIRECTION'):
            self.start_dir.SetValue(self.sc.parameters['STARTING_DIRECTION'])
        for p, textctrl in self.textCtrls.items():
            if self.sc.parameters.has_key(p):
                    textctrl.SetValue(str(self.sc.parameters[p]))

    def load_shape(self, filename):
        # Loading shapefiles is currently not supported.  -PS 2016
        # walk around char const * restriction
        sf = os.path.splitext(str(filename))[0] + '.shp'
        self.loaded['data'] = shp2lins(sf, stresscalc=self.calc)
        self.loaded['lines'] = []
        d = wx.ColourDialog(self, self.loaded['color'])
        if (d.ShowModal() == wx.ID_OK):
            self.loaded['color'] = d.GetColourData()
        self.plot()

    def EvtSetDir(self, event):
        self.sc.parameters['STARTING_DIRECTION'] = event.GetString()

    def EvtSetYeild(self, event):
        assert(float(event.GetString() > 0))
        self.sc.parameters['YIELD'] = float(event.GetString())

    def EvtSetPropStr(self, event):
        assert(float(event.GetString() > 0))
        self.sc.parameters['PROPAGATION_STRENGTH'] = float(event.GetString())

    def EvtSetPropSpd(self, event):
        assert(float(event.GetString() > 0))
        self.sc.parameters['PROPAGATION_SPEED'] = float(event.GetString())

    def EvtSetVary(self, event):
        self.sc.parameters['VARY_VELOCITY'] = self.vary.GetValue()
        if self.vary.GetValue():
            self.constant.Enable()
        else:
            self.constant.Disable()

    def EvtSetConstant(self, event):
        self.sc.parameters['k'] = float(event.GetString())

    def EvtSetUseMultiple(self, event):
        if self.use_multiple.GetValue():
            self.sc.parameters['to_plot_many_cycloids'] = True
            for ctrl in [self.parameters[p] for p in self.sc.cycloid_parameters_d]:
                ctrl.Disable()
        else:
            self.sc.parameters['to_plot_many_cycloids'] = False
            self.use_multiple.SetValue(False)
            for ctrl in [self.parameters[p] for p in self.sc.cycloid_parameters_d]:
                ctrl.Enable()

    def EvtSetStartLat(self, event):
        lat = float(event.GetString())
        assert(lat <= 90)
        assert(lat >= -90)
        self.sc.parameters['STARTING_LATITUDE'] = float(event.GetString())

    def EvtSetStartLon(self, event):
        lon = float(event.GetString())
        assert(lon <= 180)
        assert(lon >= -180)
        self.sc.parameters['STARTING_LONGITUDE'] = float(event.GetString())

    def EvtRandLat(self, event):
        # generates random lat to the 2nd decimal place (current precision of GUI)
        rand_startlat = float("%.2f" % random.uniform(-90, 90))
        # set it to parameter
        self.sc.parameters['STARTING_LATITUDE'] = rand_startlat
        # display it in textctrl
        input_startlat.SetValue('%s', rand_startlat)

    def EvtRandLon(self, event):
        rand_startlon = float("%.2f" % random.uniform(-180, 180))
        # set it to parameters
        self.sc.parameters['STARTING_LONGITUDE'] = rand_startlon
        # display in textctrl
        input_startlon.SetValue('%s', rand_startlon)
