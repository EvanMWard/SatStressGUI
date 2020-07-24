# coding=utf-8
# ===============================================================================
# POINT TAB
# ===============================================================================
from wx.lib import scrolledpanel

from SatPanel import SatPanel
from Adders import *
from WrapStaticText import WrapStaticText
import traceback
from Dialog import *
from exc import LocalError
import csv
import os

seconds_in_year = 31556926.0

class PointPanel(SatPanel):
    """
        Defines the point panel of the GUI
        """
    def params_grid(self, panel, params_d, defval, width=3, row = 1):
        pp = wx.FlexGridSizer(row, width)
        add_table_header(panel, pp, params_d)
        self.parameters.update(add_text_ctrls(panel, pp, params_d, rows = row, point = True))
        for p,d in params_d:
            self.bind_list_parameter(p)
            self.sc.parameters[p] = []
            for i in range(row):
                self.sc.parameters[p].append(defval)
        return pp

    def __init__(self, *args, **kw):
        super(PointPanel, self).__init__(*args, **kw)
        #change self.rows to change how many rows are displayed in the GUI
        self.rows = 20
        self.sc.set_parameter('point_rows',self.rows)

        #parameter name and their labels
        self.header1 = [('theta', [u'θ [°]']), ('phi', [u'φ [°]']), ('t', [u't [yrs]']), ('orbit', [u'orbital pos [°]'])]
        self.header2 = [("Ttt", [u'Stt [kPa]']), ("Tpt", [u'Spt [kPa]']), ("Tpp", [u'Spp [kPa]'])]
        self.header3 = [("s1", [u'σ1 [kPa]']), ("s3", [u'σ3 [kPa]']), ("a", [u'α [°]'])]

        sz = wx.BoxSizer(orient=wx.VERTICAL)

        sz.Add(WrapStaticText(self, label=
        u'This tab is for calculating the stress tensor at a location at the surface ' +\
        u'at a point in the orbit. It uses the Stresses tab to determine which ' +\
        u'stresses are being calculated.'), flag=wx.ALL|wx.EXPAND)

        sz.AddSpacer(20)
        self.fieldPanel = scrolledpanel.ScrolledPanel(self,-1, size=(1000,400), style=wx.SIMPLE_BORDER)
        self.fieldPanel.SetupScrolling()

        rsz = wx.BoxSizer(orient=wx.HORIZONTAL)

        p2 = wx.BoxSizer(orient=wx.VERTICAL)
        cp = wx.BoxSizer(orient=wx.HORIZONTAL)
        p0 = wx.BoxSizer(orient=wx.VERTICAL)
        p0.Add(wx.StaticText(self.fieldPanel, label=u'Time/space location'), flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.pp = self.params_grid(self.fieldPanel, self.header1, '0', width=4, row=self.rows)
        p0.Add(self.pp)
        cp.Add(p0)
        p1 = wx.BoxSizer(orient=wx.VERTICAL)
        p1.Add(wx.StaticText(self.fieldPanel, label=u'Stress Tensor at a point'), flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.tp = self.params_grid(self.fieldPanel,self.header2, '', row = self.rows)
        p1.Add(self.tp, 1, wx.ALL|wx.EXPAND)
        cp.AddSpacer(15)
        cp.Add(p1)
        p3 = wx.BoxSizer(orient=wx.VERTICAL)
        p3.Add(wx.StaticText(self.fieldPanel, label=u'principal Components'), flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.sp = self.params_grid(self.fieldPanel,self.header3, '', row = self.rows)
        p3.Add(self.sp, 1, wx.ALL|wx.EXPAND)
        cp.Add(p3)
        p2.Add(cp)

        rsz.Add(p2)
        self.fieldPanel.SetSizer(rsz)
        sz.Add(self.fieldPanel)

        bp = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.spin_value = self.rows
        self.row_ctrl = wx.SpinCtrl(self, min = 1, value = str(self.rows), style=wx.TE_PROCESS_ENTER)

        self.save_b = wx.Button(self, label=u'Save to File')
        self.b = wx.Button(self, label=u'Calculate Stress')
        self.load_b = wx.Button(self, label=u'Load from file')
        bp.Add(self.b, 1, wx.ALL|wx.EXPAND, 3)
        bp.Add(self.load_b, 1, wx.ALL|wx.EXPAND, 3)
        bp.Add(self.save_b, 1, wx.ALL|wx.EXPAND, 3)

        bp.Add(WrapStaticText(self, label=u'Rows: '), flag = wx.ALIGN_CENTER_VERTICAL)
        bp.Add(self.row_ctrl)
        sz.Add(bp)

        sz.AddSpacer(15)

        self.SetSizer(sz)

        self.row_ctrl.Bind(wx.EVT_SPINCTRL, self.spinCtrl)
        self.row_ctrl.Bind(wx.EVT_TEXT, self.spinCtrl)
        #self.row_ctrl.Bind(wx.EVT_SPIN_DOWN, lambda evt, szr = pp: self.spin_down(evt, szr))
        # Here we bind the load and save buttons to the respective events
        wx.EVT_BUTTON(self, self.b.GetId(), self.on_calc)
        self.load_b.Bind(wx.EVT_BUTTON, self.load)
        self.save_b.Bind(wx.EVT_BUTTON, self.save)
        self.update_parameters()
        self.bind_parameters()
        self.updating = False
        for i in range(1, self.rows + 1):
            self.parameters['orbit'][i].Bind(wx.EVT_KILL_FOCUS, lambda evt, row = i: self.on_orbit_update(evt, row))
            self.parameters['orbit'][i].Bind(wx.EVT_TEXT, lambda evt, row = i: self.on_orbit_update(evt, row))
            self.parameters['t'][i].Bind(wx.EVT_KILL_FOCUS, lambda evt, row = i: self.on_t_update(evt, row))
            self.parameters['t'][i].Bind(wx.EVT_TEXT, lambda evt, row = i: self.on_t_update(evt, row))


    #updates the orbit text ctrls when t is changed
    def on_t_update(self, evt, row = 1):
        self.updating = True
        try:
            self.sc.set_parameter('t', self.parameters['t'][row].GetValue(), point = row)
            sat = self.sc.get_satellite()
            o = str(float(self.sc.parameters['t'][row - 1])/sat.orbit_period()*360.0*seconds_in_year)
            self.parameters['orbit'][row].SetValue(o)
            self.sc.set_parameter('orbit', o, point = row)
        except:
            traceback.print_exc()
        self.updating = False

    #updates the t text ctrls when orbital position is changed
    def on_orbit_update(self, evt, row = 1):
        self.updating = True
        try:
            self.sc.set_parameter('orbit', self.parameters['orbit'][row].GetValue(), point = row)
            sat = self.sc.get_satellite()
            t = str(float(self.sc.parameters['orbit'][row - 1])/360.0*sat.orbit_period()/seconds_in_year)
            self.parameters['t'][row].SetValue(t)
            self.sc.set_parameter('t', t, point = row)
        except:
            traceback.print_exc()
        self.updating = False

    def on_calc(self, evt):
        try:
            self.b.SetFocus()
            self.sc.calc_tensor(self.rows)
            self.update_parameters()
        except LocalError, e:
            error_dialog(self, str(e), e.title)

    #These functions were meant to handle events generated by the spin control used to change number of points to calculate
    def spinCtrl(self, evt):
        spin_value = evt.GetEventObject().GetValue()
        if spin_value == '':
            spin_value = 1

        if (int(spin_value) > self.rows):
            self.onUp(int(spin_value))
        else:
            self.spin_down(int(spin_value))

    def onUp(self, spin_value):
        self.pp.SetRows(spin_value)
        self.tp.SetRows(spin_value)
        self.sp.SetRows(spin_value)
        for i in range(spin_value-self.rows):
            self.add_row(self.fieldPanel, self.pp, self.header1, '0')
            self.add_row(self.fieldPanel,self.tp, self.header2, '')
            self.add_row(self.fieldPanel,self.sp, self.header3, '')
        self.rows = spin_value

        self.fieldPanel.Layout()
        self.Layout()
        self.update_parameters()
        self.fieldPanel.SetupScrolling()
        self.sc.set_parameter('point_rows',self.rows)
        for p,d in self.header1+self.header2+self.header3:
            self.bind_list_parameter(p)

    def add_row(self, panel, sz, params_d, defaultval):
        for p,d in params_d:
            text = wx.TextCtrl(self.fieldPanel, style = wx.TE_PROCESS_ENTER)
            sz.Add(text, flag=wx.ALL|wx.EXPAND)
            self.parameters[p].append(text)
            self.sc.parameters[p].append(defaultval)

    def spin_down(self, spin_value):
        self.pp.SetRows(spin_value)
        self.tp.SetRows(spin_value)
        self.sp.SetRows(spin_value)
        for i in range(self.rows - spin_value):
            for p,d in self.header1+self.header2+self.header3:
                self.parameters[p][-1].Destroy()
                del self.parameters[p][-1]
                del self.sc.parameters[p][-1]
        self.rows = spin_value
        self.sc.set_parameter('point_rows',self.rows)
        self.fieldPanel.Layout()

    def load(self, evt):
        try:
            file_dialog(self,
                        message=u"Load from CSV file",
                        style=wx.OPEN,
                        wildcard='CSV files (*.csv)|*.csv',
                        action=self.load_entries)
        except Exception, e:
            traceback.print_exc()

    def set_num_rows(self,num_rows):
        self.pp.SetRows(num_rows)
        self.sp.SetRows(num_rows)
        self.tp.SetRows(num_rows)
        if (num_rows > self.rows):
            for j in range(num_rows-self.rows):
                self.add_row(self.fieldPanel,self.pp, self.header1, '0')
                self.add_row(self.fieldPanel,self.tp, self.header2, '')
                self.add_row(self.fieldPanel,self.sp, self.header3, '')
            self.update_parameters()
        else:
            for j in range(self.rows-num_rows):
                for p,d in self.header1+self.header2+self.header3:
                    self.parameters[p][-1].Destroy()
                    del self.parameters[p][-1]
                    del self.sc.parameters[p][-1]
        self.rows = num_rows
        self.row_ctrl.SetValue(num_rows)
        self.spin_value = num_rows
        self.sc.set_parameter('point_rows',self.rows)
        self.fieldPanel.Layout()
        self.fieldPanel.SetupScrolling()

    def load_entries(self, filename):
        f = open(filename)
        csvreader = csv.reader(f)
        coord = csvreader.next()  #Skip headers
        data = list(csvreader)
        self.set_num_rows(len(data))
        try:
            keys = ['theta', 'phi', 't', 'orbit']

            for i,coord in enumerate(data):

                for key in keys:
                    val = coord[keys.index(key)]
                    self.parameters[key][i+1].SetValue(val)
                    self.sc.set_parameter(key, val, point = i+1)
        except:
            traceback.print_exc()
        finally:
            f.close()
            self.fieldPanel.Layout()
            self.fieldPanel.SetupScrolling()
            self.Layout()

    #opens save dialog
    def save(self, evt):
        file_dialog(self,
                    message=u"Save to CSV file",
                    style=wx.SAVE,
                    wildcard='CSV files (*.csv)|*.csv',
                    defaultFile='untitled.csv',
                    action=self.save_pointcalc)

    #parses text ctrls and writes to csv
    def save_pointcalc(self, filename=None):
        tmp = False
        if not filename:
            filename = os.tempnam(None, 'csv')
            tmp = True
        f = open(filename, 'wb')
        writer = csv.writer(f)
        headers = [u'theta [degrees]', 'phi [degrees]', 't [yrs]', 'orbital pos [degrees]', \
                   'Stt [kPa]', 'Spt [kPa]', 'Spp [kPa]', \
                   'sigma1 [kPa]', 'sigma3 [kPa]', 'alpha [degrees]']
        writer.writerow(headers)
        keys = ['theta', 'phi', 't', 'orbit',\
        "Ttt", "Tpt", "Tpp", \
        "s1", "s3", "a"]
        for i in range(self.rows):
            row = [self.sc.parameters[key][i] for key in keys]
            writer.writerow(row)
        f.close()
