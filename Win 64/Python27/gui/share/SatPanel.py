import wx
# ===============================================================================
# SUPERCLASS
# ===============================================================================
from RadioBox2 import RadioBox2
class SatPanel(wx.Panel):
    """
    Class that serves as the superclass of all panels of GUI.
    Having all tabs share a SatelliteCalculation object under superclass SatPanel
    allows information to be passed from one panel to another.
    """
    def __init__(self, *args, **kw):
        self.sc = kw['satellite_calculation']
        del kw['satellite_calculation']
        # why are you removing the attribute?
        super(SatPanel, self).__init__(*args, **kw)
        # super(Class, self).method(arg)
        self.parameters = {}
    def bind_parameters(self):
        for k, p in self.parameters.items():
            f = self.mk_change_param(k)
            if isinstance(p, wx.TextCtrl):
                p.Bind(wx.EVT_KILL_FOCUS, f)
                p.Bind(wx.EVT_TEXT, f)
            elif isinstance(p, wx.CheckBox):
                p.Bind(wx.EVT_CHECKBOX, f)
            elif isinstance(p, wx.ComboBox):
                p.Bind(wx.EVT_COMBOBOX, f)
            elif isinstance(p, RadioBox2):
                p.Bind(wx.EVT_RADIOBOX, f)
            elif isinstance(p, list):
                pass

    def bind_list_parameter(self, param_name):
        '''
        @param - a parameter such as theta, phi which has as a value a list of text ctrls
        Parameters which have list of ctrls as values rather than just one ctrl, should call this instead. The 0th index is
        just the label for all the ctrls
        '''
        for i in range(1, len(self.parameters[param_name])):
            f = self.mk_change_list_param(param_name, i)
            self.parameters[param_name][i].Bind(wx.EVT_TEXT, f)

    def mk_change_list_param(self,param_name, i ):
        def on_change(evt):
            self.sc.set_parameter(param_name, self.parameters[param_name][i].GetValue(), i)

        return on_change

    def mk_change_param(self, k, i=0):
        def on_change(evt):
            if (not isinstance(self.parameters[k], list)):
                if (not k in self.sc.parameters or (k in self.sc.parameters and not (self.sc.parameters[k] == self.parameters[k].GetValue() ) ) ):
                    self.sc.set_parameter(k, self.parameters[k].GetValue())

            else:
                self.sc.set_parameter(k, self.parameters[k][i].GetValue(), i)
        return on_change

    def update_parameters(self):
        for p, ctrl in self.parameters.items():
            try:
                if type(ctrl) is list:
                    for i in range(1,len(ctrl)):
                        ctrl[i].SetValue(self.sc.parameters[p][i - 1])
                else:
                    if isinstance(ctrl, wx.CheckBox):
                        if self.sc.parameters[p] == 'True' or self.sc.parameters[p] == '1' or self.sc.parameters[p]:
                            ctrl.SetValue(True)
                        else:
                            ctrl.SetValue(False)
                    else:
                        ctrl.SetValue(self.sc.parameters[p])
            except KeyError:
                print str(KeyError)
    def get_panel_param(self, p):
        return self.parameters[p]
