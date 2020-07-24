import wx
# ===============================================================================
# Class & related Functions for ComboBox (dropdown choice bar) implementation
# ===============================================================================
class ComboBox2(wx.ComboBox):
    """
    Custom implementation of wx.ComboBox
    """
    def __init__(self, parent, id=-1, value='',
        pos=wx.DefaultPosition, size=wx.DefaultSize,
        choices=[], style=0, validator=wx.DefaultValidator,
        name=wx.ChoiceNameStr):
        self.__choice_map = dict([ (c,d) for c,d in choices ])
        self.__reverse_map = dict([ (d,c) for c,d in choices ])
        super(ComboBox2, self).__init__(parent, id=id, value=self.__choice_map[value],
            pos=pos, size=size, choices=[ d for c,d in choices ],
            style=style, validator=validator, name=name)

    def GetValue(self):
        return self.__reverse_map[super(ComboBox2, self).GetValue()]

    def SetValue(self, val):
        super(ComboBox2, self).SetValue(self.__choice_map[val])

# adds labels to GUI from the parameters dictionary
def add_parameters_to_sizer(parent, sz, parameters_d):
    parameters = {}
    for p, d in parameters_d:
        sz.Add(wx.StaticText(parent, label=d), flag=wx.ALIGN_CENTER_VERTICAL)
        parameters[p] = wx.TextCtrl(parent, style=wx.TE_PROCESS_ENTER)
        sz.Add(parameters[p], flag=wx.EXPAND|wx.ALL)
    return parameters

# creates an instance of ComboBox2 in the GUI
def add_combobox2_to_sizer(parent, sz, parameter, description, choices):
    parameters = {}
    sz.Add(wx.StaticText(parent, label=description), flag=wx.ALIGN_CENTER_VERTICAL)
    parameters[parameter] = ComboBox2(parent,
        value=parent.sc.parameters.get(parameter),
        choices=choices,
        style=wx.CB_DROPDOWN | wx.CB_READONLY)
    sz.Add(parameters[parameter])
    return parameters

# adds wx.Checkbox objects to GUI labeled with parameters
def add_checkboxes_to_sizer(parent, sz, parameters_d):
    parameters = {}
    for p, d in parameters_d:
        parameters[p] = wx.CheckBox(parent, label=d)
        sz.Add(parameters[p], flag=wx.ALIGN_CENTER_VERTICAL)
    return parameters