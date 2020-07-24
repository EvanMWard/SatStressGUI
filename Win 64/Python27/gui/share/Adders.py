from RadioBox2 import *
#adds an instances of RadioBox2 to the GUI, labeled with given parameters
def add_radiobox2_to_sizer(parent, sz, parameter, description, choices):
    parameters = {}
    sz.Add(wx.StaticText(parent, label=description))
    parameters[parameter] = RadioBox2(parent, choices=choices)
    sz.Add(parameters[parameter], flag=wx.EXPAND|wx.ALL)
    return parameters

#adds given paramters as static text to the GUI
def add_static_texts(parent, sz, parameters_d):
    sts = [ wx.StaticText(parent, label=d, style=wx.TE_PROCESS_ENTER) for p, d in parameters_d ]
    for st in sts:
        sz.Add(st, flag=wx.ALIGN_CENTER)
    return sts

def add_table_header(parent, sz, parameters_d):
    sts = [ wx.StaticText(parent, label=d[0], style=wx.TE_PROCESS_ENTER) for p, d in parameters_d ]
    for st in sts:
        sz.Add(st, flag=wx.ALIGN_CENTER)
    return sts

#takes a dictionary of parameters and adds it to a wx.TextCtrl object
def add_text_ctrls(parent, sz, parameters_d, rows = 1, point=False):
    parameters = {}
    if point:
        for i in range(rows):
            for p, d in parameters_d:
                if i == 0:
                    parameters[p] = d
                text = wx.TextCtrl(parent, style = wx.TE_PROCESS_ENTER)
                parameters[p].append(text)
                sz.Add(text, flag=wx.ALL|wx.EXPAND)
        return parameters
    else:
        for p, d in parameters_d:
            parameters[p] = wx.TextCtrl(parent, style = wx.TE_PROCESS_ENTER)
            sz.Add(parameters[p], flag=wx.ALL|wx.EXPAND)
        return parameters
