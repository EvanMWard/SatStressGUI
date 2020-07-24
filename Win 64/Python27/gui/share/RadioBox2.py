import wx
# ===============================================================================
# Class & Related Functions for RadioBox
# ===============================================================================
class RadioBox2(wx.Control):
    """
    Custom implementation of wx.RadioBox
    """
    def __init__(self, parent, id=-1,
        pos=wx.DefaultPosition, size=wx.DefaultSize,
        choices=[], value=None, style=0, validator=wx.DefaultValidator,
        name=wx.RadioBoxNameStr):
        self.__choice_map = dict([ (c,d) for c,d in choices ])
        super(RadioBox2, self).__init__(parent, id=id,
            pos=pos, size=size, style=style,
            validator=validator, name=name)
        sz = wx.BoxSizer(orient=wx.VERTICAL)
        self.__radiobuttons = {}
        c, d = choices[0]
        self.__radiobuttons[c] = wx.RadioButton(self, label=d, style=wx.RB_GROUP)
        for c, d in choices[1:]:
            self.__radiobuttons[c] = wx.RadioButton(self, label=d)
            if c == value:
                self.__radiobuttons[c].SetValue(True)
        for b in self.__radiobuttons.values():
            b.Bind(wx.EVT_RADIOBUTTON, self._send_radiobox_event)
        for c, d in choices:
            sz.Add(self.__radiobuttons[c])
        self.SetMinSize((20, 40))
        self.SetSizer(sz)
        self.Layout()
        self.Fit()

    def GetValue(self):
        for c, b in self.__radiobuttons.items():
            if b.GetValue():
                return c
        return None

    def _send_radiobox_event(self, evt):
        e = wx.PyEvent(self.GetId(), wx.wxEVT_COMMAND_RADIOBOX_SELECTED)
        self.value = self.GetValue()
        self.GetEventHandler().ProcessEvent(e)

    def SetValue(self, val):
        self.__radiobuttons[val].SetValue(True)
