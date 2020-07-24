import wx
# ===============================================================================
# Class, used on most panels
# ===============================================================================
class WrapStaticText(wx.StaticText):
    """
    Constrains number of characters in a line in the GUI to less than 58
    """
    def __init__(self, *args, **kw):
        super(WrapStaticText, self).__init__(*args, **kw)
        self._label = self.GetLabel()
        self._rewrap()
        wx.EVT_SIZE(self, self.change_size)
    def _rewrap(self):
        w = self.GetSize().width
        self.SetLabel(self._label)
        if w > 50:
            self.Wrap(w)
    def change_size(self, evt):
        self.Unbind(wx.EVT_SIZE)
        self._rewrap()
        self.Bind(wx.EVT_SIZE, self.change_size)