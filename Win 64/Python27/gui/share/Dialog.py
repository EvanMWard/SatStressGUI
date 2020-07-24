import wx
# ===============================================================================
# Global Functions showing dialog boxes given event
# ===============================================================================
def file_dir_dialog(parent, dialog_class, message="", style=wx.OPEN, action=None, **kw):
    fd = dialog_class(parent, message=message, style=style, **kw)
    wx.Yield()
    if (fd.ShowModal() == wx.ID_OK):
        action(fd.GetPath())
    fd.Destroy()

def file_dialog(parent, message="", style=wx.OPEN, action=None, **kw):
    file_dir_dialog(parent, wx.FileDialog, message, style, action, **kw)

def dir_dialog(parent, message="", style=wx.OPEN, action=None, **kw):
    file_dir_dialog(parent, wx.DirDialog, message, style, action, **kw)

def error_dialog(parent, e, title=u'Error'):
    d = wx.MessageDialog(parent, e, title, style=wx.ICON_ERROR|wx.OK)
    d.ShowModal()

# puts dialogs windows into layout?
def into_hbox(f):
    def f1(s, sz):
        # arranges visual elements into vert/hor line
        h = wx.BoxSizer(orient=wx.HORIZONTAL)
        f(s, h)
        sz.Add(h)
    return f1
