import wx
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
from CustomPlotToolbar import *

display_dpi = 72
class MatPlotPanel(wx.Panel):
    """
    GUI object that holds the plot area
    """
    def __init__(self, *args, **kw):
        super(MatPlotPanel, self).__init__(*args, **kw)
        self.figure = Figure(figsize=(6,5),dpi=display_dpi)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.ax = self.figure.add_subplot(111)

        toolbar = CustomPlotToolbar(self.canvas)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, flag=wx.EXPAND|wx.ALL)
        toolbar.Realize()
        sizer.Add(toolbar, flag=wx.EXPAND|wx.ALL)
        self.SetSizer(sizer)
        self.SetMinSize((625, 400))

    def get_axes(self):
        return self.ax

    def draw(self):
        self.canvas.draw()

    def colorbar(self, mappable, *a, **kw):
        #return self.figure.colorbar(mappable, ax=self.ax, *a, **kw)
        return self.figure.colorbar(mappable, *a, **kw)