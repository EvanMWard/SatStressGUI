from Dialog import dir_dialog, error_dialog
from exc import LocalError
from MatPlotPanel import *
from StepSlider import *

scale_left = 0.25
scale_bar_length = 0.38
vector_mult = 4000

class StressPlotPanel(MatPlotPanel):
    """
    Contains controls for going through the time frame dictated in "Grid" Tab.
    Specifically, the [< | < | > | >] controls
    """
    scale_y    = 0.15
    orbit_y    = 0.11
    polar_y = 0.01
    nsr_y    = 0.06
    button_l = 0.04
    bbutton_l= 0.12
    slider_h = 0.04
    slider_x = scale_left + scale_bar_length + button_l*2

    def __init__(self, *args, **kw):
        super(StressPlotPanel, self).__init__(*args, **kw)
        self.figure.subplots_adjust(bottom=0.25)
        # creates scale bar for the vectors (arrows) i.e. |-----| 91 kPa
        self.scale_ax = self.figure.add_axes([scale_left, self.scale_y, scale_bar_length, self.slider_h], frame_on=False)
        self.add_orbit()
        #self.add_polar()
        self.add_nsr()

    def get_ax_orbit(self):
        return self.figure.add_axes([scale_left, self.orbit_y, scale_bar_length, self.slider_h])

    def get_ax_nsr(self):
        return self.figure.add_axes([scale_left, self.nsr_y, scale_bar_length, self.slider_h])

    def get_ax_polar(self):
        return self.figure.add_axes([scale_left, self.polar_y, scale_bar_length, self.slider_h])

    def del_orbit(self):
        self.figure.delaxes(self.ax_orbit)
        self.del_orbit_controls()

    def del_nsr(self):
        self.figure.delaxes(self.ax_nsr)
        self.del_nsr_controls()

    def del_polar(self):
        self.figure.delaxes(self.ax_polar)
        self.del_polar_controls()

    def del_orbit_controls(self):
        for a in [self.ax_orbit_first, self.ax_orbit_prev, \
            self.ax_orbit_next, self.ax_orbit_last, self.ax_orbit_save]:
            self.figure.delaxes(a)

    def del_nsr_controls(self):
        for a in [self.ax_nsr_first, self.ax_nsr_prev, \
            self.ax_nsr_next, self.ax_nsr_last, self.ax_nsr_save]:
            self.figure.delaxes(a)

    def del_polar_controls(self):
        for a in [self.ax_polar_first, self.ax_polar_prev, \
            self.ax_polar_next, self.ax_polar_last, self.ax_polar_save]:
            self.figure.delaxes(a)

    def add_orbit(self):
        self.ax_orbit = self.get_ax_orbit()
        self.add_orbit_controls()

    def add_orbit_controls(self):
        x = self.slider_x
        self.ax_orbit_first = self.figure.add_axes([x, self.orbit_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_orbit_prev = self.figure.add_axes([x, self.orbit_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_orbit_next = self.figure.add_axes([x, self.orbit_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_orbit_last = self.figure.add_axes([x, self.orbit_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_orbit_save = self.figure.add_axes([x, self.orbit_y, self.bbutton_l, self.slider_h])

        # Note: StepSlider is custom designed class in/for gui
        self.orbit_slider = StepSlider(self.ax_orbit, 'Orbital position', 0, 1, 10, valinit=0, dragging=False)

        self.orbit_first_button = matplotlib.widgets.Button(self.ax_orbit_first, '[<')
        self.orbit_prev_button = matplotlib.widgets.Button(self.ax_orbit_prev, '<')
        self.orbit_next_button = matplotlib.widgets.Button(self.ax_orbit_next, '>')
        self.orbit_last_button = matplotlib.widgets.Button(self.ax_orbit_last, '>]')
        self.orbit_save_button = matplotlib.widgets.Button(self.ax_orbit_save, 'Save series')

        self.orbit_first_button.on_clicked(lambda e: self.orbit_slider.first())  # lambda functions
        self.orbit_prev_button.on_clicked(lambda e: self.orbit_slider.prev())    # ok, so are empirically necessary, but why?
        self.orbit_next_button.on_clicked(lambda e: self.orbit_slider.next())
        self.orbit_last_button.on_clicked(lambda e: self.orbit_slider.last())
        # hack
        self.orbit_save_button.on_clicked(lambda e: wx.CallLater(125, self.on_save_orbit_series, e))

    def add_polar(self):
        self.ax_polar = self.get_ax_polar()
        self.add_polar_controls()

    def add_polar_controls(self):
        x = self.slider_x
        self.ax_polar_first = self.figure.add_axes([x, self.polar_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_polar_prev = self.figure.add_axes([x, self.polar_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_polar_next = self.figure.add_axes([x, self.polar_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_polar_last = self.figure.add_axes([x, self.polar_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_polar_save = self.figure.add_axes([x, self.polar_y, self.bbutton_l, self.slider_h])

        # Note: StepSlider is custom designed class in/for gui
        self.polar_slider = StepSlider(self.ax_polar, 'Polar position', 0, 1, 10, valinit=0, dragging=False)

        self.polar_first_button = matplotlib.widgets.Button(self.ax_polar_first, '[<')
        self.polar_prev_button = matplotlib.widgets.Button(self.ax_polar_prev, '<')
        self.polar_next_button = matplotlib.widgets.Button(self.ax_polar_next, '>')
        self.polar_last_button = matplotlib.widgets.Button(self.ax_polar_last, '>]')
        self.polar_save_button = matplotlib.widgets.Button(self.ax_polar_save, 'Save series')

        self.polar_first_button.on_clicked(lambda e: self.polar_slider.first())
        self.polar_prev_button.on_clicked(lambda e: self.polar_slider.prev())
        self.polar_next_button.on_clicked(lambda e: self.polar_slider.next())
        self.polar_last_button.on_clicked(lambda e: self.polar_slider.last())
        self.polar_save_button.on_clicked(lambda e: wx.CallLater(125, self.on_save_polar_series, e))

    def add_nsr(self):
        self.ax_nsr = self.get_ax_nsr()
        self.add_nsr_controls()

    def add_nsr_controls(self):
        x = self.slider_x
        self.ax_nsr_first = self.figure.add_axes([x, self.nsr_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_nsr_prev = self.figure.add_axes([x, self.nsr_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_nsr_next = self.figure.add_axes([x, self.nsr_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_nsr_last = self.figure.add_axes([x, self.nsr_y, self.button_l, self.slider_h])
        x += self.button_l
        self.ax_nsr_save = self.figure.add_axes([x, self.nsr_y, self.bbutton_l, self.slider_h])
        self.nsr_slider = StepSlider(self.ax_nsr, 'NSR position', 0, 1, 10, valinit=0, dragging=False, valfmt="%.1g")
        self.nsr_first_button = matplotlib.widgets.Button(self.ax_nsr_first, '[<')
        self.nsr_prev_button = matplotlib.widgets.Button(self.ax_nsr_prev, '<')
        self.nsr_next_button = matplotlib.widgets.Button(self.ax_nsr_next, '>')
        self.nsr_last_button = matplotlib.widgets.Button(self.ax_nsr_last, '>]')
        self.nsr_first_button.on_clicked(lambda e: self.nsr_slider.first())
        self.nsr_prev_button.on_clicked(lambda e: self.nsr_slider.prev())
        self.nsr_next_button.on_clicked(lambda e: self.nsr_slider.next())
        self.nsr_last_button.on_clicked(lambda e: self.nsr_slider.last())
        self.nsr_save_button = matplotlib.widgets.Button(self.ax_nsr_save, 'Save series')
        self.nsr_save_button.on_clicked(lambda e: wx.CallLater(125, self.on_save_nsr_series, e))

    def change_slider(self, ax, slider, label=None, valmin=None, valmax=None, numsteps=None, valinit=None, valfmt=None):
        if label is None:
            label = slider.label.get_text()
        if valmin is None:
            valmin = slider.valmin
        if valmax is None:
            valmax = slider.valmax
        if numsteps is None:
            numsteps = slider.numsteps
        if valinit is None:
            valinit = slider.valinit
        if valfmt is None:
            valfmt = slider.valfmt
        f = slider.on_changed_f
        slider = StepSlider(ax, label, valmin, valmax, numsteps, valinit=valinit, dragging=False, valfmt=valfmt)
        slider.on_changed(f)
        return slider

    def change_orbit_slider(self, valmin, valmax, numsteps, valinit=None):
        if valinit is None:
            valinit = valmin
        self.figure.delaxes(self.ax_orbit)
        self.ax_orbit = self.get_ax_orbit()
        self.orbit_slider = self.change_slider(
            self.ax_orbit, self.orbit_slider, valmin=valmin, valmax=valmax, numsteps=numsteps, valinit=valinit)

    def change_nsr_slider(self, valmin, valmax, numsteps, valinit=None):
        if valinit is None:
            valinit = valmin
        self.figure.delaxes(self.ax_nsr)
        self.ax_nsr = self.get_ax_nsr()
        self.nsr_slider = self.change_slider(
            self.ax_nsr, self.nsr_slider, valmin=valmin, valmax=valmax, numsteps=numsteps, valinit=valmin, valfmt="%.1g")

    def change_polar_slider(self, valmin, valmax, numsteps, valinit=None):
        if valinit is None:
            valinit = valmin
        self.figure.delaxes(self.ax_polar)
        self.ax_polar = self.get_ax_polar()
        self.polar_slider = self.change_slider(self.ax_polar, self.polar_slider, valmin=valmin, valmax=valmax, numsteps=numsteps, valinit=valmin, valfmt="%.1g")

    def plot_scale(self, scale, valfmt):
        self.scale_ax.clear()
        while self.scale_ax.texts:
            self.scale_ax.texts.pop()
        self.scale_ax.set_xticks([])
        self.scale_ax.set_yticks([])
        self.scale_ax.text(-0.02, 0.5, 'Scale', transform=self.scale_ax.transAxes, va='center', ha='right')
        self.scale_ax.text(0.23, 0.5, valfmt % scale, transform=self.scale_ax.transAxes, va='center', ha='left')
        self.scale_ax.plot([0.00, 0.20], [0.5, 0.5], linestyle='solid', marker='|', color='black', lw=1)
        self.scale_ax.set_xlim(0.0, 1.0)

    def on_save_orbit_series(self, evt):
        try:
            dir_dialog(None,
                message=u"Save calculation series on orbit period",
                style=wx.SAVE,
                action=self.save_orbit_series)
        except LocalError, e:
            error_dialog(self, str(e), e.title)

    def on_save_nsr_series(self, evt):
        try:
            dir_dialog(self,
                message=u"Save calculation series on nsr period",
                style=wx.SAVE,
                action=self.save_nsr_series)
        except LocalError, e:
            error_dialog(self, str(e), e.title)

    def on_save_polar_series(self, evt):
        try:
            dir_dialog(self,
                       message=u"Save calculation series on nsr period",
                       style=wx.SAVE,
                       action=self.save_polar_series)
        except LocalError, e:
            error_dialog(self, str(e), e.title)