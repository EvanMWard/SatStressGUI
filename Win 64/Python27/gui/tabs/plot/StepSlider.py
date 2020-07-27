import matplotlib


class StepSlider(matplotlib.widgets.Slider):
    """
    Custom designed class for discrete slider control at bottom of plot panel to control
    satellite's orbital position.

    Used in add_orbit_controls and add_nsr_controls
    """
    def __init__(self, ax, label, valmin, valmax, numsteps, *args, **kw):
        self.steps_n = numsteps
        self.updating = False
        self.prev_val = kw.get('valinit', 0)
        matplotlib.widgets.Slider.__init__(self, ax, label, valmin, valmax, *args, **kw)
        ax.lines.remove(self.vline)

    def on_changed_f(self, val):
        pass

    def on_changed(self, f):
        def f2(val):
            if self.updating:
                return
            self.eventson = False
            self.updating = True
            val += self.valmin
            self.set_stepval(val)
            f(self.val)
            self.updating = False
            self.eventson = True
        self.on_changed_f = f
        matplotlib.widgets.Slider.on_changed(self, f2)

    def set_stepval(self, val):
        if val < self.valmin:
            self.set_val(self.valmin)
        elif val > self.valmax:
            self.set_val(self.valmax)
        elif self.valmax - self.valmin > 0 and self.steps_n > 0:
            step = float(self.valmax - self.valmin)/self.steps_n
            n0 = int((val - self.valmin)/step)
            n1 = n0 + 1
            if abs(val - self.prev_val) > 0.7*step:
                self.prev_val = round((val - self.valmin)/step)*step + self.valmin
                self.set_val(self.prev_val)
            else:
                self.set_val(self.prev_val)

    def reset(self):
        self.updating = True
        matplotlib.widgets.Slider.reset(self)
        self.updating = False

    def first(self):
        self.set_stepval(self.valmin)
        # HERE initial_split()
        # ^ I think that these comments have to do with cycloid generation -PS 2016

    def last(self):
        self.set_stepval(self.valmax)
        # HERE EVERYTHING SHOULD BE GRAPHED

    def next(self):
        step = float(self.valmax - self.valmin)/self.steps_n
        n = int((self.val - self.valmin)/step) + 1
        self.set_stepval(n*step + self.valmin)
        # ONLY GRAPH UP TO THIS POINT

    def prev(self):
        step = float(self.valmax - self.valmin)/self.steps_n
        n = int((self.val - self.valmin)/step) - 1
        self.set_stepval(n*step + self.valmin)
        # ONLY GRAPH UP TO THIS POINT