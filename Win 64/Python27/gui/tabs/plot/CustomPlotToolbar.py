from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
class CustomPlotToolbar(NavigationToolbar):
    def __init__(self, plotCanvase):
        # create default toolbar
        NavigationToolbar.__init__(self, plotCanvase)

        # remove unwanted button
        # stress plot only exists in rectangular bounds
        # may need to add in later if want movable scale bar
        # or only pan when zoommed
        POSITION_OF_PANNING_BUTTON = 3

        # remove unnecessary button (no subplots)
        POSITION_OF_CONFIGURE_SUBPLOT_BUTTON = 6
        self.DeleteToolByPos(POSITION_OF_CONFIGURE_SUBPLOT_BUTTON)