# ===============================================================================
# Class containing overhead functions for configuring, used in PlotPanel
# ===============================================================================
class Config:
    """
    Class the holds application settings --> specifically?
    """
    #default step for plots <-- what units?
    default_step = 30

    def __init__(self, configfile='config'):
        self.configfile = configfile
        self.conf = {}

    # a is optional arg
    def load(self, *a):
        try:
            c = open(self.configfile)
            self.conf = nvf2dict(c)
            c.close()
            ret = filter(lambda x: x, map(self.conf.get, a))
            if len(a) == 1 and len(ret) == 1:
                return ret[0]
            else:
                return ret
        except:
            self.conf = {}
            return []

    # **kw unpacks the extra dictionary args
    def save(self, **kw):
        for k, v in kw.items():
            self.conf[k] = v   # conf is dictionary
        try:
            c = open(self.configfile, 'w')
            c.writelines([ "%s = %s\n" % (k,v) for k,v in self.conf.items() ])
            c.close()
        except:
            pass

    def load_step(self, step_field='STEP'):
        self.load()
        try:
            return float(self.conf.get(step_field, self.default_step))
        except:
            return self.default_step

    def save_step(self, step, step_field='STEP'):
        self.conf[step_field] = step
        self.save()

#creates a global instance of config
config = Config()