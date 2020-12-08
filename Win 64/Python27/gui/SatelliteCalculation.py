# coding=utf-8
import os
# ===============================================================================
# Class containing calculations for satellite parameters and overhead functions
# ===============================================================================

import numpy

import traceback
from exc import LocalError
from gridcalc import *
from satstress import *
seconds_in_year = 31556926.0

class SatelliteCalculation(object):
    """
    Stores the parameters for a satellite, handles saving and loading of satellite and grid files.
    """
    # u'string' == unicode('string')
    # useful for representing more characters than normal ascii can, like other languages

    satellite_vars = [
        ("SYSTEM_ID", u"System ID"),
        ("PLANET_MASS", u"Planet Mass [kg]"),
        ("ORBIT_ECCENTRICITY", u"Orbit Eccentricity"),
        ("ORBIT_SEMIMAJOR_AXIS", u"Orbit Semimajor Axis [m]"),
        ("NSR_PERIOD", u"NSR Period [yrs]")]

    layer_vars_d = [
        ("LAYER_ID", u"Layer ID"),
        ("DENSITY", u"Density [kg/m3]"),
        ("YOUNGS_MODULUS", u"Young's Modulus [Pa]"),
        ("POISSONS_RATIO", u"Poisson's Ratio"),
        ("THICKNESS", u"Thickness [m]"),
        ("VISCOSITY", u"Viscosity [Pa s]")]

    satlayers_d = [
        (3, "ICE_UPPER"),
        (2, "ICE_LOWER"),
        (1, "OCEAN"),
        (0, "CORE")]

    stress_d = {
        u'Nonsynchronous Rotation': NSR,
        u'Diurnal': Diurnal,
        u'Ice Shell Thickening': IST,
        u'Obliquity': DiurnalObliquity,
        u'Polar Wander': PolarWander}

    grid_vars_d = [
        ("MIN", u'Minimum value'),
        ("MAX", u'Maximum value'),
        ("NUM", u'Number of grid points')]

    grid_parameters_d = [
        ("LAT", u'Latitude'),
        ("LON", u'Longitude'),
        ("TIME", u'Time (Periapse = 0)'),
        ("ORBIT", u'Orbital position (Periapse = 0) [°]'),
        ("NSR_PERIOD", u'NSR period'),
        ("POLE_POSITION", u'Initial North Pole Location')]

    cycloid_parameters_d = {
        'YIELD': None,
        'PROPAGATION_STRENGTH': None,
        'PROPAGATION_SPEED': None,
        'STARTING_LATITUDE': None,
        'STARTING_LONGITUDE': None,
        'STARTING_DIRECTION': None,
        'VARY_VELOCITY': None,
        'k': None}

    polarwander_coordinates = {
        'thetaRInitial': None,
        'phiRInitial': None,
        'thetaRFinal': None,
        'phiRFinal': None,
        'thetaTInitial': None,
        'phiTInitial': None,
        'thetaTFinal': None,
        'phiTFinal': None,
        'Locked': True,
        'Despinning': False}

    # These variables are used by the program to determine where to map the initial/final PW points.  -PS 2016

    def __init__(self):
        self.satellite = None
        self.satellite_changed = False
        self.satellite_save_changed = False

        self.stresses = []
        self.stresses_changed = False  # Used to determine whether or not to give a "Save Changes" message.

        self.grid = None
        self.grid_changed = False
        self.grid_save_changed = False

        # Holds all the cycloid objects and their respective parameters. Only used for multiple cycloid loading
        self.cycloids = {}
        params_for_cycloids = {}
        self.many_changed = False  # indicates if new csv file was loaded

        self.cycl_save_changed = False

        self.calc = None
        self.calc_changed = False

        self.cycloid_changed = False

        self.projection_changed = False

        self.parameters = {}
        self.parameters['NSR_PERIOD'] = 'infinity'  # initial value of NSRperiod in Satellite tab
        self.parameters['to_plot_cycloids'] = False
        self.parameters[
            'to_plot_triangles'] = True  # A flag to know whether or not to put in cycloid markers if unable to initiate/propagate.  -PS 2016
        self.parameters[
            'to_plot_pw_markers'] = True  # A flag to know whether or not to plot the PW initial/final points.  -PS 2016
        self.parameters['to_plot_many_cycloids'] = False
        self.parameters['VARY_VELOCITY'] = False
        self.parameters['k'] = 0

        self.cyc = None

    # returns boolean of whether any parameters have been changed
    def changed(self):
        return self.satellite_changed or self.stresses_changed or self.grid_changed

    # returns boolean of whether there are any changes to be saved
    def saveable_changed(self):
        return self.satellite_save_changed or self.grid_save_changed or self.cycl_save_changed

    # sets a parameter to given value and sets changed to True
    def set_parameter(self, parameter, value, point=False):
        if point:
            self.parameters[parameter][point - 1] = value
        else:
            self.parameters[parameter] = value

        if value == 'True':
            self.parameters[parameter] = True
        elif value == 'False':
            self.parameters[parameter] = False

        # if arg parameter is key of attribute satellite_vars
        if parameter in [p for p, d in self.satellite_vars]:
            self.satellite_changed = True
        # if arg parameter is key of sat_layers
        elif parameter in ["%s_%d" % (p, l) for p, d in self.layer_vars_d for l, v in self.satlayers_d]:
            self.satellite_changed = True
            self.satellite_save_changed = True
        # if one of the stress types
        elif self.stress_d.get(parameter):
            self.stresses_changed = True
        # if something to do with the grid
        elif parameter in ["%s_%s" % (p, v) for p, pd in self.grid_parameters_d for v, vd in self.grid_vars_d] + [
            'GRID_ID']:
            self.grid_changed = True
            self.grid_save_changed = True
            if parameter in ['TIME_MAX', 'TIME_MIN', 'TIME_NUM'] and self.parameters.has_key('nsr_time') \
                    and self.parameters.has_key('TIME_NUM'):
                self.parameters['TIME_MAX'] = self.get_parameter(float, 'TIME_MIN', 0) + \
                                              self.get_parameter(float, 'nsr_time', 0) * \
                                              self.get_parameter(float, 'TIME_NUM', 0)
            if not str(value):
                del self.parameters[parameter]

        elif parameter in self.cycloid_parameters_d.keys():
            if not parameter in ['STARTING_DIRECTION', 'VARY_VELOCITY']:
                self.parameters[parameter] = float(value)
            else:
                self.parameters[parameter] = value


        # if NSR related, grid would automatically update?
        elif parameter.startswith('nsr_'):
            self.grid_changed = True
            if parameter == 'nsr_time':
                self.parameters['TIME_MAX'] = self.get_parameter(float, 'TIME_MIN', 0) + self.get_parameter(float,
                                                                                                            'nsr_time',
                                                                                                            0)
            if not str(value):
                del self.parameters[parameter]
        # if projection (has to do with plot) change
        elif parameter == 'projection':
            self.projection_changed = True
        # stress parameter
        elif parameter == 'delta_tc':
            self.stresses_changed = True
            # in km
            self.stress_d['Ice Shell Thickening'].delta_tc = self.get_parameter(float, 'delta_tc', 0) * 1000
        # stress parameter
        elif parameter == 'diffusivity':
            self.stresses_changed = True
            self.stress_d['Ice Shell Thickening'].diffusivity = self.get_parameter(float, 'diffusivity', 0)
        # stress parameter
        elif parameter == 'obliquity':
            self.stresses_changed = True
            self.stress_d['Obliquity'].eps = \
                numpy.radians(self.get_parameter(float, 'obliquity', 0))
        elif parameter == 'periapsis_arg':
            self.stresses_changed = True
            self.stress_d['Obliquity'].periapsis_arg = \
                numpy.radians(self.get_parameter(float, 'periapsis_arg', 0))
        # elif parameter in [k for k,v in self.cycloid_parameters_d]:
        # TODO
        self.cycloid_changed = True

    # accessor for parameters
    def get_parameter(self, f, parameter, default_value=None):
        try:
            return f(self.parameters[parameter])
        except:
            return default_value

    # constructs a satellite from a file if not all parameters are inputted
    def get_satellite(self):
        if self.satellite_changed or not self.satellite:
            self.mk_satellite()
        return self.satellite

    # constructs a satellite object from a file
    def mk_satellite(self):
        filename, tmp = self.save_satellite()
        try:
            self.load_satellite(filename)
        finally:
            os.unlink(filename)
        return self.satellite

    # converts the nsr period from seconds to years
    def nsr_period_seconds2years(self):
        self.parameters['NSR_PERIOD'] = "%.4f" % (float(self.parameters['NSR_PERIOD']) / seconds_in_year)
        if self.parameters['NSR_PERIOD'] == 'inf':
            self.parameters['NSR_PERIOD'] = 'infinity'

    # converts the nsr period from years to seconds
    def nsr_period_years2seconds(self):
        self.parameters['NSR_PERIOD'] = "%.4f" % (float(self.parameters['NSR_PERIOD']) * seconds_in_year)
        if self.parameters['NSR_PERIOD'] == 'inf':
            self.parameters['NSR_PERIOD'] = 'infinity'

    # opens a satellite file and parses it for relevant variables
    def load_satellite(self, filename):
        f = open(filename)
        try:
            self.satellite = Satellite(f)
            for p, v in self.satellite.satParams.items():
                try:
                    self.parameters[p] = v
                except KeyError:
                    pass
            self.satellite_changed = True
            self.satellite_save_changed = False
            self.nsr_period_seconds2years()
        except Exception, e:
            self.satellite = None
            traceback.print_exc()
            raise LocalError(e, u'Satellite Error')
        finally:
            f.close()

    # builds a string with relevant satellite parameters
    def dump_satellite(self):
        try:
            s = "\n".join(["%s = %s" % (v, self.parameters[v]) for v, d in self.satellite_vars])
            s += "\n\n"
            for l, ld in self.satlayers_d:
                for v, d in self.layer_vars_d:
                    p = "%s_%d" % (v, l)
                    s += "%s = %s\n" % (p, self.parameters[p])
                s += "\n"
            return s
        except KeyError, e:
            raise LocalError('Satellite parameter %s is not defined' % str(e), u'Satellite Error')
        except Exception, e:
            raise LocalError(str(e), u'Satellite Error')

    # writes the string created by dump_satellite to a file
    def save_satellite(self, filename=None):
        tmp = False
        if not filename:
            filename = os.tempnam(None,
                                  'sat')  # Saves to the temporary directory.  This is necessary b/c of how satstress.py handles reading in the variables.  -PS 2016
            tmp = True
        f = open(filename, 'w')
        t = self.parameters['NSR_PERIOD']
        self.nsr_period_years2seconds()
        f.write(self.dump_satellite())
        self.parameters['NSR_PERIOD'] = t  # why does it do this? -PS 2016
        f.close()
        if not tmp:
            self.satellite_save_changed = False
        return filename, tmp

    # updates stresses if changed or does not exist
    def get_stresses(self):
        if self.stresses_changed or self.satellite_changed or not self.stresses:
            sat = self.get_satellite()
            self.stresses = [self.stress_d[v](sat) for v in
                             filter(lambda v: self.parameters.get(v), self.stress_d.keys())]
        return self.stresses

    # updates grid or makes one if one does not already exist
    def get_grid(self):
        if self.grid_changed or not self.grid:
            self.mk_grid()
        return self.grid

    # converts years to seconds in the parameters
    def parameter_yrs2secs(self, p):
        v = self.get_parameter(float, p)
        if v:
            self.parameters[p] = "%g" % (v * seconds_in_year)

    # converts seconds to years in the paramters
    def parameter_secs2yrs(self, p):
        v = self.get_parameter(float, p)
        if v:
            self.parameters[p] = "%g" % (float(self.parameters[p]) / seconds_in_year)

    # saves grid parameters to file
    def save_grid(self, filename=None):
        tmp = False
        if filename is None:
            filename = os.tempnam(None, 'grid')
            tmp = True
        f = open(filename, 'w')
        try:
            t_min = self.parameters['TIME_MIN']
            t_max = self.parameters['TIME_MAX']
            self.parameter_yrs2secs('TIME_MIN')
            self.parameter_yrs2secs('TIME_MAX')
        except:
            pass
        f.write(self.dump_grid())
        try:
            self.parameters['TIME_MIN'] = t_min
            self.parameters['TIME_MAX'] = t_max
        except:
            pass
        f.close()
        if not tmp:
            self.grid_save_changed = False
        return filename, tmp

    # opens .grid files generateed by this program
    def load_grid(self, filename):
        f = open(filename)
        try:
            for p, d in self.grid_parameters_d:
                for v, dv in self.grid_vars_d:
                    self.set_parameter("%s_%s" % (p, v), '')
            for p, v in nvf2dict(f).items():
                try:
                    self.set_parameter(p, v)
                except:
                    pass
            try:
                self.parameter_secs2yrs('TIME_MIN')
                self.parameter_secs2yrs('TIME_MAX')
            except:
                pass
            try:
                self.set_parameter('nsr_time',
                                   str(self.get_parameter(float, 'TIME_MAX') - self.get_parameter(float, 'TIME_MIN')))
            except:
                pass
            self.grid_save_changed = False
            self.grid_changed = True
        except Exception, e:
            print e.__class__.__name__, e
            raise LocalError(e, u'Grid Error')
        finally:
            f.close()

    # reads a loaded grid file for relevant parameters
    def mk_grid(self):
        filename, tmp = self.save_grid()
        f = open(filename)
        try:
            self.grid = Grid(f, self.get_satellite())
        except Exception, e:
            raise LocalError(e, u'Grid Error')
        finally:
            f.close()
            os.unlink(filename)
        return self.grid

    # builds the grid for use with save_sgrid
    def dump_grid(self):
        try:
            s = "GRID_ID = %s\n" % self.parameters['GRID_ID']
            self.parameters['NSR_PERIOD_MIN'] = self.parameters['NSR_PERIOD_MAX'] = self.get_satellite().nsr_period
            self.parameters['NSR_PERIOD_NUM'] = 1
            saved = {}
            for rv, rd in self.grid_parameters_d:
                for cv, cd in self.grid_vars_d:
                    p = "%s_%s" % (rv, cv)
                    v = str(self.parameters.get(p, ''))
                    if v and not saved.get(p):
                        s += "%s = %s\n" % (p, v)
                        saved[p] = True
            return s
        except KeyError, e:
            raise LocalError('Grid parameter %s is not defined' % str(e), u'Grid Error')
        except Exception, e:
            raise LocalError(str(e), u'Grid Error')

    # writes calculated love numbers to file
    def save_love(self, filename):
        f = open(filename, 'w')
        for c in self.get_stresses():
            try:
                f.write(str(c))
                f.write("\n")
            except:
                pass
        f.close()

    # updates the calculations and changes the state of self.getstress
    def calculate(self):
        try:
            self.calc = StressCalc(self.get_stresses())
            self.satellite_changed = self.grid_changed = self.stresses_changed = False
            self.calc_changed = True
            return self.calc
        except Exception, e:
            print e.__class__.__name__, str(e)
            if self.satellite and not self.stresses:
                traceback.print_exc()
                raise LocalError(u'Stresses are not defined', u'Calc Error')
            else:
                traceback.print_exc()
                raise LocalError(str(e), u'Calc Error')

    # updates calculations
    def get_calc(self, k=None):
        if self.changed() or self.calc is None:
            self.calculate()
        return self.calc

    # calculates tensor stresses
    def calc_tensor(self, rows=1):
        for i in range(rows):
            theta, phi, t = [float(self.parameters[p][i]) for p in ['theta', 'phi', 't']]
            t *= seconds_in_year
            theta, phi = map(numpy.radians, [theta, phi])
            calc = self.get_calc()
            Ttt, Tpt, Tpp = ["%g" % (x / 1000.) for x in calc.tensor(numpy.pi / 2 - theta, phi, t)]
            self.parameters["Ttt"][i] = Ttt
            self.parameters["Tpt"][i] = Tpt
            self.parameters["Tpp"][i] = Tpp  # http://pirlwww.lpl.arizona.edu/~hoppa/science.html
            s1, a1, s3, a3 = calc.principal_components(numpy.pi / 2 - theta, phi, t)
            self.parameters['s1'][i] = "%g" % (s1 / 1000.)
            self.parameters['s3'][i] = "%g" % (s3 / 1000.)
            self.parameters['a'][i] = "%.2f" % numpy.degrees(a1)

    # saves netcdf files from parameters
    # This function currently won't work because we have commented out the netCDF3 import.
    # The netCDF functionality isn't really necessary, and the netCDF library doesn't play well on Windows. -PS 2016
    def save_netcdf(self, filename):
        try:
            sat = self.get_satellite()
            # unchangable to satisfy write_netcdf
            calc = StressCalc([NSR(sat), Diurnal(sat)])
            grid = self.get_grid()
            gridcalc = GridCalc(grid, calc)
            gridcalc.write_netcdf(filename)
            sat = self.get_satellite()
            nc = netCDF3.Dataset(filename, 'a')
            nc.satellite_nsr_period = sat.nsr_period
            nc.close()
            self.mk_satellite()
            self.satellite_save_changed = self.grid_save_changed = False
        except Exception, e:
            raise LocalError(str(e), 'Export Error')

    # helper function for load_netcdf_satellite
    def write_var(self, fs, nc, v):
        fs.write("%s = %s\n" % (v.upper(), getattr(nc, v)))

    # constructs a netcdf file
    def load_netcdf_satellite(self, nc):
        satfilename = os.tempnam(None, 'sat')
        fs = open(satfilename, 'w')
        for v in ['system_id',
                  'planet_mass',
                  'orbit_eccentricity',
                  'orbit_semimajor_axis']:
            self.write_var(fs, nc, v)
        fs.write("\n")
        for l in range(3, -1, -1):
            for v in ['layer_id',
                      'density',
                      'lame_mu',
                      'lame_lambda',
                      'thickness',
                      'viscosity']:
                self.write_var(fs, nc, "%s_%d" % (v, l))
            fs.write("\n")
        try:
            fs.write("NSR_PERIOD = %e\n" % float(nc.satellite_nsr_period))
        except:
            fs.write("NSR_PERIOD = %e\n" % float(nc.variables['nsr_period'][0]))
        fs.close()
        try:
            self.load_satellite(satfilename)
        except:
            pass
        os.unlink(satfilename)

    # helper function for load_netcdf_grid
    def write_var2(self, fs, nc, v, fv):
        fs.write("%s_MIN = %g\n" % (fv, float(nc.variables[v][0])))
        try:
            fs.write("%s_MAX = %g\n" % (fv, float(nc.variables[v][-1])))
        except:
            fs.write("%s_MAX = %g\n" % (fv, float(nc.variables[v][0])))
        fs.write("%s_NUM = %d\n" % (fv, len(nc.dimensions[v])))

    # writes netcdf grids with write_var
    def load_netcdf_grid(self, nc):
        gridfilename = os.tempnam(None, 'grid')
        fs = open(gridfilename, 'w')
        self.write_var(fs, nc, 'grid_id')
        for v, fv in [
            ('latitude', 'LAT'),
            ('longitude', 'LON'),
            ('nsr_period', 'NSR_PERIOD')]:
            self.write_var2(fs, nc, v, fv)
        if nc.variables['time'].units == 'seconds':
            self.write_var2(fs, nc, 'time', 'TIME')
        elif nc.variables['time'].units == 'degrees':
            self.write_var2(fs, nc, 'time', 'ORBIT')
        fs.close()
        try:
            self.load_grid(gridfilename)
        except:
            pass
        os.unlink(gridfilename)

    def load_netcdf(self, filename):
        nc = netCDF3.Dataset(filename, 'r')
        self.load_netcdf_satellite(nc)
        self.load_netcdf_grid(nc)
        self.set_parameter('Nonsynchronous Rotation', True)
        self.set_parameter('Diurnal', True)
        self.satellite_save_changed = self.grid_save_changed = False
        nc.close()
