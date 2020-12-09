"""
A framework for calculating the surface stresses at a particular place and time
on a satellite resulting from one or more tidal potentials.

1 Input and Output
==================
Because C{satstress} is a "library" module, it doesn't do a lot of input and
output - it's mostly about doing calculations.  It does need to read in the
specification of a L{Satellite} object though, and it can write the same kind
of specification out.  To do this, it uses name-value files, and a function
called L{nvf2dict}, which creates a Python dictionary (or "associative array").

A name-value file is just a file containing a bunch of name-value pairs, like::

  ORBIT_ECCENTRICITY = 0.0094   # e must be < 0.25

It can also contain comments to enhance human readability (anything following a
'#' on a line is ignored, as with the note in the line above).

2 Satellites
============
Obviously if we want to calculate the stresses on the surface of a satellite,
we need to define the satellite, this is what the L{Satellite} object does.

2.1 Specifying a Satellite
--------------------------
  In order to specify a satellite, we need:

    - an ID of some kind for the planet/satellite pair of interest
    - the charactaristics of the satellite's orbital environment
    - the satellite's internal structure and material properties
    - the forcings to which the satellite is subjected

  From a few basic inputs, we can calculate many derived characteristics, such
  as the satellite's orbital period or the surface gravity.

  The internal structure and material properties are specified by a series of
  concentric spherical shells (layers), each one being homogeneous throughout
  its extent.  Given the densities and thicknesses of the these layers, we can
  calculate the satellite's overall size, mass, density, etc.

  Specifying a tidal forcing may be simple or complex.  For instance, the
  L{Diurnal} forcing depends only on the orbital eccentricity (and other
  orbital parameters already supplied), and the L{NSR} forcing requires only
  the addition of the non-synchronous rotation period of the shell.  Specifying
  an arbitrary true polar wander trajectory would be much more complex.

  For the moment, becuase we are only including simple forcings, their
  specifying parameters are read in from the satellite definition file.  If
  more, and more complex forcings are eventually added to the model, their
  specification will probably be split into a separate input file.

2.2 Internal Structure and Love Numbers
---------------------------------------
  C{satstress} treats the solid portions of the satellite as U{viscoelastic
  Maxwell solids <http://en.wikipedia.org/wiki/Maxwell_material>}, that respond
  differently to forcings having different frequencies (S{omega}).  Given the a
  specification of the internal structure and material properties of a
  satellite as a series of layers, and information about the tidal forcings the
  body is subject to, it's possible to calculate appropriate Love numbers,
  which describe how the body responds to a change in the gravitational
  potential.

  Currently the calculation of Love numbers is done by an external program
  written in Fortran by John Wahr and others, with roots reaching deep into the
  Dark Ages of computing.  As that code (or another Love number code) is more
  closely integrated with the model, the internal structure of the satellite
  will become more flexible, but for the moment, we are limited to assuming
  a 4-layer structure:

    - B{C{ICE_UPPER}}: The upper portion of the shell (cooler, stiffer)
    - B{C{ICE_LOWER}}: The lower portion of the shell (warmer, softer)
    - B{C{OCEAN}}: An inviscid fluid decoupling the shell from the core.
    - B{C{CORE}}: The silicate interior of the body.

3 Stresses
==========

C{satstress} can calculate the following stress fields:

  1. B{L{Diurnal}}: stresses arising from an eccentric orbit, having a
     forcing frequency equal to the orbital frequency.

  2. B{L{NSR}}: stresses arising due to the faster-than-synchronous rotation
     of a floating shell that is decoupled from the satellite's interior by a
     fluid layer (an ocean).

The expressions defining these stress fields are derived in "Modeling Stresses
on Satellites due to Non-Synchronous Rotation and Orbital Eccentricity Using
Gravitational Potential Theory" (U{preprint, 15MB PDF
<http://satstress.googlecode.com/files/Wahretal2008.pdf>}) by Wahr et al.
(submitted to I{Icarus}, in March, 2008).

3.1 Stress Fields Live in L{StressDef} Objects
----------------------------------------------
  Each of the above stress fields is defined by a similarly named L{StressDef}
  object.  These objects contain the formulae necessary to calculate the
  surface stress.  The expressions for the stresses depend on many parameters
  which are defined within the L{Satellite} object, and so to create a
  L{StressDef} object, you need to provide a L{Satellite} object.

  There are many formulae which are identical for both the L{NSR} and
  L{Diurnal} stress fields, and so instead of duplicating them in both classes,
  they reside in the L{StressDef} I{base class}, from which all L{StressDef}
  objects inherit many properties.

  The main requirement for each L{StressDef} object is that it must define the
  three components of the stress tensor S{tau}:

    - C{Ttt} (S{tau}_S{theta}S{theta}) the north-south (latitudinal) component

    - C{Tpt} (S{tau}_S{phi}S{theta} = S{tau}_S{theta}S{phi}) the shear
      component

    - C{Tpp} (S{tau}_S{phi}S{phi}) the east-west (longitudinal) component

3.2 Stress Calculations are Performed by L{StressCalc} Objects
--------------------------------------------------------------
  Once you've I{instantiated} a L{StressDef} object, or several of them (one
  for each stress you want to include), you can compose them together into a
  L{StressCalc} object, which will actually do calculations at given points on
  the surface, and given times, and return a 2x2 matrix containing the
  resulting stress tensor (each component of which is the sum of all of the
  corresponding components of the stress fields that were used to instantiated
  the L{StressCalc} object).

  This is (hopefully) easier than it sounds.  With the following few lines, you
  can construct a satellite, do a single calculation on its surface, and see
  what it looks like:

  >>> from satstress.satstress import *
  >>> the_sat = Satellite(open("input/Europa.satellite"))
  >>> the_stresses = StressCalc([Diurnal(the_sat), NSR(the_sat)])
  >>> Tau = the_stresses.tensor(theta=pi/4.0, phi=pi/3.0, t=10000)
  >>> print(Tau)

  The C{test} program included in the satstress distribution shows a slightly
  more complex example, which should be enough to get you started using the
  package.

3.3 Extending the Model  
-----------------------
  Other stress fields can (and hopefully will!), be added easily, so long as
  they use the same mathematical definition of the membrane stress tensor
  (S{tau}), as a function of co-latitude (S{theta}) (measured south from the
  north pole), east-positive longitude (S{phi}), measured from the meridian
  on the satellite which passes through the point on the satellite directly
  beneath the parent planet (assuming a synchronously rotating satellite),
  and time (B{M{t}}), defined as seconds elapsed since pericenter.

  This module could also potentially be extended to also calculate the
  surface strain (S{epsilon}) and displacement (B{M{s}}) fields, or to
  calculate the stresses at any point within the satellite.

  @group Exceptions (error handling classes): *Error
"""

# regular expressions
import re

# required for command line parsing
import sys

# Scientific functions... like complex exponentials
import scipy, scipy.special
import numpy

# Physical constants, like the Newton's Gravitational Constant Big G
import physcon as pc

# New Love number calculator
import thirdLove as love3

##############################
#      HELPER FUNCTIONS      #
##############################
def nvf2dict(nvf, comment='#'):
    """
    Reads from a file object listing name value pairs, creating and returning a
    corresponding Python dictionary.

    The file should contain a series of name value pairs, one per line
    separated by the '=' character, with names on the left and values on the
    right.  Blank lines are ignored, as are lines beginning with the comment
    character (assumed to be the pound or hash character '#', unless otherwise
    specified).  End of line comments are also allowed.  String values should
    not be quoted in the file.  Names are case sensitive.

    Returns a Python dictionary that uses the names as keys and the values as
    values, and so all Python limitations on what can be used as a dictionary
    key apply to the name fields.

    Leading and trailing whitespace is stripped from all names and values, and
    all values are returned as strings.

    @param nvf: an open file object from which to read the name value pairs
    @param comment: character which begins comments
    @type nvf: file
    @type comment: str
    @return: a dictionary containing the name value pairs read in from C{nvf}.
    @rtype: dict
    @raise NameValueFileParseError: if a non-comment input line does not
    contain an '=' character, or if a non-comment line has nothing but
    whitespace preceeding or following the '=' character.
    @raise NameValueFileDuplicateNameError: if more than one instance of the
    same name is found in the input file C{nvf}.

    """

    params = {}
    for line in nvf:

        line = line.strip()
        # strip comments
        line = re.sub("\s*%s.*" % (comment), "", line)
        # if the line is all whitespace, ignore it
        if re.match("^\s*$", line):
            continue

        # if we don't have two non-space strings separated by '=', Error!
        if not (re.match("[^\s]*.*=.*[^\s]", line)):
            raise NameValueFileParseError(nvf, line)
        # split the NAME=VALUE pair into two variables if possible.
        try:
            (name,value) = re.split("=",line,1)
        except ValueError:
            raise NameValueFileParseError(nvf, line)

        # remove leading and trailing whitespace
        name  = name.strip()
        value = value.strip()

        # If the file we're parsing has the same name specified more than once,
        # then we're going to end up clobbering the earlier value here, which
        # is probably not what we want to do.
        if params.has_key(name):
            raise NameValueFileDuplicateNameError(nvf, name)
        else:
            params[name] = value
    return params
# end nvf2dict()

def eigen2(a, b, c):
    """
    Calculate the normalized eigenvectors and eigenvalues for a 2x2 symmetric
    matrix, where a, b, and c are:

    A   B

    B   C

    """
    sqrt_thing = numpy.sqrt(abs(a*a + 4.0*b*b - 2.0*a*c + c*c))

    lambda1   = (0.5)*(a + c - sqrt_thing)
    eig1theta = (-a + c + sqrt_thing)/(-2.0*b)
    eig1phi   = 1.0

    lambda2   = (0.5)*(a + c + sqrt_thing)
    eig2theta = (-a + c - sqrt_thing)/(-2.0*b)
    eig2phi   = 1.0

    mag1 = numpy.sqrt(eig1theta**2 + eig1phi**2)
    mag2 = numpy.sqrt(eig2theta**2 + eig2phi**2)

    eig1theta /= mag1
    eig1phi   /= mag1
    eig2theta /= mag2
    eig2phi   /= mag2

    return(numpy.array([lambda1, eig1theta, eig1phi, lambda2, eig2theta, eig2phi]))

##############################
#          CLASSES           #
##############################
class Satellite(object):
    """An object describing the physical structure and context of a satellite.
    
    Defines a satellite's material properties, internal structure, orbital
    context, and the tidal forcings to which it is subjected.

    @ivar sourcefilename: the file object which was read in order to create
    the L{Satellite} instance.
    @type sourcefilename: str
    @ivar satParams: dictionary containing the name value pairs read in from
    the input file.
    @type satParams: dict
    @ivar system_id: string identifying the planet/satellite system,
    corresponds to C{SYSTEM_ID} in the input file.
    @type system_id: str
    @ivar orbit_semimajor_axis: semimajor axis of the satellite's orbit [m],
    corresponds to C{ORBIT_SEMIMAJOR_AXIS} in the input file.
    @type orbit_semimajor_axis: float
    @ivar orbit_eccentricity: the satellite's orbital eccentricity,
    corresponds to C{ORBIT_ECCENTRICITY} in the input file.
    @type orbit_eccentricity: float
    @ivar planet_mass: the mass of the planet the satellite orbits [kg],
    corresponds to C{PLANET_MASS} in the input file.
    @type planet_mass: float
    @ivar nsr_period: the time it takes for the decoupled ice shell to
    complete one full rotation [s], corresponds to C{NSR_PERIOD} in the input
    file.  May also be set to infinity (inf, infinity, INF, INFINITY).
    @type nsr_period: float
    @ivar num_layers: the number of layers making up the satellite, as
    indicated by the number of keys within the satParams dictionary contain
    the string 'LAYER_ID'.  Currently this must equal 4.
    @type num_layers: int
    @ivar layers: a list of L{SatLayer} objects, describing the layers making
    up the satellite.  The layers are ordered from the center of the satellite
    outward, with layers[0] corresponding to the core.
    @type layers: list
    """

    def __init__(self, satFile=None): #
        """
        Construct a Satellite object from a satFile

        Required input file parameters:
        ===============================
        The Satellite is initialized from a name value file (as described under
        L{nvf2dict}).  The file must define the following parameters, all of which
        are specified in SI (MKS) units.

          - B{C{SYSTEM_ID}}:  A string identifying the planetary system, e.g.
            JupiterEuropa.

          - B{C{PLANET_MASS}}:  The mass of the planet the satellite orbits
            [kg].

          - B{C{ORBIT_ECCENTRICITY}}:  The eccentricity of the satellite's
            orbit.  Must not exceed 0.25.

          - B{C{ORBIT_SEMIMAJOR_AXIS}}:  The semimajor axis of the satellite's
            orbit [m].

          - B{C{NSR_PERIOD}}:  The time it takes for the satellite's icy shell
            to undergo one full rotation [s].  If you don't want to have any
            NSR stresses, just put INFINITY here.

        An example satFile is included with the satstress package::

            satstress-X.Y.Z/input/Europa.satellite

        Where X.Y.Z refers to the version numbers.

        @param satFile: Open file object containing name value pairs specifying
        the satellite's internal structure and orbital context, and the tidal
        forcings to which the satellite is subjected.
        @type satFile: file
        @return: a Satellite object corresponding to the proffered input file.
        @rtype: L{Satellite}

        @raise NameValueFileError: if parsing of the input file fails.
        @raise MissingSatelliteParamError: if a required input field is not
        found within the input file.
        @raise NonNumberSatelliteParamError: if a required input which is of
        a numeric type is found within the file, but its value is not
        convertable to a float.
        @raise LoveLayerNumberError: if the number of layers specified is not
        exactly 4.
        @raise GravitationallyUnstableSatelliteError: if the layers are found not
        to decrease in density from the core to the surface.
        @raise ExcessiveSatelliteMassError: if the mass of the satellite's parent
        planet is not at least 10 times larger than the mass of the satellite.
        @raise LargeEccentricityError: if the satellite's orbital eccentricity is
        greater than 0.25
        @raise NegativeNSRPeriodError: if the NSR period of the satellite is less
        than zero.
        """

        # Make a note of the input file defining this satellite:
        self.sourcefilename = satFile.name

        # Slurp the proffered parameters into a dictionary:
        try:
            self.satParams = nvf2dict(satFile, comment='#')
        except NameValueFileError:
            raise

        # Try to assign our data members to the values we read in from
        # the file.  If we didn't get all the values we need, raise an
        # exception saying so.

        try:
            self.system_id = self.satParams['SYSTEM_ID']
        except KeyError:
            raise MissingSatelliteParamError(self, 'SYSTEM_ID')

        try:
            self.orbit_semimajor_axis = float(self.satParams['ORBIT_SEMIMAJOR_AXIS'])
        except KeyError:
            raise MissingSatelliteParamError(self, 'ORBIT_SEMIMAJOR_AXIS')
        except ValueError:
            raise NonNumberSatelliteParamError(self, 'ORBIT_SEMIMAJOR_AXIS')

        try:
            self.orbit_eccentricity = float(self.satParams['ORBIT_ECCENTRICITY'])
        except KeyError:
            raise MissingSatelliteParamError(self, 'ORBIT_ECCENTRICITY')
        except ValueError:
            raise NonNumberSatelliteParamError(self, 'ORBIT_SEMIMAJOR_AXIS')

        try:
            self.planet_mass = float(self.satParams['PLANET_MASS'])
        except KeyError:
            raise MissingSatelliteParamError(self, 'PLANET_MASS')
        except ValueError:
            raise NonNumberSatelliteParamError(self, 'PLANET_MASS')

        try:
            self.nsr_period = float(self.satParams['NSR_PERIOD'])
        except KeyError:
            raise MissingSatelliteParamError(self, 'NSR_PERIOD')
        except ValueError:
            raise NonNumberSatelliteParamError(self, 'NSR_PERIOD')

        

        # Figure out how many layers are listed in the input file, and give the
        # user an error if it's not exactly 4 (the current Love number code can
        # only deal with 4):
        self.num_layers = 0
        for key in self.satParams.keys():
            if re.match('LAYER_ID', key): self.num_layers += 1
    
        if self.num_layers != 4:
            raise LoveLayerNumberError(self)
 
        # construct the layer objects:
        self.layers = []
        n = 0
        while (n < self.num_layers):
            self.layers.append(SatLayer(self,layer_n=n))

            # Verify that layers increase in density toward the center of the
            # satellite.  If they don't the satellite will be gravitationally
            # unstable, and the Love number code will blow up:
            if (n > 0):
                if (self.layers[n].density > self.layers[n-1].density): 
                    raise GravitationallyUnstableSatelliteError(self, n)
            n = n+1

        if self.planet_mass < 2*self.mass():
            raise ExcessiveSatelliteMassError(self)

        if self.orbit_eccentricity > 0.25:
            raise LargeEccentricityError(self)

        if self.nsr_period < 0:
            raise NegativeNSRPeriodError(self)

        self.userDefLove = False

    # end __init__

    # Methods defining derived quantities: 
    def mass(self):
        """Calculate the mass of the satellite. (the sum of the layer masses)"""
        mass   = 0.0
        radius = 0.0
        for layer in self.layers:
            mass += ((4.0/3.0)*scipy.pi)*((radius+layer.thickness)**3 - (radius**3))*layer.density
            radius += layer.thickness
        return(mass)

    def radius(self):
        """Calculate the radius of the satellite (the sum of the layer thicknesses)."""
        radius = 0.0
        for layer in self.layers:
            radius += layer.thickness
        return(radius)

    def density(self):
        """Calculate the mean density of the satellite in [kg m^-3]."""
        return(self.mass() / ((4.0/3.0)*scipy.pi*(self.radius()**3)))

    def surface_gravity(self):
        """Calculate the satellite's surface gravitational acceleration in [m s^-2]."""
        return(pc.G*self.mass() / (self.radius()**2))

    def orbit_period(self):
        """Calculate the satellite's Keplerian orbital period in seconds."""
        return(2.0*scipy.pi*scipy.sqrt( (self.orbit_semimajor_axis**3) / (pc.G*self.planet_mass) ))

    def mean_motion(self):
        """Calculate the orbital mean motion of the satellite [rad s^-1]."""
        return(2.0*scipy.pi / (self.orbit_period()))

    #  end derived quantities

    def __str__(self): #
        """Output a satellite definition file equivalent to the object."""

        myStr = """#
# Satellite system definition file for use with the Python satstress package.
# All quantities are in SI (meters, kilograms, seconds) units.  For more
# information, see:
#
# http://code.google.com/p/satstress

############################################################
# Basic Satellite parameters:
############################################################

# Satellite system name
SYSTEM_ID = %s

# System orbital parameters:
PLANET_MASS          = %g
ORBIT_ECCENTRICITY   = %g
ORBIT_SEMIMAJOR_AXIS = %g

# Additional parameters required to calculate Love numbers:
NSR_PERIOD = %g # seconds (== %g yr)

############################################################
# Derived properties of the Satellite:
############################################################

# RADIUS          = %g
# MASS            = %g
# DENSITY         = %g
# SURFACE_GRAVITY = %g
# ORBIT_PERIOD    = %g
# MEAN_MOTION     = %g
""" % (self.system_id,\
       self.planet_mass,\
       self.orbit_eccentricity,\
       self.orbit_semimajor_axis,\
       self.nsr_period,\
       self.nsr_period/(31556926),\
       self.radius(),\
       self.mass(),\
       self.density(),\
       self.surface_gravity(),\
       self.orbit_period(),\
       self.mean_motion())

        myStr += """
############################################################
# Layering structure of the Satellite:
############################################################
"""
        n = 0
        while (n < self.num_layers):
            myStr += self.layers[n].ordered_str(n)
            n += 1

        return(myStr)
    # end __str__
# end class Satellite

def avg_body_density(satellite):
    ''' finds average density of satellite body '''
    r = 0
    m = 0
    # for l in satellite.layers[:-1]:  <-- means all execpt the last,
    # so if these are equivalent, then there are 2 items in list?
    # list[0:1] is the same as a copy of list[0]....
    # ok satellite.layers[3] is used, so the two are not equivalent... ck
    for l in satellite.layers[0:1]:
        v = 4.0 / 3.0 * numpy.pi * (numpy.power(r + l.thickness, 3) - numpy.power(r, 3))
        m += l.density * v
        r += l.thickness
    return m / (4.0 / 3.0 * numpy.pi * numpy.power(r, 3))

class SatLayer(object): #
    """
    An object describing a uniform material layer within a satellite.

    Note that a layer by itself has no knowledge of where within the satellite
    it resides.  That information is contained in the ordering of the list of
    layers within the satellite object.

    @ivar layer_id: a string identifying the layer, e.g. C{CORE}, or C{ICE_LOWER}
    @type layer_id: str
    @ivar density: the density of the layer, at zero pressure [kg m^-3],
    @type density: float
    @ivar lame_mu: the layer's Lame parameter, S{mu} (the shear modulus) [Pa].
    @type lame_mu: float
    @ivar lame_lambda: the layer's Lame parameter, S{lambda} [Pa].
    @type lame_lambda: float
    @ivar thickness: the radial thickness of the layer [m].
    @type thickness: float
    @ivar viscosity: the viscosity of the layer [Pa s].
    @type viscosity: float
    """

    def __init__(self, sat, layer_n=0): #
        """
        Construct an object representing a layer within a L{Satellite}.

        Gets values from the L{Satellite.satParams} dictionary for the
        layer that corresponds to the value of C{layer_n}.

        Each layer is defined by seven parameter values, and each layer has a
        unique numeric identifier, appended to the end of all the names of its
        parameters.  Layer zero is the core, with the number increasing as the
        satellite is built up toward the surface.  In the below list the "N" at
        the end of the parameter names should be replaced with the number of
        the layer (C{layer_n}).  Currently, because of the constraints of the
        Love number code that we are using, you must specify 4 layers (C{CORE},
        C{OCEAN}, C{ICE_LOWER}, C{ICE_UPPER}).

          - B{C{LAYER_ID_N}}: A string identifying the layer, e.g. C{OCEAN} or
            C{ICE_LOWER}

          - B{C{DENSITY_N}}: The density of the layer at zero pressure [m
            kg^-3]

          - B{C{LAME_MU_N}}: The real-valued Lame parameter S{mu} (shear
            modulus) [Pa]

          - B{C{LAME_LAMBDA_N}}: The real-valued Lame parameter S{lambda} [Pa]

          - B{C{THICKNESS_N}}: The thickness of the layer [m]

          - B{C{VISCOSITY_N}}: The viscosity of the layer [Pa s]

        Not all layers necessarily require all parameters in order for the
        calculation to succeed, but it is required that they be provided.
        Parameters that will currently be ignored:

          - B{C{VISCOSITY}} of the ocean and the core.

          - B{C{LAME_MU}} of the ocean, assumed to be zero.

        @param sat: the L{Satellite} object to which the layer belongs.
        @type sat: L{Satellite}
        @param layer_n: layer_n indicates which layer in the satellite is being
        defined, with n=0 indicating the center of the satellite, and
        increasing outward.  This is needed in order to select the appropriate
        values from the L{Satellite.satParams} dictionary.
        @type layer_n: int

        @raise MissingSatelliteParamError: if any of the seven input parameters
        listed above is not found in the L{Satellite.satParams} dictionary 
        belonging to the C{sat} object passed in.
        @raise NonNumberSatelliteParamError: if any numeric instance variable
        list above has a value that cannot be converted to a float.
        @raise LowLayerDensityError: if a layer is specififed with a
        density of less than 100 [kg m^-3]
        @raise LowLayerThicknessError: if a layer is specified with a thickness
        of less than 100 [m].
        @raise NegativeLayerParamError: if either of the Lame parameters or the
        viscosity of a layer is less than zero.

        @return: a L{SatLayer} object having the properties specified in the
        @rtype: L{SatLayer}

        """

        # Assign all of the SatLayer's data attributes, converting those
        # that should be numbers into floats, and making sure that if we
        # didn't get a required parameter, an exception is raised.
        try:
            self.layer_id = sat.satParams['LAYER_ID_' + str(layer_n)]
        except KeyError:
            raise MissingSatelliteParamError(sat, 'LAYER_ID_' + str(layer_n))
            
        try:
            self.density = float(sat.satParams['DENSITY_' + str(layer_n)])
        except KeyError:
            raise MissingSatelliteParamError(sat, 'DENSITY_' + str(layer_n))
        except ValueError:
            raise NonNumberSatelliteParamError(sat, 'DENSITY_' + str(layer_n))

        try:
            self.youngs_modulus = float(sat.satParams['YOUNGS_MODULUS_' + str(layer_n)])
        except KeyError:
            raise MissingSatelliteParamError(sat, 'YOUNGS_MODULUS_' + str(layer_n))
        except ValueError:
            raise NonNumberSatelliteParamError(sat, 'YOUNGS_MODULUS_' + str(layer_n))

        try:
            self.poissons_ratio = float(sat.satParams['POISSONS_RATIO_' + str(layer_n)])
        except KeyError:
            raise MissingSatelliteParamError(sat, 'POISSONS_RATIO_' + str(layer_n))
        except ValueError:
            raise NonNumberSatelliteParamError(sat, 'POISSONS_RATIO_' + str(layer_n))

        try:
            self.thickness = float(sat.satParams['THICKNESS_' + str(layer_n)])
        except KeyError:
            raise MissingSatelliteParamError(sat, 'THICKNESS_' + str(layer_n))
        except ValueError:
            raise NonNumberSatelliteParamError(sat, 'THICKNESS_' + str(layer_n))

        try:
            self.viscosity = float(sat.satParams['VISCOSITY_' + str(layer_n)])
        except KeyError:
            raise MissingSatelliteParamError(sat, 'VISCOSITY_' + str(layer_n))
        except ValueError:
            raise NonNumberSatelliteParamError(sat, 'VISCOSITY_' + str(layer_n))

        # Test the SatLayer data attributes to make sure they have nominally
        # reasonable values:

        if self.density <= 100.0:
            raise LowLayerDensityError(sat, layer_n)

        if self.thickness <= 100.0:
            raise LowLayerThicknessError(sat, layer_n)

        if self.youngs_modulus < 0.0:
            raise NegativeLayerParamError(sat, 'YOUNGS_MODULUS_' + str(layer_n))

        if self.poissons_ratio < 0.0:
            raise NegativeLayerParamError(sat, 'POISSONS_RATIO_' + str(layer_n))

        if self.viscosity < 0.0:
            raise NegativeLayerParamError(sat, 'VISCOSITY_' + str(layer_n))

        self.lame_mu = None
        self.lame_lambda = None
        self.bulk_modulus = None
        self.set_bulk_modulus(self.youngs_modulus, self.poissons_ratio)
        self.set_lame_mu(self.youngs_modulus, self.poissons_ratio)
        self.set_lame_lambda(self.bulk_modulus, self.lame_mu)
    # end __init__()

    def __str__(self):
        """
        Output a human and machine readable text description of the layer.
        
        Note that because the layer object does not know explicitly where it is
        within the stratified L{Satellite} (that information is contained in
        the ordering of Satellite.layers list).  

        Derived quantities are also output for clarity, but are preceeded by
        the comment character '#' to avoid accidental use in re-constituting a
        satellite.
        """
        return(self.ordered_str())
    # end __str__()

    def ordered_str(self, layer_n=None):
        """
        Return a string representation of the L{SatLayer}, including information
        about where in the L{Satellite} the layer lies.

        When being output as part of a L{Satellite} object, the L{SatLayer} needs
        to communicate which layer it is.

        @param layer_n: Layer location within satellite.  Zero is the core.  Defaults
        to None, in which case no ordering information is conveyed in output.
        @type layer_n: int
        """
        if layer_n is not None:
            order_str = "_"+str(layer_n)
        else:
            order_str = ""

        myStr = """
LAYER_ID%s    = %s
DENSITY%s     = %g
# LAME_MU%s     = %g
# LAME_LAMBDA%s = %g
THICKNESS%s   = %g
VISCOSITY%s   = %g
# MAXWELL_TIME%s    = %g
# BULK_MODULUS%s    = %g
YOUNGS_MODULUS%s  = %g
POISSONS_RATIO%s  = %g
# P_WAVE_VELOCITY%s = %g
""" % (order_str, str(self.layer_id),\
       order_str, self.density,\
       order_str, self.lame_mu,\
       order_str, self.lame_lambda,\
       order_str, self.thickness,\
       order_str, self.viscosity,\
       order_str, self.maxwell_time(),\
       order_str, self.bulk_modulus,\
       order_str, self.youngs_modulus,\
       order_str, self.poissons_ratio,\
       order_str, self.p_wave_velocity())
        return(myStr)
    # Methods for calculating quantities that depend on data attributes 
    """def bulk_modulus(self):
        #defaults to 2e+9 if Poisson's Ratio is 0.5, a fluid
        if (self.poissons_ratio == 0.5):
            return (2e+9)
        return(self.youngs_modulus / (3.0 * (1.0 - 2.0 * self.poissons_ratio)))"""

    def set_bulk_modulus(self,young,poisson):
        """Calculate the bulk modulus (S{kappa}) of the layer [Pa]."""
        #defaults to 2e+9 if Poisson's Ratio is 0.5, a fluid
        if (poisson == 0.5):
            self.bulk_modulus = 2e+9
        else:
            self.bulk_modulus = (young / (3.0 * (1.0 - 2.0 * poisson)))

    def set_lame_mu(self, young, poisson):
        """Calculates and sets the shear modulus (mu) of the layer given Young's Modulus and Poisson's Ratio."""
        self.lame_mu = float(young / (2.0 * (1.0 + poisson)))

    def set_lame_lambda(self, bulk, shear):
        """Calculates the lame parameter (lambda) of the layer given the bulk and shear moduli."""
        self.lame_lambda = float(bulk - (2.0 / 3.0) * shear)

    def set_youngs_modulus(self, bulk, shear):
        """Calculates the Young's modulus (E) of the layer given the bulk and shear moduli."""
        self.youngs_modulus = float(9.0 * bulk * shear / (3.0 * bulk + shear))

    def set_poissons_ratio(self, bulk, shear):
        """Calculates the Poisson's ratio (nu) of the layer given the bulk and shear moduli."""
        self.poissons_ratio = float((3.0 * bulk - 2.0 * shear) / (2.0 * (3.0 * bulk + shear)))

    def maxwell_time(self):
        """Calculate the Maxwell relaxation time of the layer [s] (viscosity/lame_mu)."""
        if(self.viscosity==0.0 or self.viscosity=="NONE" or self.lame_mu==0.0):
            return(0.0)
        else:
            return(self.viscosity / self.lame_mu)

    #def youngs_modulus(self):
    #    """Calculate the Young's modulus (E) of the layer [Pa]."""
    #    return(self.lame_mu * (3.0*self.lame_lamda() + 2.0*self.lame_mu) / (self.lame_lamda() + self.lame_mu))

    #def poissons_ratio(self):
    #    """Calculate poisson's ratio (S{nu}) of the layer [Pa]."""
    #    return(self.lame_lamda() / (2.0 * (self.lame_lamda() + self.lame_mu) ))

    def p_wave_velocity(self):
        """Calculate the velocity of a compression wave in the layer [m s^-1]"""
        return(scipy.sqrt(self.bulk_modulus/self.density))
    # end derived quantities

# end class SatLayer

class LoveNum(object): #
    """A container class for the complex Love numbers: h2, k2, and l2.

    @ivar h2: the degree 2 complex, frequency dependent Love number h.
    @type h2: complex
    @ivar k2: the degree 2 complex, frequency dependent Love number k.
    @type k2: complex
    @ivar l2: the degree 2 complex, frequency dependent Love number l.
    @type l2: complex
    """

    def __init__(self, h2_real=0, h2_imag=0, k2_real=0, k2_imag=0, l2_real=0, l2_imag=0):
        """Using the real and imaginary parts, create complex values."""
        self.h2 = h2_real + h2_imag*1.0j
        self.k2 = k2_real + k2_imag*1.0j
        self.l2 = l2_real + l2_imag*1.0j

    def __str__(self):
        """Return a human readable string representation of the Love numbers"""
        return("h2 = %s\nk2 = %s\nl2 = %s" % (self.h2, self.k2, self.l2))

    def update_h2(self, (real, imag)):
      self.h2 = real + imag*1.0j

    def update_k2(self, (real, imag)):
      self.k2 = real + imag*1.0j

    def update_l2(self, (real, imag)):
      self.l2 = real + imag*1.0j
#end class LoveNum

class PW_Variables(object): #
  """
  This class contains all of the variables needed to calculate Polar Wander.
  This class was created following the format of LoveNum.
  In the future, all of these variables should be stored elsewhere and just passed into this file.
  But because of how it currently works, by reading in the satellite parameters from a file,
  we didn't have enough time to implement such a change across the whole program.

  -Peter Sinclair (2016)
  """

  def __init__(self, thetaRi=0, phiRi=0, thetaRf=0, phiRf=0, thetaTi=0, phiTi=0, thetaTf=0, phiTf=0, locked=True, despinning=False, InitialSpinPeriod=0, FinalSpinPeriod=0):
    self.thetaRInitial = scipy.pi/2 - numpy.radians(thetaRi) # Converts from latitude to colatitude
    self.phiRInitial = numpy.radians(phiRi)
    self.thetaRFinal = scipy.pi/2 - numpy.radians(thetaRf)
    self.phiRFinal = numpy.radians(phiRf)

    self.locked = locked
    #If the satellite is tidally locked, the rotational axis and tidal axis will always be at a right angle to each other.
    if self.locked:
      self.thetaTInitial = self.thetaRInitial + scipy.pi/2
      self.thetaTFinal = self.thetaRFinal + scipy.pi/2
      self.phiTInitial = self.phiRInitial
      self.phiTFinal = self.phiRFinal
    else:
      self.thetaTInitial = scipy.pi/2 - numpy.radians(thetaTi)
      self.phiTInitial = numpy.radians(phiTi)
      self.thetaTFinal = scipy.pi/2 - numpy.radians(thetaTf)
      self.phiTFinal = numpy.radians(phiTf)

    self.despinning = despinning
    if self.despinning:
      self.InitialSpin = (2*scipy.pi) / (InitialSpinPeriod * 60 * 60)
      self.FinalSpin = (2*scipy.pi) / (FinalSpinPeriod * 60 * 60)

  def __str__(self):
    return("thetaRInitial = %s, phiRInitial = %s\nthetaRFinal = %s, phiRFinal = %s\nthetaTInitial = %s, phiTInitial = %s\nthetaTFinal = %s, phiTFinal = %s, locked = %s, despinning = %s, InitialSpin = %s, FinalSpin = %s" %
      (self.thetaRInitial, self.phiRInitial, self.thetaRFinal, self.phiRFinal, self.thetaTInitial, self.phiTInitial, self.thetaTFinal, self.phiTFinal, self.locked,
        self.despinning, self.InitialSpin, self.FinalSpin))
  
  # Updates the various parameters from the GUI in the same way that LoveNum updates.
  def update_thetaRi(self, thetaRi):
    self.thetaRInitial = scipy.pi/2 - numpy.radians(thetaRi)
    if self.locked:
      self.thetaTInitial = self.thetaRInitial + scipy.pi/2

  def update_phiRi(self, phiRi):
    self.phiRInitial = numpy.radians(phiRi)
    if self.locked:
      self.phiTInitial = self.phiRInitial

  def update_thetaRf(self, thetaRf):
    self.thetaRFinal = scipy.pi/2 - numpy.radians(thetaRf)
    if self.locked:
      self.thetaTFinal = self.thetaRFinal + scipy.pi/2

  def update_phiRf(self, phiRf):
    self.phiRFinal = numpy.radians(phiRf)
    if self.locked:
      self.phiTFinal = self.phiRFinal

  def update_thetaTi(self, thetaTi):
    self.thetaTInitial = scipy.pi/2 - numpy.radians(thetaTi)

  def update_phiTi(self, phiTi):
    self.phiTInitial = numpy.radians(phiTi)

  def update_thetaTf(self, thetaTf):
    self.thetaTFinal = scipy.pi/2 - numpy.radians(thetaTf)

  def update_phiTf(self, phiTf):
    self.phiTFinal = numpy.radians(phiTf)

  def lock_body(self, locked):
    self.locked = locked

  def spin_change(self, despinning):
    self.despinning = despinning

  def update_InitialSpin(self, InitialSpinPeriod):
    self.InitialSpin = 2 * scipy.pi / (InitialSpinPeriod * 60 * 60)

  def update_FinalSpin(self, FinalSpinPeriod):
    self.FinalSpin = 2 * scipy.pi / (FinalSpinPeriod * 60 * 60)
# end class PW_Variables

class StressDef(object): #{{{
    """A base class from which particular tidal stress field objects descend.

    Different tidal forcings are specified as sub-classes of this superclass (one
    for each separate forcing).

    In the expressions of the stress fields, the time M{t} is specified in
    seconds, with zero occuring at periapse, in order to be compatible with the
    future inclusion of stressing mechanisms which may have explicit time
    dependence instead of being a function of the satellite's orbital position
    (e.g. a true polar wander trajectory).

    Location is specified within a polar coordinate system having its origin at the
    satellite's center of mass, using the following variables:

      - co-latitude (S{theta}): The arc separating a point on the surface of
        the satellite from the north pole (0 < S{theta} < S{pi}).

      - longitude (S{phi}): The arc separating the meridian of a point and the
        meridian which passes under the average location of the primary
        (planet) in the sky over the course of an orbit (0 < S{phi} < 2S{pi}).
        B{East is taken as positive.}

    Each subclass must define its own version of the three components of the
    membrane stress tensor, C{Ttt}, C{Tpp}, and C{Tpt} (the north-south,
    east-west, and shear stress components) as methods.

    @cvar satellite: the satellite which the stress is being applied to.
    @type satellite: L{Satellite}
    @cvar omega: the forcing frequency associated with the stress.
    @type omega: float
    @cvar love: the Love numbers which result from the given forcing frequency
    and the specified satellite structure.
    @type love: L{LoveNum}
    @cvar dependson: the set of input parameters (as represented by the strings
    used to define them in input files) upon which the stress field depends.
    @type dependson: set
    """

    # These variables (objects) are used in the superclass and in the subclasses,
    # and so need to be defined here.
    omega = 0.0
    satellite = None
    love = LoveNum(0,0,0,0,0,0)
    dependson = set(['DENSITY',\
                     'YOUNGS_MODULUS',\
                     'POISSONS_RATIO',\
                     'THICKNESS',\
                     'VISCOSITY'])
    # The variables of each subclass 
    # indicates whether to use user_inputted love numbers or calcLoveNew()
    useUser = False
    loveUser = LoveNum()
    UserCoordinates = PW_Variables()

    # Common StressDef Methods: 
    def __str__(self):
        """
        Output information about this stress field, including frequency
        dependent parameters.
        """
        myStr = """
%s_SURFACE_DELTA  = %g
%s_OMEGA          = %g
%s_FORCING_PERIOD = %g
%s_MU_TWIDDLE     = %g + %gj
%s_LAMBDA_TWIDDLE = %g + %gj
%s_LOVE_H2        = %s
%s_LOVE_K2        = %s
%s_LOVE_L2        = %s
""" % (self.__name__, self.Delta(),\
       self.__name__, self.omega,\
       self.__name__, self.forcing_period(),\
       self.__name__, self.mu_twiddle().real, self.mu_twiddle().imag,\
       self.__name__, self.lambda_twiddle().real, self.lambda_twiddle().imag,\
       self.__name__, self.love.h2,\
       self.__name__, self.love.k2,\
       self.__name__, self.love.l2)

        return(myStr)

    def forcing_period(self):
        """
        Calculate the forcing period based on the forcing frequency (omega).
        """
        return(2*scipy.pi/self.omega)

    def calcLove(self): #
        """
        Calculate the Love numbers for the satellite and the given forcing.
        If an infinite forcing period is given, return zero valued Love
        numbers.
        This is a wrapper function, which can be used to call different Love
        number codes in the future.

        @raise InvalidLoveNumberError: if the magnitude of the imaginary part
        of any Love number is larger than its real part, if the real part is
        ever less than zero, or if the real coefficient of the imaginary part
        is ever positive.
        """
        # I am separating the Love number calculation out into a separate
        # method because I anticipate having more than one option here,
        # possibly allowing the user to explicitly set their Love numbers in
        # the input file, and possibly wrapping the Fortran code in Python,
        # or translating James Roberts' code into pure python, so I don't have
        # to worry about this thing compiling everywhere...

        # If the period is infinite (omega = 0) then we aren't going to get any
        # stresses from the forcing at all, and our Love number code will probably
        # blow up anyway, so for the moment we'll just set the Love numbers to be
        # null, and deal with them appropriately in the stress calculation elsewhere.

        if self.omega == 0.0:
            self.calcLoveInfinitePeriod()
        else:
            self.calcLoveWahr4LayerExternal()
        
        # Do a little sanity checking on the Love numbers to make sure that we have
        # realistic values:
        # - Imaginary parts should always have smaller magnitudes than real parts
        # - Real parts should be greater than zero
        # - Imaginary parts should be less than zero
        if (abs(self.love.h2.real) < abs(self.love.h2.imag)) or\
           (abs(self.love.l2.real) < abs(self.love.l2.imag)) or\
           (self.love.h2.real < 0) or\
           (self.love.l2.real < 0) or\
           ((self.love.h2.imag*-1j).real > 0) or\
           ((self.love.l2.imag*-1j).real > 0):
            raise InvalidLoveNumberError(self, self.love)
    # end calcLove()

    def calcLoveNew(self):

      if self.omega == 0.0:
          self.calcLoveInfinitePeriod()
      else:
          param = numpy.zeros((6,4),'complex')
          param[0,0] = self.satellite.layers[0].thickness
          param[0,1] = param[0,0] + self.satellite.layers[1].thickness
          param[0,2] = param[0,1] + self.satellite.layers[2].thickness
          param[0,3] = param[0,2] + self.satellite.layers[3].thickness
          param[0,:] = param[0,:]/1000.0

          param[1,:] = numpy.array([1, 0, 1, 1])
          param[2,:] = numpy.array([self.satellite.layers[0].density, self.satellite.layers[1].density, self.satellite.layers[2].density, self.satellite.layers[3].density])
          #Treat the first two layers as elastic.
          param[3,:] = numpy.array([self.satellite.layers[0].lame_mu, self.satellite.layers[1].lame_mu, self.mu_twiddle(2), self.mu_twiddle(3)])
          #Specify liquid layer lambda as 2e9 because liquid inputs do not 
          #always have a Poisson's ratio of 0.5 on input, so bulk modulus 
          #is not computed correctly.
          param[4,:] = numpy.array([self.satellite.layers[0].lame_lambda, 2e+9, self.lambda_twiddle(2), self.lambda_twiddle(3)])
          
          [h2Re,h2Im,k2Re,k2Im,l2Re,l2Im] = love3.thirdLove(10,param,self.omega)
          self.love = LoveNum(h2Re,h2Im,k2Re,k2Im,l2Re,l2Im)

    def calcLoveInfinitePeriod(self): # 
        """
        Return a set of zero Love numbers constructed statically.
        
        This method is included so we don't have to worry about whether the
        Love number code can deal with being given an infinite period.  All
        stresses will relax to zero with an infinite period (since the shear
        modulus S{mu} goes to zero), so it doesn't really matter what we set
        the Love numbers to here.
        """
        # Yes, I know this actually corresponds to a rigid body.
        self.love = LoveNum(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        
    def calcLoveWahr4LayerExternal(self): #
        """Use John Wahr's Love number code to calculate h, k, and l.
        
        At the moment, the code is fairly limited in the kind of input it can
        take.  The specified satellite must:
          - use a Maxwell rheology
          - have a liquid water ocean underlying the ice shell
          - have a 4-layer structure (ice_upper, ice_lower, ocean, core)
        Eventually the Love number code will be more closely integrated with
        this package, allowing more flexibility in the interior structure of
        the satellite.

        A temporary directory named lovetmp-XXXXXXX (where the X's are a random
        hexadecimal number) is created in the current working directory, within
        which the Love number code is run.  The directory is deleted
        immediately following the calculation.

        @raise LoveExcessiveDeltaError: if L{StressDef.Delta}() > 10^9 for
        either of the ice layers.
        """
        # random number generator for creating random temporary filenames:
        import random
        # methods for interfacing with the operating system:
        import os

        # make sure that Delta isn't so big that the Love number code will
        # fail.  For John's Love number code, things work up to about Delta =
        # 10^9 Note that we need to check for both of the ice layers.  Negative
        # array indices count from the END of the list, so, the surface of the
        # satellite is -1, the next-to-surface layer is -2, etc.

        for layer_n in (-1, -2):
            if self.Delta(layer_n) > 1e9:
                raise LoveExcessiveDeltaError(self, layer_n)
       
        orig_dir = os.getcwd()

        # write an input file that John's Love number code can read
        # Because we're doing a bunch of stuff outside of this package's
        # confines here... and we really have no idea what the OS is like
        # we do this inside a try/except structure... and do our best to
        # clean up after ourselves.
        try:
            lovetmp = "lovetmp-%s" % hex(int(random.random()*2**32))[2:-1]
            os.mkdir(lovetmp)

            # We need to be reasonably sure that the Love number code is
            # actually in the current PATH.  If the package has been installed
            # in the normal way, (using setup.py) then we're fine (the Love
            # number code ought to have been installed elsewhere on the
            # system), but if we're in the code repository, or if we're just
            # doing the tests pre-installation, then the code is in a funny
            # spot - so let's add that spot to our PATH explicitly here.
            # Unfortunately, this will only work on Unix machines...
            dirname = os.path.dirname(os.path.abspath(__file__))
            os.environ['PATH'] += ":%s/love/john_wahr" % (dirname,)

            os.chdir(lovetmp)
            love_infile = open("in.love", 'w')

            love_infile.write("""%g		Mean Density of Satellite (g/cm^3)
1		Rheology (1=Maxwell, 0=elastic)
%g		Forcing period (earth days, 86400 seconds)
%g		Viscosity of upper ice layer (Pa sec)
%g		Viscosity of lower ice layer (Pa sec)
%g		Young's modulus for the rocky core
%g		Poisson's ratio for the rocky core
%g		Density of the rocky core (g/cm^3)
%g		Young's modulus for ice
%g		Poisson's ratio for ice
%g		Density of ice
1		Decoupling fluid layer (e.g. global ocean)? (1=yes, 0=no)
%g		Thickness of fluid layer (km)
%g		Density of fluid layer (g/cm^3)
%g		P-wave velocity in fluid layer (km/sec)
%g		Total radius of satellite (km)
%g		Thickness of upper (cold) ice layer (km)
%g		Thickness of lower (warm) ice layer (km)
165		Total number of calculation nodes
3		Number of density discontinuities
154		Node number of innermost boundary
161		Node number of 2nd innermost boundary
163		Node number of 3rd innermost boundary""" % (self.satellite.density()/1000.0,\
                                                    self.forcing_period()/86400,\
                                                    self.satellite.layers[3].viscosity,\
                                                    self.satellite.layers[2].viscosity,\
                                                    self.satellite.layers[0].youngs_modulus,\
                                                    self.satellite.layers[0].poissons_ratio,\
                                                    self.satellite.layers[0].density/1000.0,\
                                                    self.satellite.layers[3].youngs_modulus,\
                                                    self.satellite.layers[3].poissons_ratio,\
                                                    self.satellite.layers[3].density/1000.0,\
                                                    self.satellite.layers[1].thickness/1000.0,\
                                                    self.satellite.layers[1].density/1000.0,\
                                                    self.satellite.layers[1].p_wave_velocity(),\
                                                    self.satellite.radius()/1000.0,\
                                                    self.satellite.layers[3].thickness/1000.0,\
                                                    self.satellite.layers[2].thickness/1000.0) )

            love_infile.close()

            # Actually call the code (it looks at in.love and writes to out.love
            # in our current working directory)
            os.popen('calcLoveWahr4Layer')

            # read in the output file that it creates
            # take the last line and strip() off the newline at the end
            # then split() it by whitespace
            luv_out = open("out.love", 'r')
            luv = list(luv_out)[-1].strip('\n').split()
            luv_out.close()

            # extract those fields that correspond to the love numbers
            (h2_real, h2_imag, k2_real, k2_imag, l2_real, l2_imag) = luv[3:5]+luv[6:8]+luv[9:11]
            # generate a LoveNum Object made of our new love numbers, and return it.
            self.love = LoveNum(float(h2_real), float(h2_imag), float(k2_real), float(k2_imag), float(l2_real), float(l2_imag))

        except:
            # If the above attemped to run the love2 program did cause some
            # kind of error, we should just raise the exception up the chain.
            raise

        finally:
            # Get back to our original working directory, cleanup after the
            # love2 program.  We need to do these things even if our Love number
            # calculation failed somehow...
            os.chdir(orig_dir)
            try:
                os.remove(os.path.join(lovetmp, "in.love"))
                os.remove(os.path.join(lovetmp, "out.love"))
                os.remove(os.path.join(lovetmp, "out.model"))
                os.rmdir(lovetmp)
            except OSError:
                # If the Love number calculation DID fail, then the
                # above attempt to get rid of the files associated with it
                # will probably raise exceptions (directory not found, etc)
                # and that's okay - we just do nothing.
                pass

    # end calcLoveJohnWahr4Layer()

    def Delta(self, layer_n=-1):
        """
        Calculate S{Delta}, a measure of how viscous the layer's response is.

        @param layer_n: indicates which satellite layer Delta should be
        calculated for, defaulting to the surface (recall that layer 0 is the
        core)
        @type layer_n: int
        @return: S{Delta}= S{mu}/(S{omega}*S{eta})
        @rtype: float
        """
        return(self.satellite.layers[layer_n].lame_mu/(self.omega*self.satellite.layers[layer_n].viscosity))

    def Z(self):
        """
        Calculate the value of Z, a constant that sits in front of many terms
        in the potential defined by Wahr et al. (2008).

        @return: Z, a common constant in many of the Wahr et al. potential terms.
        @rtype: float
        """
        numerator = 3.0*pc.G*self.satellite.planet_mass*self.satellite.radius()**2
        denominator = 2.0*self.satellite.orbit_semimajor_axis**3
        return (numerator/denominator)
    
    def mu_twiddle(self, layer_n=-1):
        """
        Calculate the frequency-dependent Lame parameter S{mu} for a Maxwell
        rheology.

        @param layer_n: number of layer for which we want to calculate S{mu}, defaults
        to the surface (with the core being layer zero).
        @return: the frequency-dependent Lame parameter S{mu} for a Maxwell rheology
        @rtype: complex
        """
        return(self.satellite.layers[layer_n].lame_mu*(1.0/(1-1.0j*self.Delta(layer_n))))

    def lambda_twiddle(self, layer_n=-1):
        """
        Calculate the frequency-dependent Lame parameter S{lambda} for a
        Maxwell rheology.

        @param layer_n: number of layer for which we want to calculate S{mu}, defaults
        to the surface (with the core being layer zero).
        @return: the frequency-dependent Lame parameter S{lambda} for a Maxwell
        rheology.
        @rtype: complex
        """
        lame_lambda = self.satellite.layers[layer_n].lame_lambda
        lame_mu = self.satellite.layers[layer_n].lame_mu

        numerator   = (1.0 - 1j*self.Delta(layer_n)*((2.0*lame_mu+3.0*lame_lambda)/(3.0*lame_lambda)))
        denominator = (1.0 - 1j*self.Delta(layer_n))
        return(lame_lambda * (numerator/denominator))

    def alpha(self):
        """
        Calculate the coefficient alpha twiddle for the surface layer (see
        Wahr et al. 2008).

        @return: Calculate the coefficient alpha twiddle for the surface layer
        (see Wahr et al. 2008).
        @rtype: complex
        """
        numerator   = 3.0*self.lambda_twiddle()+2.0*self.mu_twiddle()
        denominator = self.lambda_twiddle()+2.0*self.mu_twiddle()
        return(numerator/denominator)

    def Gamma(self):
        """
        Calculate the coefficient capital Gamma twiddle for the surface layer
        (see Wahr et al. 2008).

        @return: the coefficient capital Gamma twiddle for the surface layer
        (see Wahr et al. 2008).
        @rtype: complex
        """
        return(self.mu_twiddle()*self.love.l2)

    def b1(self):
        """
        Calculate the coefficient beta one twiddle for the surface layer (see Wahr et al. 2008).

        @return: the coefficient beta one twiddle for the surface layer (see Wahr et al. 2008).
        @rtype: complex
        """
        return(self.mu_twiddle()*(self.alpha()*(self.love.h2-3.0*self.love.l2)+3.0*self.love.l2))

    def g1(self):
        """
        Calculate the coefficient gamma one twiddle for the surface layer (see Wahr et al. (2008)).

        @return: the coefficient gamma one twiddle for the surface layer (see Wahr et al. (2008)).
        @rtype: complex
        """
        return(self.mu_twiddle()*(self.alpha()*(self.love.h2-3.0*self.love.l2)-self.love.l2))
    
    def b2(self):
        """
        Calculate the coefficient beta two twiddle for the surface layer (see Wahr et al. (2008)).

        @return: the coefficient beta two twiddle for the surface layer (see Wahr et al. (2008)).
        @rtype: complex
        """
        return(self.mu_twiddle()*(self.alpha()*(self.love.h2-3.0*self.love.l2)-3.0*self.love.l2))

    def g2(self):
        """
        Calculate the coefficient gamma two twiddle for the surface layer (see Wahr et al. (2008)).

        @return: the coefficient gamma two twiddle for the surface layer (see Wahr et al. (2008)).
        @rtype: complex
        """
        return(self.mu_twiddle()*(self.alpha()*(self.love.h2-3.0*self.love.l2)+self.love.l2))

    # end Common StressDef Methods 

    def Ttt(self, theta, phi, t):
        """
        Calculates the S{tau}_S{theta}S{theta} (north-south) component of the
        stress tensor.
        
        In the base class, this is a purely virtual method - it must be defined
        by the subclasses that describe particular tidal stresses.

        @param theta: the co-latitude of the point at which to calculate the stress [rad].
        @type theta: float
        @param phi: the east-positive longitude of the point at which to
        calculate the stress [rad].
        @type phi: float
        @param t: the time, in seconds elapsed since pericenter, at which to perform the
        stress calculation [s].
        @type t: float
        @return: the S{tau}_S{theta}S{theta} component of the 2x2 membrane stress tensor.
        @rtype: float
        """
        pass

    def Tpp(self, theta, phi, t):
        """
        Calculates the S{tau}_S{phi}S{phi} (east-west) component of the stress tensor.
        
        In the base class, this is a purely virtual method - it must be defined
        by the subclasses that describe particular tidal stresses.

        @param theta: the co-latitude of the point at which to calculate the stress [rad].
        @type theta: float
        @param phi: the east-positive longitude of the point at which to
        calculate the stress [rad].
        @type phi: float
        @param t: the time, in seconds elapsed since pericenter, at which to perform the
        stress calculation [s].
        @type t: float
        @return: the S{tau}_S{phi}S{phi} component of the 2x2 membrane stress tensor.
        @rtype: float
        """
        pass

    def Tpt(self, theta, phi, t):
        """
        Calculates the S{tau}_S{phi}S{theta} (off-diagonal) component of the
        stress tensor.
        
        In the base class, this is a purely virtual method - it must be defined
        by the subclasses that describe particular tidal stresses.

        @param theta: the co-latitude of the point at which to calculate the stress [rad].
        @type theta: float
        @param phi: the east-positive longitude of the point at which to
        calculate the stress [rad].
        @type phi: float
        @param t: the time in seconds elapsed since pericenter, at which to perform the
        stress calculation [s].
        @type t: float
        @return: the S{tau}_S{phi}S{theta} component of the 2x2 membrane stress tensor.
        @rtype: float
        """
        pass
#}}} end class StressDef

class NSR(StressDef): #{{{
    """An object defining the stress field which arises from the
    non-synchronous rotation (NSR) of a satellite's icy shell.
    
    NSR is a subclass of L{StressDef}.  See the derivation and detailed
    discussion of this stress field in in Wahr et al. (2008).
    """

    def __init__(self, satellite):
        """Initialize the definition of the stresses due to NSR of the ice shell.
        
        The forcing frequency S{omega} is the frequency with which a point on
        the surface passes through a single hemisphere, because the NSR stress
        field is degree 2 (that is, it's 2x the expected S{omega} from a full
        rotation)

        Because the core is not subject to the NSR forcing (it remains tidally
        locked and synchronously rotating), all stresses within it are presumed
        to relax away, allowing it to deform into a tri-axial ellipsoid, with
        its long axis pointing toward the parent planet.  In order to allow for
        this relaxation the shear modulus (S{mu}) of the core is set to an
        artificially low value for the purpose of the Love number calculation.
        This increases the magnitude of the radial deformation (and the Love
        number h2) significantly.  See Wahr et al. (2008) for complete
        discussion.

        @param satellite: the satellite to which the stress is being applied.
        @type satellite: L{Satellite}
        @return: an object defining the NSR stresses for a particular satellite.
        @rtype: L{NSR}
        """

        self.__name__ = 'NSR'
        self.satellite = satellite
        # Set the core's shear modulus to be low, so we get a core that responds to
        # the NSR forcing nearly as a fluid:
        self.satellite.layers[0].lame_mu = self.satellite.layers[0].lame_mu / 200.0
        self.satellite.layers[0].set_youngs_modulus(self.satellite.layers[0].bulk_modulus, self.satellite.layers[0].lame_mu)
        self.satellite.layers[0].set_poissons_ratio(self.satellite.layers[0].bulk_modulus, self.satellite.layers[0].lame_mu)

        # For NSR, the forcing period is only half the NSR period (or the forcing
        # frequency is twice the NSR frequency) because the shell goes through an
        # an entire oscillation in half a rotation.  That's why this has 4*pi instead
        # of only 2*pi.
        self.omega = 4.0*scipy.pi/self.satellite.nsr_period
        
        if not self.useUser:
          self.calcLoveNew()
        else:
          self.love = self.loveUser

        # Restore the satellite's core to normal...
        self.satellite.layers[0].lame_mu = self.satellite.layers[0].lame_mu*200.0
        self.satellite.layers[0].set_youngs_modulus(self.satellite.layers[0].bulk_modulus, self.satellite.layers[0].lame_mu)
        self.satellite.layers[0].set_poissons_ratio(self.satellite.layers[0].bulk_modulus, self.satellite.layers[0].lame_mu)

        self.dependson |= set(['PLANET_MASS',\
                               'ORBIT_SEMIMAJOR_AXIS',\
                               'NSR_PERIOD'])

    def Ttt(self, theta, phi, t):
        """
        Calculates the S{tau}_S{theta}S{theta} (north-south) component of the
        stress tensor.
        """
        
        if self.omega == 0:
            return(0.0)
        TttN = (self.b1()-self.g1()*scipy.cos(2.0*theta))*scipy.exp(1j*(2.0*phi+self.omega*t))
        TttN = TttN.real
        TttN = TttN * (self.Z()/(2.0*self.satellite.surface_gravity()*self.satellite.radius()))
        return(TttN)

    def Tpp(self, theta, phi, t):
        """
        Calculates the S{tau}_S{phi}S{phi} (east-west) component of the stress tensor.
        """
        if self.omega == 0:
            return(0.0)
        TppN = (self.b2()-self.g2()*scipy.cos(2.0*theta))*scipy.exp(1j*(2.0*phi+self.omega*t))
        TppN = TppN.real
        TppN = TppN * (self.Z()/(2.0*self.satellite.surface_gravity()*self.satellite.radius()))
        return(TppN)

    def Tpt(self, theta, phi, t):
        """
        Calculates the S{tau}_S{phi}S{theta} (off-diagonal) component of the
        stress tensor.
        """
        
        if self.omega == 0:
            return(0.0)
        TptN = self.Gamma()*1j*scipy.exp(1j*(2.0*phi+self.omega*t))*scipy.cos(theta)
        TptN = TptN.real
        TptN = TptN * (2.0*self.Z()/(self.satellite.surface_gravity()*self.satellite.radius()))
        return(TptN)
#}}} end class NSR

class Diurnal(StressDef): #{{{
    """
    An object defining the stress field that arises on a satellite due to an
    eccentric orbit.

    Diurnal is a subclass of L{StressDef}.  See the derivation and detailed
    discussion of this stress field in in Wahr et al. (2008).
    """

    def __init__(self, satellite):
        """
        Sets the object's satellite and omega attributes; calculates Love numbers.

        @param satellite: the satellite to which the stress is being applied.
        @type satellite: L{Satellite}
        @return: an object defining the NSR stresses for a particular satellite.
        @rtype: L{NSR}
        """

        self.__name__ = 'Diurnal'
        self.satellite = satellite
        self.omega = self.satellite.mean_motion()

        if not self.useUser:
          self.calcLoveNew()
        else:
          self.love = self.loveUser

        self.dependson |= set(['ORBIT_ECCENTRICITY',\
                               'ORBIT_SEMIMAJOR_AXIS',\
                               'PLANET_MASS'])

    def Ttt(self, theta, phi, t):
        """
        Calculates the S{tau}_S{theta}S{theta} (north-south) component of the
        stress tensor.
        """
        Ttt1 = 3.0*(self.b1()-self.g1()*scipy.cos(2.0*theta))*scipy.exp(1j*self.omega*t)*scipy.cos(2.0*phi)
        Ttt2 = -1.0*(self.b1()+3.0*self.g1()*scipy.cos(2.0*theta))*scipy.exp(1j*self.omega*t)
        Ttt3 = -4.0*(self.b1()-self.g1()*scipy.cos(2.0*theta))*1j*scipy.exp(1j*self.omega*t)*scipy.sin(2.0*phi)
        TttD = Ttt1 + Ttt2 + Ttt3
        TttD = TttD.real*self.satellite.orbit_eccentricity*self.Z()/(2.0*self.satellite.surface_gravity()*self.satellite.radius())
        return(TttD)

    def Tpp(self, theta, phi, t):
        """
        Calculates the S{tau}_S{phi}S{phi} (east-west) component of the stress tensor.
        """
        Tpp1 = 3.0*(self.b2()-self.g2()*scipy.cos(2.0*theta))*scipy.exp(1j*self.omega*t)*scipy.cos(2.0*phi)
        Tpp2 = -1.0*(self.b2()+3.0*self.g2()*scipy.cos(2.0*theta))*scipy.exp(1j*self.omega*t)
        Tpp3 = -4.0*(self.b2()-self.g2()*scipy.cos(2.0*theta))*1j*scipy.exp(1j*self.omega*t)*scipy.sin(2.0*phi)
        TppD = Tpp1 + Tpp2 + Tpp3
        TppD = TppD.real*self.satellite.orbit_eccentricity*self.Z()/(2.0*self.satellite.surface_gravity()*self.satellite.radius())
        return(TppD)

    def Tpt(self, theta, phi, t):
        """
        Calculates the S{tau}_S{phi}S{theta} (off-diagonal) component of the
        stress tensor.
        """
        Tpt1 = -4.0*self.Gamma()*1j*scipy.exp(1j*(self.omega*t))*scipy.cos(theta)*scipy.cos(2.0*phi)
        Tpt2 = -3.0*self.Gamma()*scipy.exp(1j*self.omega*t)*scipy.cos(theta)*scipy.sin(2.0*phi)
        TptD = Tpt1 + Tpt2
        TptD = TptD.real*2.0*self.satellite.orbit_eccentricity*self.Z()/(self.satellite.surface_gravity()*self.satellite.radius())
        return(TptD)
#}}} end class Diurnal

# This DiurnalObliquity class originated from satstressgui.py but satstress.py was a more appropriate location
class DiurnalObliquity(StressDef): #{{{
    """
    Diurnal tidal stress with obliquity
    """

    # these values may change
    eps = 0 # obliquity
    periapsis_arg = 0 # argument of periapsis - starting orbit position

    def __init__(self, satellite):
        """
        Sets the object's satellite and omega attributes; calculates Love numbers.

        @param satellite: the satellite to which the stress is being applied.
        @type satellite: L{Satellite}
        @return: an object defining the diurnal stresses for a particular satellite.
        @rtype: L{Diurnal}
        """

        self.__name__ = 'DiurnalObliquity'    # attribute class name
        self.satellite = satellite
        self.omega = self.satellite.mean_motion()
        self.calcLoveNew()

        if not self.useUser:
          self.calcLoveNew()
        else:
          self.love = self.loveUser

    def Ttt(self, theta, phi, t):
        """
        Calculates the S{tau}_S{theta}S{theta} (north-south) component of the
        stress tensor.
        """
        
        #Non-obliquity diurnal stresses, total is desired.
        """Ttt1 = 3.0*(self.b1()-self.g1()*scipy.cos(2.0*theta))*scipy.exp(1j*self.omega*t)*scipy.cos(2.0*phi)
        Ttt2 = -1.0*(self.b1()+3.0*self.g1()*scipy.cos(2.0*theta))*scipy.exp(1j*self.omega*t)
        Ttt3 = -4.0*(self.b1()-self.g1()*scipy.cos(2.0*theta))*1j*scipy.exp(1j*self.omega*t)*scipy.sin(2.0*phi)
        TttD = Ttt1 + Ttt2 + Ttt3
        TttD = TttD.real*self.satellite.orbit_eccentricity*self.Z()/(2.0*self.satellite.surface_gravity()*self.satellite.radius())"""

        TttDO = self.g1()*scipy.exp(1j*(self.periapsis_arg + self.omega*t))*scipy.sin(2.0*theta)*scipy.cos(phi)
        TttDO = TttDO.imag*scipy.sin(2.0*self.eps)*self.Z()/(self.satellite.surface_gravity()*self.satellite.radius())

        return(TttDO)

    def Tpp(self, theta, phi, t):
        """
        Calculates the S{tau}_S{phi}S{phi} (east-west) component of the stress tensor.
        """

        #Non-obliquity diurnal stresses, total is desired.
        """Tpp1 = 3.0*(self.b2()-self.g2()*scipy.cos(2.0*theta))*scipy.exp(1j*self.omega*t)*scipy.cos(2.0*phi)
        Tpp2 = -1.0*(self.b2()+3.0*self.g2()*scipy.cos(2.0*theta))*scipy.exp(1j*self.omega*t)
        Tpp3 = -4.0*(self.b2()-self.g2()*scipy.cos(2.0*theta))*1j*scipy.exp(1j*self.omega*t)*scipy.sin(2.0*phi)
        TppD = Tpp1 + Tpp2 + Tpp3
        TppD = TppD.real*self.satellite.orbit_eccentricity*self.Z()/(2.0*self.satellite.surface_gravity()*self.satellite.radius())"""

        TppDO = self.g2()*scipy.exp(1j*(self.periapsis_arg + self.omega*t))*scipy.sin(2.0*theta)*scipy.cos(phi)
        TppDO = TppDO.imag*scipy.sin(2.0*self.eps)*self.Z()/(self.satellite.surface_gravity()*self.satellite.radius())

        return(TppDO)

    def Tpt(self, theta, phi, t):
        """
        Calculates the S{tau}_S{phi}S{theta} (off-diagonal) component of the
        stress tensor.
        """

        #Non-obliquity diurnal stresses, total is desired.
        """Tpt1 = -4.0*self.Gamma()*1j*scipy.exp(1j*(self.omega*t))*scipy.cos(theta)*scipy.cos(2.0*phi)
        Tpt2 = -3.0*self.Gamma()*scipy.exp(1j*self.omega*t)*scipy.cos(theta)*scipy.sin(2.0*phi)
        TptD = Tpt1 + Tpt2
        TptD = TptD.real*2.0*self.satellite.orbit_eccentricity*self.Z()/(self.satellite.surface_gravity()*self.satellite.radius())"""

        TptDO = self.Gamma()*scipy.exp(1j*(self.periapsis_arg + self.omega*t))*scipy.sin(theta)*scipy.sin(phi)
        TptDO = TptDO.imag*scipy.sin(2.0*self.eps)*2.0*self.Z()/(self.satellite.surface_gravity()*self.satellite.radius())

        return(TptDO)
#}}} end class Diurnal with Obliquity


class PolarWander(StressDef): #{{{
    """
    Stresses due to polar wander.
    
    Uses Matsuyama and Nimmo 2007 ("Tectonic patterns on reoriented and despun planetary bodies").

    The stresses appear to be correct, in magnitude and location.  However, the vectors do not
    appear to be pointing in the correct direction.  It looks like they're actually mirrored along the y axis.
    
    Added by Andre Ismailyan and Peter Sinclair, 2016
    """

    def __init__(self, satellite):
        """
        Sets the object's satellite and omega attributes; calculates Love numbers.

        @param satellite: the satellite to which the stress is being applied.
        @type satellite: L{Satellite}
        @return: an object defining the polar wander stresses for a particular satellite.
        @rtype: L{PW}
        """

        self.__name__ = 'PolarWander'    # attribute class name
        self.satellite = satellite

        self.omega = self.satellite.mean_motion()

        self.calcLoveNew()

        if not self.useUser:
          self.calcLoveNew()
        else:
          self.love = self.loveUser
          
        self.PWCoordinates = self.UserCoordinates
        self.LoveNumber = self.love

        if self.PWCoordinates.despinning:
          """
          By changing the rotation rate, the Polar Wander stress can calculate despinning.
          This could be changed to be a separate stress, but it is an elastic model, and combining
          despinning with Polar Wander could be difficult because PW relies on rotation rate as well.
          Currently, the user can calculate despinning individually simply by putting in the same initial
          and final coordinates.
          """
          self.FlatteningInitial = ((self.satellite.radius())**3 * self.PWCoordinates.InitialSpin**2) / (2.0 * pc.G * self.satellite.mass()) * self.love.h2.real
          self.PWCi = self.FlatteningInitial * self.satellite.layers[3].lame_mu * (1.0 + self.satellite.layers[3].poissons_ratio)/(5.0 + self.satellite.layers[3].poissons_ratio)
          self.FlatteningFinal = ((self.satellite.radius())**3 * self.PWCoordinates.FinalSpin**2) / (2.0 * pc.G * self.satellite.mass()) * self.love.h2.real
          self.PWCf = self.FlatteningFinal * self.satellite.layers[3].lame_mu * (1.0 + self.satellite.layers[3].poissons_ratio)/(5.0 + self.satellite.layers[3].poissons_ratio)
        else:
          #For a tidally locked body in synchronous rotation, the rotation rate will be the same as the orbit rate.
          self.Flattening = ((self.satellite.radius())**3 * self.satellite.mean_motion()**2) / (2.0 * pc.G * self.satellite.mass()) * self.love.h2.real
          self.PWCi = self.Flattening * self.satellite.layers[3].lame_mu * (1.0 + self.satellite.layers[3].poissons_ratio)/(5.0 + self.satellite.layers[3].poissons_ratio)
          self.PWCf = self.PWCi
        """
        PWC is a constant which is used in all of the stress calculations.  I just took it out to be calculated separately.
        The model from Matsuyama and Nimmo (2007) is an elastic one.
        The real portion of the h2 Love number should suffice in this case, but this requires checking.
        """
        
    def CosGamma(self, theta, phi, theta_R, phi_R):
        """
        One of the spherical transformations for calculating Polar Wander Stress
        """
        return scipy.cos(theta) * scipy.cos(theta_R) + scipy.sin(theta) * scipy.sin(theta_R) * scipy.cos(phi - phi_R)
        
    def SinGammaCosPsi(self, theta, phi, theta_R, phi_R):
        """
        One of the spherical transformations for calculating Polar Wander Stress
        """
        return scipy.cos(theta_R) * scipy.sin(theta) - scipy.sin(theta_R) * scipy.cos(theta) * scipy.cos(phi - phi_R)
        
    def SinGammaSinPsi(self, theta, phi, theta_R, phi_R):
        """
        One of the spherical transformations for calculating Polar Wander Stress
        """
        return scipy.sin(theta_R) * scipy.sin(phi - phi_R)

    def Ttt(self, theta, phi, t):
        """
        Calculates the S{tau}_S{theta}S{theta} (north-south) component of the
        stress tensor.
        """

        RDi = self.PWCi * ((6 * (self.SinGammaCosPsi(theta, phi, self.PWCoordinates.thetaRInitial, self.PWCoordinates.phiRInitial))**2) + 
          (9 * self.CosGamma(theta, phi, self.PWCoordinates.thetaRInitial, self.PWCoordinates.phiRInitial)**2) - 5) #Calculates the stress from the initial position of the pole

        RDf = self.PWCf * ((6 * (self.SinGammaCosPsi(theta, phi, self.PWCoordinates.thetaRFinal, self.PWCoordinates.phiRFinal))**2) + 
          (9 * self.CosGamma(theta, phi, self.PWCoordinates.thetaRFinal, self.PWCoordinates.phiRFinal)**2) - 5) #Calculates the stress from the reorientation of the pole

        TDi = self.PWCi * ((6 * (self.SinGammaCosPsi(theta, phi, self.PWCoordinates.thetaTInitial, self.PWCoordinates.phiTInitial))**2) + 
          (9 * self.CosGamma(theta, phi, self.PWCoordinates.thetaTInitial, self.PWCoordinates.phiTInitial)**2) - 5) #Calculates the stress from the initial tidal deformation

        TDf = self.PWCf * ((6 * (self.SinGammaCosPsi(theta, phi, self.PWCoordinates.thetaTFinal, self.PWCoordinates.phiTFinal))**2) + 
          (9 * self.CosGamma(theta, phi, self.PWCoordinates.thetaTFinal, self.PWCoordinates.phiTFinal)**2) - 5) #Calculates the stress from the final tidal deformation

        TttPW = (2.0/3.0) * (RDi - RDf - 3.0 * TDi + 3.0 * TDf)

        """
        TttPW = self.PWC * (1.0/2.0) * scipy.sin(theta_R) * 4.0 * scipy.cos(theta_R) * scipy.sin(2.0 * theta) * scipy.cos(phi) - scipy.sin(theta_R) * scipy.cos(2.0 * theta) * (3.0 + scipy.cos(2.0 * phi)) + 10.0 * (scipy.sin(phi))**2
        
        The equation here is in Appendix C of Matsuyama and Nimmo.
        This equation is extremely simplified, and assumes:
        - a reference frame with the z-axis directed along the final north pole
        - x-axis directed along the longitude of the initial north pole
        - amount of reorientation given by theta_R
        - ignore tidal deformation
        It is unused, but has been left here for reference.
        -PS
        """
        return(TttPW)

    def Tpp(self, theta, phi, t):
        """
        Calculates the S{tau}_S{phi}S{phi} (east-west) component of the stress tensor.
        """

        RDi = self.PWCi * ((-6.0 * (self.SinGammaCosPsi(theta, phi, self.PWCoordinates.thetaRInitial, self.PWCoordinates.phiRInitial))**2) + 
          (3.0 * (self.CosGamma(theta, phi, self.PWCoordinates.thetaRInitial, self.PWCoordinates.phiRInitial))**2) + 1.0) #Calculates the stress from the initial position of the pole

        RDf = self.PWCf * ((-6.0 * (self.SinGammaCosPsi(theta, phi, self.PWCoordinates.thetaRFinal, self.PWCoordinates.phiRFinal))**2) + 
          (3.0 * (self.CosGamma(theta, phi, self.PWCoordinates.thetaRFinal, self.PWCoordinates.phiRFinal))**2) + 1.0) #Calculates the stress from the reorientation of the pole

        TDi = self.PWCi * ((-6.0 * (self.SinGammaCosPsi(theta, phi, self.PWCoordinates.thetaTInitial, self.PWCoordinates.phiTInitial))**2) + 
          (3.0 * (self.CosGamma(theta, phi, self.PWCoordinates.thetaTInitial, self.PWCoordinates.phiTInitial))**2) + 1.0) #Calculates the stress from the initial tidal deformation

        TDf = self.PWCf * ((-6.0 * (self.SinGammaCosPsi(theta, phi, self.PWCoordinates.thetaTFinal, self.PWCoordinates.phiTFinal))**2) + 
          (3.0 * (self.CosGamma(theta, phi, self.PWCoordinates.thetaTFinal, self.PWCoordinates.phiTFinal))**2) + 1.0) #Calculates the stress from the final tidal deformation   
        
        TppPW = (2.0/3.0) * (RDi - RDf - 3.0 * TDi + 3.0 * TDf)

        """
        TppPW = self.PWC * (1.0/2.0) * scipy.sin(theta_R) * (scipy.sin(phi))**2 + 12.0 * scipy.cos(theta_R) * scipy.sin(2.0*theta) * scipy.cos(phi) - 3.0 * scipy.sin(theta_R) * scipy.cos(2.0 * theta) * (3 + scipy.cos(2.0 * phi))
        
        The equation here is in Appendix C of Matsuyama and Nimmo.
        This equation is extremely simplified, and assumes:
        - a reference frame with the z-axis directed along the final north pole
        - x-axis directed along the longitude of the initial north pole
        - amount of reorientation given by theta_R
        - ignore tidal deformation
        It is unused, but has been left here for reference.
        -PS
        """
        return(TppPW)

    def Tpt(self, theta, phi, t):
        """
        Calculates the S{tau}_S{phi}S{theta} (off-diagonal) component of the
        stress tensor.
        """

        RDi = self.PWCi * 2.0 * ((self.SinGammaSinPsi(theta, phi, self.PWCoordinates.thetaRInitial, self.PWCoordinates.phiRInitial)) * 
          (self.SinGammaCosPsi(theta, phi, self.PWCoordinates.thetaRInitial, self.PWCoordinates.phiRInitial))) #Calculates the stress from the initial position of the pole

        RDf = self.PWCf * 2.0 * ((self.SinGammaSinPsi(theta, phi, self.PWCoordinates.thetaRFinal, self.PWCoordinates.phiRFinal)) * 
          (self.SinGammaCosPsi(theta, phi, self.PWCoordinates.thetaRFinal, self.PWCoordinates.phiRFinal))) #Calculates the stress from the reorientation of the pole

        TDi = self.PWCi * 2.0 * ((self.SinGammaSinPsi(theta, phi, self.PWCoordinates.thetaTInitial, self.PWCoordinates.phiTInitial)) * 
          (self.SinGammaCosPsi(theta, phi, self.PWCoordinates.thetaTInitial, self.PWCoordinates.phiTInitial))) #Calculates the stress from the initial tidal deformation

        TDf = self.PWCf * 2.0 * ((self.SinGammaSinPsi(theta, phi, self.PWCoordinates.thetaTFinal, self.PWCoordinates.phiTFinal)) * 
          (self.SinGammaCosPsi(theta, phi, self.PWCoordinates.thetaTFinal, self.PWCoordinates.phiTFinal))) #Calculates the stress from the final tidal deformation   
        # The stress here is slightly different from the one given in Matsuyama and Nimmo (2008) because
        # I have used a double angle identity to convert it to a simpler form.
        
        TptPW = (-2.0) * (RDi - RDf - 3.0 * TDi + 3.0 * TDf) 

        """
        TptPW = self.PWC * 4.0 * scipy.sin(theta_R) * scipy.sin(phi) * (scipy.sin(theta_R) * scipy.cos(theta) * scipy.cos(phi) - scipy.cos(theta_R) * scipy.sin(theta))
        
        The equation here is in Appendix C of Matsuyama and Nimmo.
        This equation is extremely simplified, and assumes:
        - a reference frame with the z-axis directed along the final north pole
        - x-axis directed along the longitude of the initial north pole
        - amount of reorientation given by theta_R
        - ignore tidal deformation
        It is unused, but has been left here for reference.
        -PS
        """
        return(TptPW)
#}}} end class Polar Wander


#This IST class originated from satstressgui.py but satstress.py was a more appropriate location.
class IST(StressDef): #{{{
    """
    Ice Shell Thickening
    There is another way to induce stress on the surface this is through ice shell thickening (IST).
    This stress is going to be the same everywhere so the calculation only has to be done once and
    then aggregated over the entire grid.
    """
    """it may change, diffusivity and lam_1 are parameters that become relevant in the viscoelastic
    regime, as descrbed in Nimmo (2004)"""
    delta_tc = 0
    diffusivity = 0
    lam_1 = 0.65 #Stefan parameter. 

    #IST does not have any associated Love numbers.
    useUser = False
    loveUser = None

    def __init__(self, satellite):
        self.satellite = satellite

    def Ttt(self, theta, phi, t):
        """
        Sigma_t = (E*(delta tc)*(delta rho))/((1-nu)*R*(rho w))

        R - satellite radius 

        Young's Modulus
        E = Shear modulus * ((3*lame parameter + 2*shear modulus)/(shear modulus + lame parameter)

        Possion's ratio
        nu = lame parameter / (2*(lame + shear modulus))

        delta rho - density contrast between water and ice

        delta tc - change in thickness (ask)
        """
        delta_rho = self.satellite.layers[1].density - self.satellite.layers[3].density

        if self.diffusivity == 0:
          return numpy.ones(numpy.shape(theta)) * self.satellite.layers[3].youngs_modulus * self.delta_tc * \
              (self.satellite.layers[1].density - self.satellite.layers[3].density) / \
              ((1 - self.satellite.layers[3].poissons_ratio) * self.satellite.radius() * self.satellite.layers[1].density)
        else:
            return numpy.ones(numpy.shape(theta)) * delta_rho / ( self.satellite.layers[3].density * self.satellite.radius()) * \
                scipy.sqrt( 6.0 * self.satellite.layers[3].youngs_modulus * self.satellite.layers[3].viscosity / ( 1.0 - self.satellite.layers[3].poissons_ratio ) ) * \
                2.0*self.lam_1*scipy.sqrt(self.diffusivity) * \
                scipy.special.dawsn( self.delta_tc / ( 2.0 * self.lam_1 * scipy.sqrt(self.diffusivity) ) * scipy.sqrt( self.satellite.layers[3].youngs_modulus / ( 6.0 * self.satellite.layers[3].viscosity * ( 1.0 - self.satellite.layers[3].poissons_ratio ) ) ) )
    
    def Tpp(self, theta, phi, t):
        return self.Ttt(phi, theta, t)

    def Tpt(self, theta, phi, t):
        return numpy.zeros(numpy.shape(theta))
#}}} end class IST

class StressCalc(object): #{{{
    """
    An object which calculates the stresses on the surface of a L{Satellite}
    that result from one or more stress fields.

    @ivar stresses: a list of L{StressDef} objects, corresponding to the stresses which are to
    be included in the calculations done by the L{StressCalc} object.
    @type stresses: list
    """

    def __init__(self, stressdefs):
        """
        Defines the list of stresses which are to be calculated at a given point.

        @param stressdefs: a list of L{StressDef} objects, corresponding to the
        stresses which are to be included in the calculation.
        @type stressdefs: list
        @return:
        @rtype: L{StressCalc}
        """

        self.stresses = stressdefs
    def tensor(self, theta, phi, t): #{{{2
        """
        Calculates surface stresses and returns them as the elements of a
        symmetric 2x2 membrane stress tensor, but in linear form:
        (Ttt,Tpt,Tpp).

        @param theta: the co-latitude of the point at which to calculate the stress [rad].
        @type theta: float
        @param phi: the east-positive longitude of the point at which to
        calculate the stress [rad].
        @type phi: float
        @param t: the time in seconds elapsed since pericenter, at which to perform the
        stress calculation [s].
        @type t: float
        @return: The 3 elements of the symmetric 2x2 membrane stress tensor S{tau}
        @rtype: Numpy.array
        """

        # First calculate all of the stresses, at all locations:
        Ttt = numpy.array( [ stress.Ttt(theta, phi, t) for stress in self.stresses ] )
        Tpp = numpy.array( [ stress.Tpp(theta, phi, t) for stress in self.stresses ] )
        Tpt = numpy.array( [ stress.Tpt(theta, phi, t) for stress in self.stresses ] )

        # Now, if there are multiple stresses (e.g. NSR and Diurnal), sum
        # across them to get the total stress from all fields:
        Ttt = Ttt.sum(axis=0)
        Tpt = Tpt.sum(axis=0)
        Tpp = Tpp.sum(axis=0)

        return(Ttt,Tpt,Tpp)
    # }}}2 end tensor

    def principal_components(self, theta, phi, t): #{{{2
        """
        Calculates the principal components of the surface stresses and returns
        them as a tuple: (tens_mag, tens_az, comp_mag, comp_az), i.e. the
        magnitudes and azimuths (orientation as an angle clockwise from north)
        for the more tensile (more positive) and less tensile (more
        compressive, more negative) stresses.  Azimuths are in radians, always
        between 0 and pi, and are always separated by pi/2 radians.

        !! Note input parameter theta is co-latitude, not latitude
        """
        Ttt, Tpt, Tpp = self.tensor(theta,phi,t)

        eigval_A, theta_A, phi_A, eigval_B, theta_B, phi_B = eigen2(Ttt, Tpt, Tpp)

        az_A = numpy.mod(numpy.arctan2(phi_A,-theta_A), numpy.pi)
        az_B = numpy.mod(numpy.arctan2(phi_B,-theta_B), numpy.pi)

        tens_mag = numpy.where(eigval_A >  eigval_B, eigval_A, eigval_B)
        tens_az  = numpy.where(eigval_A >  eigval_B, az_A, az_B)
        comp_mag = numpy.where(eigval_A <= eigval_B, eigval_A, eigval_B)
        comp_az  = numpy.where(eigval_A <= eigval_B, az_A, az_B)

        return(tens_mag, tens_az, comp_mag, comp_az)

    #}}}2 end principal_components

    def mean_global_stressmag(self, num_samples=10000, time_sec=0.0): #{{{2
        """
        Calculate the stresses on the surface of the satellite at num_samples
        randomly chosen locations, evenly distributed over the sphere, and
        return the mean values of both the more and less tensile stresses.
        """

        rand_thetas, rand_phis = random_loncolatpoints(num_samples)
        pc = self.principal_components(rand_thetas, rand_phis, time_sec)
        return(pc[0].mean(),pc[2].mean())

    #}}}2 end mean_global_stressmag

    def mean_global_stressdiff(self, num_samples=10000, time_sec=0.0): #{{{2
        """
        Calculate the stresses on the surface of the satellite at num_samples
        randomly chosen locations, evenly distributed over the sphere, and
        return the mean values of both the more and less tensile stresses.
        """

        rand_thetas, rand_phis = random_loncolatpoints(num_samples)
        pc = self.principal_components(rand_thetas, rand_phis, time_sec)
        return((pc[0]-pc[2]).mean())

    #}}}2 end mean_global_stressdiff

# end class StressCalc #}}}

def random_loncolatpoints(N): #{{{
    """
    Generate N evenly distributed random points on a sphere.
    """
    return(2*numpy.pi*numpy.random.rand(N), numpy.arccos(2*numpy.random.rand(N)-1))
#}}}

# Exception classes (for error handling): #{{{
# 
# We're being pretty strict here: if we get anything dubious, we just raise the
# exception and allow it to kill the program.

class Error(Exception):
    """Base class for errors within the satstress module."""
    pass

class NameValueFileError(Error): # 
    """Base class for errors related to NAME=VALUE style input files."""

class NameValueFileParseError(NameValueFileError): # 
    """Indicates a poorly formatted NAME=VALUE files."""
    def __init__(self, nvf, line):
        """Stores the file and line that generated the parse error.
        
        The file object (nvf) and the contents of the poorly formed line
        (badline) are stored within the exception, so we can print an error
        message with useful debugging information to the user."""

        self.nvf  = nvf
        self.line = line

    def __str__(self):
        return("""
Poorly formatted NAME=VALUE pair found in the following file:

%s

Each non-comment line must consist of a non-whitespace string, followed by \
an '=', followed by another non-whitespace string.  For example: \
    
MY_PARAMETER = 12345.678e+09
    
Here's the line I got stuck on:

%s""" % (self.nvf.name, self.line))
    #  end NameValueFileParseError

class NameValueFileDuplicateNameError(NameValueFileError): # 
    """Indicates multiple copies of the same name in an input file."""
    def __init__(self, nvf, name):
        """Stores the file and the NAME that was found to be multiply defined.
        
        NAME is the key, which has been found to be multiply defined in
        the input file, nvf."""
        self.nvf  = nvf
        self.name = name

    def __str__(self):
        return("""
The name %s was defined multiple times in the input file:

%s

Each NAME in a NAME=VALUE file must be defined only once, otherwise it's
unclear which VALUE ought to be stored.
""" % (self.name, self.nvf.name))
    #  end NameValueFileDuplicateNameError

#  end NameValueFileError classes

class SatelliteParamError(Error): # 
    """Indicates a problem with the Satellite initialization."""
    pass

class MissingSatelliteParamError(SatelliteParamError):
    """Indicates a required parameter was not found in the input file."""
    def __init__(self, sat, missingname):
        self.sat = sat
        self.missingname = missingname
    def __str__(self):
        return("""
The required parameter %s was not found in the satellite definition file:

%s
""" % (self.missingname, self.sat.sourcefilename))

class InvalidSatelliteParamError(SatelliteParamError):
    """Raised when a required parameter is not found in the input file."""
    def __init__(self, sat):
        """
Default initialization of an InvalidSatelliteParamError

Simply sets self.sat = sat (a Satellite object).  Most errors can be well
described using only the parameters stored in the Satellite object.
"""
        self.sat = sat

class LargeEccentricityError(InvalidSatelliteParamError):
    """Raised when satellite orbital eccentricity is > 0.25"""
    def __str__(self):
        return("""
The specified orbital eccentricity (e = %g) found in the file

%s

is too big.  The model implemented by this module, and described in Wahr et al. \
(2008) requires that orbital eccentricity be relatively small (< 0.25), for \
larger values of eccentricity, some of the mathematical approximations used \
break down.
""" % (self.sat.orbit_eccentricity, self.sat.sourcefilename))

class NegativeNSRPeriodError(InvalidSatelliteParamError):
    """Raised if the satellite's NSR period is less than zero"""
    def __str__(self):
        return("""
The input file:

%s

specifies an NSR_PERIOD < 0.0:

NSR_PERIOD = %g

but NSR_PERIOD can't be negative...
""" % (self.sat.sourcefilename, self.sat.nsr_period))

class ExcessiveSatelliteMassError(InvalidSatelliteParamError):
    """
    Raised if the satellite's parent planet is less than 10x as massive as the
    satellite."""

    def __str__(self):
        return("""
The mass of the planet you specified:

PLANET_MASS = %g kg

is not that much larger than the mass of the satellite:

SATELLITE_MASS = %g kg

That seems unlikely to be correct.  Maybe you have a units problem, or \
accidentally set the density of one of the layers in the satellite to be \
unrealistically high? \
""" % (self.sat.planet_mass, self.sat.mass()))

class LoveLayerNumberError(InvalidSatelliteParamError):
    """
    Raised if the number of layers specified in the satellite definition file
    is incompatible with the Love number code.
    """
    def __str__(self):
        return("""
You specified a satellite with %d layers in the file:

%s

Unfortunately, at the moment the Love number code that the model uses (based on \
Dahlen, 1976) can only deal with a 4 layer satellite (core, ocean, and a \
2-layer ice shell).  You can, however, set both of the ice layers to have \
identical properties, and vary the thicknesses of the layers significantly to \
minimize their effects if you want.
""" % (self.sat.num_layers, self.sat.sourcefilename))

class InvalidLoveNumberError(InvalidSatelliteParamError):
    """Raised if the Love numbers are found to be suspicious."""
    def __init__(self, stress, love):
        self.stress = stress

    def __str__(self):
        return("""
Sanity check failed on %s Love numbers:

%s

Real parts must have larger magnitudes than imaginary parts.  Real parts must \
be positive.  Imaginary parts must be negative.  Satellite specification was \
read in from the following file: \

%s

""" % (stress.__name__, stress.love, stress.satellite.sourcefilename))
class LoveExcessiveDeltaError(InvalidSatelliteParamError):
    """Raised when S{Delta} > 10^9 for any of the ice layers, at which point
    the Love number code becomes unreliable."""
    def __init__(self, stress, layer_n):
        self.stress = stress
        self.n = layer_n

    def __str__(self):
        return("""
The value of Delta (= lame_mu/(viscosity*omega), a measure of how viscously \
a layer responds to an applied forcing frequency) was found to be too large: \

Delta = %g

in layer %d (LAYER_ID = %s), under the %s forcing.

specified in the input file:

%s

For the Love number code to produce reliable results, Delta should be less than \
about 10^9.  Either reduce the forcing period or increase the viscosity.  If \
you want an infinite forcing period (i.e. if you want all stresses due to the \
forcing to relax to zero) use: \

FORCING_PERIOD = infinity
""" % (self.stress.Delta(self.n),\
       self.n,\
       self.stress.satellite.layers[self.n].layer_id,\
       self.stress.__name__,\
       self.stress.satellite.sourcefilename))

class GravitationallyUnstableSatelliteError(InvalidSatelliteParamError):
    """
    Raised if the density of layers is found not to decrease as you move toward the
    surface from the center of the satellite.
    """
    def __init__(self, sat, layer_n):
        """Overrides the base InvalidSatelliteParamError.__init__() function.

We also need to know which layers have a gravitationally unstable arrangement.
"""
        self.sat = sat
        self.n   = layer_n

    def __str__(self):
        return("""
You have specified a satellite which is gravitationally unstable because \
one of the layers is on top of another having lower density:

Density of layer %d (%s) = %g [kg m^-3]
Density of layer %d (%s) = %g [kg m^-3]
    
This will cause the Love number code to blow up.  Maybe you have a typo \
in this input file:

%s

Or maybe you ordered the layers incorrectly.  Layer 0 is the core, and the \
layer number increases with increasing distance from the core (it's like the \
radial dimension in spherical coordinates, not like "depth").
""" % (self.n, self.sat.layers[self.n].layer_id, self.sat.layers[self.n].density,\
       self.n-1, self.sat.layers[self.n-1].layer_id, self.sat.layers[self.n-1].density,\
       self.sat.sourcefilename))

class NonNumberSatelliteParamError(InvalidSatelliteParamError):
    """Indicates that a non-numeric value was found for a numerical parameter."""
    def __init__(self, sat, badname):
        self.sat     = sat
        self.badname = badname
    def __str__(self):
        return("""
Found a non-numeric value for a numeric parameter:

%s = %s

in the input file:

%s
""" % (self.badname, self.sat.satParams[self.badname], self.sat.sourcefilename))

class LowLayerDensityError(InvalidSatelliteParamError):
    """Indicates that a layer has been assigned an unrealistically low density."""
    def __init__(self, sat, layer_n):
        self.sat   = sat
        self.n     = layer_n

    def __str__(self):
        return("""
Found an unrealistically low layer density:

DENSITY_%d = %g

specified in the input file:

%s

Recall that all input parameters are in SI (meters, kilograms, seconds) units,
and so, for example, water has a density of 1000 [kg m^-3], not 1.0 [g cm^-3].
""" % (self.n, float(self.sat.satParams['DENSITY_'+str(self.n)]), self.sat.sourcefilename))

class LowLayerThicknessError(InvalidSatelliteParamError):
    """Indicates that a layer has been given too small of a thickness"""
    def __init__(self, sat, layer_n):
        self.sat   = sat
        self.n     = layer_n

    def __str__(self):
        return("""
Found a layer thickness of less than 100 meters:

THICKNESS_%d = %g

specified in the input file:

%s

Recall that all input parameters are in SI (meters, kilograms, seconds) units.
""" % (self.n, float(self.sat.satParams['THICKENSS_'+str(self.n)]), self.sat.sourcefilename))

class NegativeLayerParamError(InvalidSatelliteParamError):
    """Indicates a layer has been given an unphysical material property."""
    def __init__(self, sat, badparam):
        self.sat      = sat
        self.badparam = badparam

    def __str__(self):
        return("""
Found an unexpectedly negative value for %s:

%s = %g

specified in the input file:

%s
""" % (self.badparam, self.badparam, float(self.sat.satParams[self.badparam]), self.sat.sourcefilename))
#}}}
#GNU Terry Pratchett