# Class definitions and related functions for generating Cycloids in SatStressGUI

# possibly import pytest in order to use assertions also as error checkers
# brought up because asserts issues automatically come up as an error window
import satstress
import random
from numpy import *
from pylab import *
from mpl_toolkits import basemap
# for object serialization : http://pythontips.com/2013/08/02/what-is-pickle-in-python/
import pickle
from matplotlib import pyplot as plt
# Open Source Geospacial libraries
from osgeo import ogr
from osgeo import osr
import os

from lineament import spherical_reckon

CycloidCoords = [] # A global array variable used to store the cycloid coordinates for saving to a shapefile. -PS 2016

class Cycloid(object):
    """
    A one dimensional feature on the surface of a spherical satellite similar to Lineament
    """
    def __init__(self, stresscalc, cycloid_name, threshold, propagation_strength, propagation_speed,
                       start_lat, start_lon, start_dir='West', varyVelocity=False, k=0, maxorbit=360, degree_step=0.5):
        """
        Generate a cycloid given a {stresscalc} object which defines which stresses to 
        calculate, the yield {threshold} value to initiate the crack, {propagation_strength},
        {propagation_speed}, and starting point (lon, lat), starting direction, 
        and stepping interval.

        There also includes the option of varying the crack propagation velocity as a function
        of stress magnitude. For this option, in addition to giving a propogation speed, the
        user must also give a proportionality constant 'k' such that:
        
        velocity = user_given_propogation_speed + k(tensile_str - propogation_strength).
     
        stresscalc is a L{satstress.StressCalc} object defining the field we
        are going to be comparing the lineament to.
        """
        
        ###### check validity of arguments #####
        assert start_dir is not None
        # check that k is valid if varyVelocity option is set
        if varyVelocity == True:
            assert k != 0
            assert(isinstance(k, (float, int, long)))

        ##### initialize #####
        self.cycloid_name = cycloid_name
        self.stresscalc = stresscalc
        self.tensile_str = float(threshold)*1000 # kPa
        self.propagation_strength = float(propagation_strength)*1000   # kPa
        self.propagation_speed = float(propagation_speed)         # m/s
        self.start_lat = float(start_lat)
        self.start_lon = float(start_lon)
        self.dir = start_dir
        self.varyVelocity = varyVelocity
        self.propconst = float(k)

        # lists of lat, lon points of the cycloid in format [pt1, pt2, None, pt2, pt3, None...]
        # for faster plotting of many line segments EXPLAIN why we do so...
        # ^^^^ No idea why they save things like this.  -PS (2016)
        self.lat_cycloidpts = []
        self.lon_cycloidpts = []
        self.current_orbit_pos = 0
        self.time = 0

        self.initiated = False
        self.not_propagated = False
        # list of tuples (lat, long, principle stress comps) summarizing important info 
        # for export/display purposes 
        self.data = []
    
        # keeps track of how far around the satellite the cycloid has traveled
        self.traveled = 0 

        # set the ang_dist traveled each interval i.e. angular distance traveled by cycloid
        # going at propagation_speed every x degrees satellite travels in orbit
        # for which a point for a cycloid plot is generated.
        self.degree_step = degree_step
        self.time_step = self.set_interval(self.degree_step)
        self.ang_dist = self.calc_ang_dist()
    
        self.initiate_split(maxorbit)
    # BEGIN: helper functions
    def set_interval(self, degree):
        """
        Calculates the self.time interval in seconds for finding angular distance
        """
        return degree / 360 * self.stresscalc.stresses[0].satellite.orbit_period()

    def calc_ang_dist(self):
        """
        Calculating angular distance is by finding 
        linear_distance = propagation_speed * self.time_step,
        and then converting to angular distance by
        angular_distance = linear_distance / satellite.radius
        """
        return self.propagation_speed * self.time_step / self.stresscalc.stresses[0].satellite.radius()

    def vary_ang_dist(self, stressmag):
        """
        Used to caculate the velocity to the next step when varying velocity option is set.
        
        Note: There are two separate functions instead of replacing calc_ang_dist with 
        vary_ang_dist w/ k = 0 because vary_ang_dist needs to be calculated every step
        while calc_ang_dist only need to be done once.
        """
        newspeed = self.propagation_speed + self.propconst * (stressmag - self.tensile_str)
        return newspeed * self.time_step / self.stresscalc.stresses[0].satellite.radius()

    def get_stresses(self, lat, lon):
        """
        Retrieves stress components at given point lat, lon (which are in radians?) 
        s1, s1az = tensile stress magnitude, tensile stress direction as azimuth
        s3, s3az = compressive stress magnitude, compressive stress direction as azimuth
        """
        return self.stresscalc.principal_components(pi/2 -lat, lon, self.time)

    def reset_cycl(self):
        """
        Used to reset parameters to begin anew. May also just create a new cycloid
        and save the old ones in some file.
        """
        self.lat_cycloidpts = []
        self.lon_cycloidpts = []
        self.data = []
    # END helper functionsm

    # BEGIN main functions
    def initiate_split(self, orbitmax):
        """
        Initiates, or fails to initiate a crack at the starting location. {Orbitmax} is included
        for passing into propagate_split(). 
        Since the value of the stress components depends on the {self.time} elapsed since periapsis,
        this value is also included as a parameter.
        """

        # get stress components
        s1, s1az, s3, s3az = self.get_stresses(radians(self.start_lat), radians(self.start_lon))

        # if tensile stress > tensile strength of ice
        if s1 > self.tensile_str:
            # update lists
            self.initiated = True
            self.lat_cycloidpts.append(self.start_lat)
            self.lon_cycloidpts.append(self.start_lon)
            self.pt_data = (self.start_lat, self.start_lon, s1, s1az, s3, s3az)
            self.data.append(self.pt_data)
            self.propagate_split(self.pt_data)
            CycloidCoords.append((self.start_lon, self.start_lat)) # Appends to the shapaefile saving variable. -PS 2016

        else:
            # indicate, possibly through popup window that initial point was
            # did not have enough tensile to crack
            print "Unable to initiate cycloid at (%.2f, %.2f)." % (self.start_lat, self.start_lon)
            
    def propagate_split(self, prevpt_data):
        """
        Propagates cycloid from most recent location ({prevlon}, {prevlat}) given in degrees,
        until either the stress levels get too low or reached the user-defined {orbitmax} parameter.
        Since the value of the stress components depends on the {self.time} elapsed since periapsis,
        this value is also included as a parameter.
        """
        
        while self.traveled <= self.current_orbit_pos:
            # passing in as arg reduces need to recalculate and speeds things up
            # unpack stress component of previous point
            prevlat, prevlon, s1, s1az, s3, s3az = self.pt_data

            ##### check arg validity #####
            assert(len(self.lat_cycloidpts) == len(self.lon_cycloidpts))
            assert(len(self.lat_cycloidpts) >= 1)
            # check that prevlat and prevlon inputs are consistent with cycloidpts list
           
            ##### Find next point #####

            # if varyVelocity option is set
            if self.varyVelocity == True:
                self.ang_dist = self.vary_ang_dist(s1)

            if self.dir == 'East':
                if s1az < pi / 2.0:
                    crack_az = s1az + pi/2.0
                else:
                    crack_az = s1az - pi/2.0
            else:
                if s1az < pi/2.0:
                    crack_az = s1az - pi/2.0
                else:
                    crack_az = s1az + pi/2.0

            newlon, newlat = spherical_reckon(radians(prevlon), radians(prevlat), crack_az, self.ang_dist)

            ##### Find new principal component and continue #####
            # newlon, newlat already in radians, so no conversion needed
            (new_s1, new_s1az, new_s3, new_s3az) = self.get_stresses(newlat, newlon)

            # propagate if sufficient strength
            
            if new_s1 > self.propagation_strength:
                # convert newlat, newlon back into degrees
                dnewlat = degrees(newlat)
                dnewlon = LonInBounds(degrees(newlon))

                # append to lat list
                self.lat_cycloidpts.append(dnewlat)
                self.lat_cycloidpts.append(None)
                self.lat_cycloidpts.append(dnewlat)
                # append to lon list
                self.lon_cycloidpts.append(dnewlon)
                self.lon_cycloidpts.append(None)
                self.lon_cycloidpts.append(dnewlon)
                # update data list
                self.pt_data = (dnewlat, dnewlon, new_s1, new_s1az, new_s3, new_s3az)
                self.data.append(self.pt_data)

                CycloidCoords.append((dnewlat, dnewlon))

                # increment and reassign local (to the function) variables
                self.time += self.time_step
                self.traveled += self.degree_step

            else: # new_s1 <= self.propagation_strength
                # failed to propagate at location {prevlat, prevlon} at {self.time},
                # but may propagate at same location at future {self.time}, so
                # increment self.time and orbit traveled, but not lat, lon
                self.time += self.time_step
                self.traveled += self.degree_step

        ######## when broken out of while loop, i.e. self.traveled > orbitmax
        # clean up lists 
        if (len(self.lat_cycloidpts) != 1 and len(self.lon_cycloidpts) != 1):
            # remove last element from list so that lists end with None
            self.lat_cycloidpts.pop(-1)
            self.lon_cycloidpts.pop(-1)
        else:
            # either let user know, or do nothing
            print "Cycloid was initiated at starting point but not propagated"
            self.not_propagated = True

    def plotcoordsonbasemap(self,basemap_ax, figure, orbit_pos, triangles, show_names):
        """
            To see if basemap stress values can be accessed and drawn proportionally
            To see if basemap can be accessed and drawn upon ok
            """
        self.current_orbit_pos = orbit_pos

        if show_names:
            x,y = basemap_ax (self.start_lon, self.start_lat)
            figure.annotate(self.cycloid_name, xy =(x,y) , xycoords='data')
        if self.initiated:
            self.propagate_split(self.pt_data)
        else:
            if triangles:
                basemap_ax.plot(self.start_lon, self.start_lat, 'k^', markersize=10)
        if self.not_propagated and triangles:
            basemap_ax.plot(self.start_lon, self.start_lat, 'w^', markersize=10)
        x,y = basemap_ax(self.lon_cycloidpts[:int(3*orbit_pos/self.degree_step)], self.lat_cycloidpts[:int(3*orbit_pos/self.degree_step)])
        #3 because thats how many points are added per degree_step
        basemap_ax.plot(x, y, 'ko', markersize=0.5, zorder=20)

########### Global functions ########### 

def plotCycloid(cyc, orbit_pos, basemap_ax):
    """
    Plots onto a given basemap the cycloids generated from periapsis until {self.time}.
    In order to be efficient and save self.time from calculations, the cycloid calculation
    will be done in 'real self.time'. Recalculation will be avoided. Specifically, cycloid
    points that are already generated will not be recalculated, only plotted. If the
    cycloid points have been calculated only up to self.time t_o, and {self.time} > t_o, then 
    only the remaining points will be calculated. After all points are know, the cycloid
    will be plotted.
    """
    def plot2basemap():
        #convert coordinates
        x, y = basemap_ax(cyc.lat_cycloidpts, cyc.lon_cycloidpts)
        # plot onto basemap
        basemap_ax.plot(x, y, 'ko', markersize=0.5)

    ### MAIN FCT ###

    assert(len(cyc.lat_cycloidpts) == len(cyc.lon_cycloidpts))

    # find the ideal length of lat, lon lists at specified orbit position
    # this assumes the {orbit_pos} is a multiple of self.degree_step
    # makes cau;gir 
    l = (orbit_pos / self. degree_step + 1) * 3

    if (len(cyc.lat_cycloidpts) == 0) and (orbit_pos != 0):
        # initate and propogate until given orbit position
        cyc.initiate_split(orbit_pos)
        plot2basemap()
    elif len(cyc.lat_cycloidpts) < l:
        # continue propogating until given orbit position
        t = cyc.calc_interval(orbit_pos)
        cyc.propagate_split(cyc.lat_cycloidpts[-1], cyc.lon_cycloidpts[-1], orbit_pos)
        plot2basemap()
    else: # len(cyc.lat_cycloidpts) >= l, so all points are already generated
        plot2basemap()

def LonInBounds(deg):
    """
    Assists spherical_reckon() function
    Function to help make sure the radian->degrees conversion and calculated
    longitude remain in bounds (-180, 180)
    """
    upper = 180
    lower = -180

    # works correctly even if deg is negative
    remainder = deg % 360

    if remainder <= 180:
        return remainder
    else:
        return remainder - 360

def LatInBounds(deg):
    """
    Unclear whether this function is necessary
    Function to help make sure the radian->degrees conversion and calculated
    latitude remain in bounds (-90, 90)
    """
    pass

def plotcoordsonbasemap(stresscalc, basemap_ax,
                        threshold,
                        propagation_strength,
                        propagation_speed,
                        start_lat,
                        start_lon, start_dir, vary_velocity,k,maxorbit):
    """
    To see if basemap stress values can be accessed and drawn proportionally
    To see if basemap can be accessed and drawn upon ok
    """

    cyc = Cycloid(stresscalc, threshold, propagation_strength, propagation_speed, start_lon, start_lat, start_dir, 0.1, vary_velocity, k)
    cyc.initiate_split(maxorbit)
    print 'cycloid made \n Attempting to initiate crack'
    x,y = basemap_ax(cyc.lon_cycloidpts, cyc.lat_cycloidpts)
    basemap_ax.plot(x, y, 'ko', markersize=0.5, zorder=20)

def SaveCycloidAsShape(filename):
    """
    Saves a cycloid as a .shp file without georeferencing, to be used in ArcGIS software.

    Peter Sinclair and Andre Ismailyan (2016)
    """
    print "Saving Cycloid at ", filename
    cycloidline = ogr.Geometry(ogr.wkbLineString)
    for point in CycloidCoords:
        x, y = point
        print point
        cycloidline.AddPoint(x, y)

    outDriver = ogr.GetDriverByName('ESRI Shapefile')
    # Remove output shapefile if it already exists
    if os.path.exists(str(filename)):
        outDriver.DeleteDataSource(str(filename))

    # Create the output shapefile
    outDataSource = outDriver.CreateDataSource(str(filename))
    outLayer = outDataSource.CreateLayer(str(filename), geom_type=ogr.wkbLineString)
    featureDefn = outLayer.GetLayerDefn()
    outFeature = ogr.Feature(featureDefn)
    outFeature.SetGeometry(cycloidline)
    outLayer.CreateFeature(outFeature)

    # Close DataSource
    outDataSource.Destroy()

    print "Done Saving"

#GNU Terry Pratchett
