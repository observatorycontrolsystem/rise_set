#!/usr/bin/env python

'''rise_set/angle.py - multi-unit Angle representation.

This package provides the Angle() class, which provides a simple way to specify
an angle in degrees, radians, RA or Dec, and to transform between these.

Authors: Eric Saunders (esaunders@lcogt.net)
         Michelle Becker (mbecker@lcogt.net)

May 2010
April 2011
'''

# Required for true (non-integer) division
from __future__ import division
from builtins import str
from builtins import map
from builtins import object

# Standard library imports
import math
import re

# Usually one would use the six.string_types for the following,
# but since this is a small project we can use this to remove
# the dependency on six.
try:
    basestring
except NameError:
    basestring = str


class Angle(object):
    '''A multi-unit Angle implementation.
       Usage: Angle(degrees|radians)
       e.g.
            a = Angle(degrees=45)
            a = Angle(radians=12 20 34.5)
    '''

    def __init__(self, degrees = None, radians = None, units = 'arc'):
        self.units   = units
        self.degrees = 0.0

        # Check if the units entered were in time or arc
        if self.units not in ['arc', 'time']:
            msg = (self.units + " not a valid unit. Try 'time' or 'arc'.")
            raise AngleConfigError(msg)

        # Complain if none or both are specified
        if not (bool(degrees) ^ bool(radians)) and (degrees != 0) and (radians != 0):
            msg = ("Specify an angle in either degrees or radians.")
            raise AngleConfigError(msg)

        # Must enter either degrees or radians
        if degrees:
            self.from_degrees(degrees)

        elif radians:
            self.units = 'arc'
            self.from_radians(radians)



    def from_degrees(self, degrees):
        '''Set the Angle using a value provided in degrees.'''

        if isinstance(degrees, basestring):
            self.degrees = self.from_sexegesimal(degrees)
        else:
            if self.units == 'time':
                self.degrees = degrees * 360/24
            else:
                self.degrees = degrees



    def from_radians(self, radians):
        '''Set the Angle using a value provided in radians.'''

        if isinstance(radians, basestring):
            radians = self.from_sexegesimal(radians)

        self.degrees = math.degrees(radians)



    def from_sexegesimal(self, sexegesimal):
        '''Convert from sexegesimal into degrees'''

       # Match sexegesimal with almost arbitrary delimiters
        match = re.search(r'^([+-])?(\d*)[^-\d]+(\d\d?)[^-\d]+(\d\d?(?:[.]\d+)?)$'
                          ,sexegesimal)

        # Check we extracted three numbers
        if ( match and len(match.groups()) == 4 ):
            sign = match.groups()[0]
            hrs, mins, secs = list(map(float, match.groups()[1:4]))

        else:
            error  = "Invalid sexegesimal format '%s': " % sexegesimal
            error += "Try colon or space delimiters instead (e.g. -12:34:56)"
            raise InvalidAngleError(error)

        if self.units == 'arc':
            #Then we know its in seconds of arc
            decimal_value = hrs + (mins/60) + (secs/60/60)

        else:
            #It must be in time
            decimal_value = hrs*(360/24)  +  mins*(360/24/60)  +  secs*(360/24/60/60)

        if sign == '-':
            decimal_value *= -1

        return decimal_value



    def in_degrees(self):
        '''Return the value of the angle in degrees.'''
        return self.degrees



    def in_radians(self):
        '''Return the value of the angle in radians.'''

        return math.radians(self.degrees)


    def in_hours(self):
        '''Return the value of the angle in hours of time.'''

        return self.degrees / 15


    def in_sexegesimal(self, radians = False):
        "Convert from degrees to sexegesimal degrees"

        negative = False
        decimal_value = self.degrees

        if radians:
            self.units = 'arc'
            decimal_value = self.in_radians()

        if decimal_value < 0:
            # Make it positive, and apply negative at the end
            decimal_value *= -1
            negative = True

        if self.units == 'arc':
            total, deg_hrs  = decimal_value, int(decimal_value)

        else:
            # Convert into units of time
            total, deg_hrs = decimal_value * (24/360), int(decimal_value * (24/360))

        # Construct the minutes and seconds from the decimal part
        all_min  = (total - deg_hrs) * 60
        mins     = int(all_min)
        secs     = (all_min - mins) *60

        if negative:
            deg_hrs = "-" + str(deg_hrs)

        return "%s %s %.12g" % (deg_hrs, mins, secs)


    def __repr__(self):
        return "Angle(%s degrees)" % self.in_degrees()

    def __key(self):
        return (self.degrees, self.units)

    def __eq__(self, other):
        return self.__key() == other.__key()

    def __hash__(self):
        return hash(self.__key())



class InvalidAngleError(Exception):
    '''Error for out-of-range angles provided to the Angle class.'''
    pass


class AngleConfigError(Exception):
    '''Error for invalid constructor arguments to the Angle class.'''
    pass
