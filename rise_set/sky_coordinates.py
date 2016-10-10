#!/usr/bin/env python

'''rise_set/sky_coordinates.py
Validated representations of RA and Dec.

Author: Michelle Becker (mbecker@lcogt.net)
April 2011
'''

# Required for true (non-integer) division
from __future__ import division
from builtins import str


# Internal imports
from rise_set.angle import Angle, InvalidAngleError


class RightAscension(Angle):
    def __init__(self, time=None, degrees=None, radians=None, units='time'):

        angle_degrees = time
        angle_radians = radians
        angle_units   = units

        if degrees != None:
            angle_degrees = degrees
            angle_units   = 'arc'

        Angle.__init__(self, angle_degrees, angle_radians, angle_units)

        # Check input
        self.validate_ra()


    def validate_ra(self):
        '''Check the right ascension is valid (0 <= ra < 24 00 00).'''

        # RA cannot be greater than one circle in size
        if ( self.degrees >= 360 or self.degrees < 0 ) :
            msg  = "'%s' is not within allowable RA range" % str(self.degrees)
            msg += " (0 <= ra < 360 degrees)."
            raise InvalidAngleError(msg)

        return True


    def __repr__(self):
        return "RightAscension(%s degrees)" % self.in_degrees()



class Declination(Angle):
    def __init__(self, degrees=None, radians=None, units='arc'):
        Angle.__init__(self, degrees, radians, units)

        # Check input
        self.validate_dec()


    def validate_dec(self):
        '''Check the declination is valid (-90 <= dec <= +90).'''

        # Treat the special case of the north and south poles
        if ( self.degrees > 90   or   self.degrees < -90 ):
            msg  = "'%s' is not within allowable Dec range" % str(self.degrees)
            msg += " (-90 <= dec <= +90 degrees)."
            raise InvalidAngleError(msg)

        return True


    def __repr__(self):
        return "Declination(%s degrees)" % self.in_degrees()
