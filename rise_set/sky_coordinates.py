#!/usr/bin/python

'''rise_set/sky_coordinates.py 
Validates angles to store the RA and Dec values in the Angle class.
Michelle Becker for LCOGT
April, 2011'''
#___________________________________________________________
#Imports 
# Required for true (non-integer) division
from __future__ import division

#The angle class
from angle import Angle, InvalidAngleError, AngleConfigError

# Standard library imports
from math import degrees, radians
import re
import sys
#___________________________________________________________


class RightAscension(Angle):
    
    def __init__(self, sign = '+', degrees = None, radians = None, units = 'time'):
        Angle.__init__(self, degrees, radians, units)
        self.ra = self.in_sexegesimal()
        self.validate_ra(sign, self.ra[:2], self.ra[3:5], self.ra[6:])

         
    def validate_ra(self, sign, hr, min, sec):
        '''Check the right ascension is valid (0 <= ra < 24 00 00).'''

        # Forbid negative values
        if ( sign == '-' ):
            raise InvalidAngleError("RA cannot be negative")

        # Hours must range from 0 to 23 inclusive
        if ( int(hr) > 23 ):
            raise InvalidAngleError("'%d' is not a valid hour in RA" % hr)

        # Minutes must range from 0 to 59 inclusive
        if ( int(min) > 60 ):
            raise InvalidAngleError("'%d' is not a valid minute in RA" % min)

        # Seconds must range from 0 to 59 inclusive
        if ( float(sec) > 60 ):
            raise InvalidAngleError("'%d' is not a valid second in RA" % sec)

        return True
        
class Declination(Angle):
    
    def __init__(self, degrees = None, radians = None, units = 'arc'):
        Angle.__init__(self, degrees, radians, units)
        self.dec = self.in_sexegesimal()
        self.validate_dec(self.dec[:2], self.dec[3:5], self.dec[6:])


    def validate_dec(self, deg, min, sec):
        '''Check the declination is valid (-90 <= dec <= +90).'''

        # Treat the special case of the north and south poles
        if ( deg == 90   or   deg == -90 ):
            if ( min > 0   or   sec > 0 ):
                raise InvalidAngleError("'%d %d %d' exceeds maximum allowable declination" % (deg, min, sec))

        return True
