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
    
    def __init__(self, degrees = None, radians = None, units = 'time'):
        Angle.__init__(self, degrees, radians, units)

        # Check input
        self.validate_ra()

         
    def validate_ra(self):
        '''Check the right ascension is valid (0 <= ra < 24 00 00).'''
        
        # RA cannot be greater than one circle in size
        if ( self.degrees >= 360 or self.degrees < 0 ) :
            raise InvalidAngleError("'%s' is not within allowable RA range(0 <= ra < 360 degrees)." % (str(self.degrees)))
            
        return True
        
        
        
class Declination(Angle):
    
    def __init__(self, degrees = None, radians = None, units = 'arc'):
        Angle.__init__(self, degrees, radians, units)
        
        # Check input
        self.validate_dec()
        
        
        
    def validate_dec(self):
        '''Check the declination is valid (-90 <= dec <= +90).'''
        
        # Treat the special case of the north and south poles
        if ( self.degrees > 90   or   self.degrees < -90 ):
            raise InvalidAngleError("'%s' is not within allowable Dec range(-90 <= dec <= +90 degrees)." % (str(self.degrees)))

        return True
