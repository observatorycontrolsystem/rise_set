#!/usr/bin/python

'''rise_set/declination.py 
Validates angles to store the Dec values in the Angle class.
Michelle Becker for LCOGT
April, 2011'''
#___________________________________________________________
#Imports 
# Required for true (non-integer) division
from __future__ import division
from angle import Angle



# Standard library imports
from math import degrees, radians
import re
import sys
#___________________________________________________________


class Declination(Angle):
    
    def __init__(self, degrees = None, radians = None, units = 'arc', **kwargs):
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
        

