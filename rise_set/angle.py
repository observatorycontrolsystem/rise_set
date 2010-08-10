#!/usr/bin/python

'''rise_set/angle.py - multi-unit Angle representation.

This package provides the Angle() class, which provides a simple way to specify
an angle in degrees, radians, RA or Dec, and to transform between these.

Author: Eric Saunders (esaunders@lcogt.net)

May 2010
'''

# Required for true (non-integer) division
from __future__ import division     

# Standard library imports
from math import degrees, radians
import re
import sys


class Angle(object):
    '''A multi-unit Angle implementation. 
       Usage: Angle(degrees|radians|ra|dec=[value])
       e.g. 
            a = Angle()
            a = Angle(degrees=45)
            a = Angle(ra=12 20 34.5)
    '''

    def __init__(self, **kwargs):
        self.unit_mapping = {
                                'degrees'     : self.from_degrees,
                                'radians'     : self.from_radians,
                                'ra'          : self.from_sexegesimal,
                                'dec'         : self.from_sexegesimal,
                             }

        self.degrees = None
        try:
            for key in kwargs:
                if ( key == 'ra' ):
                    self.unit_mapping[key](kwargs[key], ra=True)
                elif ( key == 'dec' ):
                    self.unit_mapping[key](kwargs[key], dec=True)
                else:
                    self.unit_mapping[key](kwargs[key])
        except KeyError as e:
            msg = ("Error constructing Angle: " + str(e)
                  + " is not a valid unit. Try 'degrees', 'radians', 'ra',"
                  + "or 'dec' instead.")

            raise AngleConfigError(msg)
    


    def from_degrees(self, degrees):
        '''Set the Angle using a value provided in degrees.'''
    
        self.degrees = degrees
    
    
    
    def from_radians(self, radians):
        '''Set the Angle using a value provided in radians.'''
        self.degrees = degrees(radians)
    
    
    
    def from_sexegesimal(self, sexegesimal, ra=False, dec=False):
        '''Set the Angle using a provided right ascension or declination.'''
    
        # Check we have been told the type of the sexegesimal number
        if ( not (ra or dec) ):
            raise InvalidAngleError("Sexegesimal number must be qualified"
                                        + " with either 'ra' or 'dec'")
       
        # Match sexegesimal with almost arbitrary delimiters
        match = re.search('^([+-])?(\d\d)[^-\d]*(\d\d)[^-\d]*(\d\d(?:[.]\d+)?)$'
                          ,sexegesimal)

        # Check we extracted three numbers
        if ( match and len(match.groups()) == 4 ):
            sign = match.groups()[0]
            hr, min, sec = map(float, match.groups()[1:4])
        else:
            error = ("Invalid sexagesimal format: Try colon or space delimiters"
                     + " instead (e.g. -12:34:56)")
            raise InvalidAngleError(error)


        # If the input is a right ascension
        if ( ra ):        
            self.validate_ra(sign, hr, min, sec)

            # We have a valid RA, so do the conversion
            self.degrees = self.ra_to_degrees(hr, min, sec)

 
        # If the input is a declination
        if ( dec ):        
            self.validate_dec(hr, min, sec)
 
            # We have a valid declination, so do the conversion
            self.degrees = self.dec_to_degrees(sign, hr, min, sec)


         
    def in_degrees(self):
        '''Return the value of the angle in degrees.'''
        return self.degrees



    def in_radians(self):
        '''Return the value of the angle in radians.'''
        return radians(self.degrees)



    def validate_ra(self, sign, hr, min, sec):
        '''Check the right ascension is valid (0 <= ra < 24 00 00).'''

        # Forbid negative values
        if ( sign == '-' ):
            raise InvalidAngleError("RA cannot be negative")
            
        # Hours must range from 0 to 23 inclusive
        if ( hr > 23 ):
            raise InvalidAngleError("'%d' is not a valid hour in RA" % hr)

        # Minutes must range from 0 to 59 inclusive
        if ( min > 60 ):
            raise InvalidAngleError("'%d' is not a valid minute in RA" % min)

        # Seconds must range from 0 to 59 inclusive
        if ( sec > 60 ):
            raise InvalidAngleError("'%d' is not a valid second in RA" % sec)

        return True



    def ra_to_degrees(self, hr, min, sec):
        '''Convert an RA value (hr/min/sec) to degrees.'''

        return hr*(360/24)  +  min*(360/24/60)  +  sec*(360/24/60/60)



    def validate_dec(self, deg, min, sec):
        '''Check the declination is valid (-90 <= dec <= +90).'''        
    
        # Treat the special case of the north and south poles
        if ( deg == 90   or   deg == -90 ):
            if ( min > 0   or   sec > 0 ):
                raise InvalidAngleError("'%d %d %d' exceeds maximum allowable declination" % (deg, min, sec))

        return True



    def dec_to_degrees(self, sign, deg, min, sec):
        '''Convert a declination value (+-deg/min/sec) to degrees.'''
    
        if ( sign == '-' ):
            deg = -1 * deg
            min = -1 * min
            sec = -1 * sec

        return deg + (min/60) + (sec/60/60)



class InvalidAngleError(Exception):
    '''Error for out-of-range angles provided to the Angle class.'''
    
    def __init__(self, value):
        self.value = value
        
    def __str__(self):
        return self.value



class AngleConfigError(Exception):
    '''Error for invalid constructor arguments to the Angle class.'''

    def __init__(self, value):
        self.value = value
        
    def __str__(self):
        return self.value

