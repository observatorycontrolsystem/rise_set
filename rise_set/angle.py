#!/usr/bin/python

'''rise_set/angle.py - multi-unit Angle representation.

This package provides the Angle() class, which provides a simple way to specify
an angle in degrees, radians, RA or Dec, and to transform between these.

Author: Eric Saunders (esaunders@lcogt.net)
Co-Author: Michelle Becker

May 2010
April 2011
'''

# Required for true (non-integer) division
from __future__ import division

# Standard library imports
from math import degrees, radians
import re
import sys


class Angle(object):
    '''A multi-unit Angle implementation.
       Usage: Angle(degrees|radians)
       e.g.
            a = Angle(degrees=45)
            a = Angle(radians=12 20 34.5)
    '''

    def __init__(self, degrees = None, radians = None, units = 'arc', **kwargs):
        self.unit_mapping = {
                                'degrees'     : self.from_degrees,
                                'radians'     : self.from_radians,
                             }
        self.units   = units            
        self.degrees = 0.0
        
        #Check if the units entered were in time or arc
        if self.units not in ['arc', 'time']:
            msg = (self.units + " not a valid unit. Please enter 'time' or 'arc'")
            raise AngleConfigError(msg)
        
        #We expect to only get one item in the dictionary
        # if we get more this means someone messed up, 
        # only take the first on entered and run some validation
        
        if degrees:
            self.measurement = degrees
            key         = 'degrees'
        elif radians:
            self.measurement = radians
            key         = 'radians'
        else:
            msg = ("Please specify an angle in either degrees or radians")
            raise AngleConfigError(msg)
            
        if type(self.measurement) == str:
            self.from_sexegesimal(self.measurement)
        else:
            try:
                self.unit_mapping[key](self.measurement)

            except KeyError as e:
                msg = ("Error constructing Angle: " + str(e)
                      + " is not a valid unit. Try 'degrees', 'radians' instead.")

                raise AngleConfigError(msg)
        

    def from_degrees(self, degrees):
        '''Set the Angle using a value provided in degrees.'''

        self.degrees = degrees


    def from_radians(self, radians):
        '''Set the Angle using a value provided in radians.'''
        self.degrees = degrees(radians)


    def from_sexegesimal(self, sexegesimal):
        '''Convert from sexegesimal into degrees'''

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

        if self.units == 'arc': 
            #Then we know its in seconds of arc
            self.degrees = hr + (min/60) + (sec/60/60)
            
        else: 
            #It must be in time
            self.degrees = hr*(360/24)  +  min*(360/24/60)  +  sec*(360/24/60/60) 
        
        return self.degrees


    def in_degrees(self):
        '''Return the value of the angle in degrees.'''
        return self.degrees


    def in_radians(self):
        '''Return the value of the angle in radians.'''
        return radians(self.degrees)    
    
    
    def in_sexegesimal(self):
        "Convert from degrees to sexegesimal"
        
        if self.units == 'arc':
            total, degHrs  = self.degrees, int(self.degrees)
            
        else:
            total, degHrs = self.degrees * (24/360), int(self.degrees * (24/360))
            
        allMin = (total - degHrs) * 60
        min    = int(allMin)
        sec    = (allMin - min) *60
        
        digits = [degHrs, min, sec]
        
        #Put in format to match the regex valid for from_sexegesimal
        for i, item in enumerate(digits):
            if ('.' in str(item)) and len(str(item).split('.')[0])<2 :
                digits[i] = '0' + str(item)
            if len(str(item)) < 2:
                digits[i] = '0' + str(item)
        
        return "%s %s %s" %(digits[0], digits[1], digits[2]) 
            

        
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

