#!/usr/bin/env python

'''
rates.py - convenience objects for rates (e.g. proper motion).

TODO: description

Author: Eric Saunders
March 2011
'''

# Required for true (non-integer) division
from __future__ import division
from angle import Angle

class ProperMotion(object):

    def __init__(self, component, units = 'year'):

        # Check if the units are valid
        possible_units= ['year', 'century']

        if units not in possible_units:
            msg = ("Insert either year, decade, or century")
            raise RatesConfigError(msg)

        self.units     = units
        self.component = component



    def in_radians_per_year(self):
        # Get the angle in radians
        radians = self.component.in_radians()

        # Convert the timescale if necessary
        if self.units != 'year':
            radians /= 100

        # Return it
        return radians



    def in_degrees_per_year(self):
        # Get the angle in degrees
        degrees = self.component.in_degrees()

        # Convert the timescale if necessary
        if self.units != 'year':
            degrees /= 100

        # Return it
        return degrees



    def in_radians_per_century(self):
        # Get the angle in radians
        radians = self.component.in_radians()

        # Convert the timescale if necessary
        if self.units != 'century':
            radians *= 100

        # Return it
        return radians



    def in_degrees_per_century(self):
        # Get the angle in degrees
        degrees = self.component.in_degrees()

        # Convert the timescale if necessary
        if self.units != 'century':
            degrees *= 100

        # Return it
        return degrees




class RatesConfigError(Exception):
    '''Error for invalid constructor arguments to the Angle class.'''

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value
