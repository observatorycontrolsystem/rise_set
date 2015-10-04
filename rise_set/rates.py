#!/usr/bin/env python

'''
rates.py - convenience objects for rates (e.g. proper motion).

Proper motion is an angle traversed over a period of time. This class provides a
simple way to define such a quantity, and get useful transformations. The default
time period is a year, in line with the most common astronomical convention.

Authors: Eric Saunders
         Michelle Becker

March 2011
'''

# Required for true (non-integer) division
from __future__ import division
from builtins import object


class ProperMotion(object):

    def __init__(self, component, time = 'year'):
        '''ProperMotion takes an angular component, which can be either
        an Angle, RightAscension, or Declination object, and an optional
        unit of time, either 'year' or 'century'.'''

        # Check if the units are valid
        possible_units_of_time = ['year', 'century']

        if time not in possible_units_of_time:
            msg = ("Valid units of time are either 'year' or 'century'")
            raise RatesConfigError(msg)

        self.time     = time
        self.component = component


    def in_radians_per_year(self):
        # Get the angle in radians
        radians = self.component.in_radians()

        # Convert the timescale if necessary
        if self.time == 'century':
            radians /= 100

        return radians


    def in_degrees_per_year(self):
        # Get the angle in degrees
        degrees = self.component.in_degrees()

        # Convert the timescale if necessary
        if self.time == 'century':
            degrees /= 100

        return degrees


    def in_radians_per_century(self):
        # Get the angle in radians
        radians = self.component.in_radians()

        # Convert the timescale if necessary
        if self.time == 'year':
            radians *= 100

        return radians


    def in_degrees_per_century(self):
        # Get the angle in degrees
        degrees = self.component.in_degrees()

        # Convert the timescale if necessary
        if self.time == 'year':
            degrees *= 100

        return degrees


    def __key(self):
        return (self.time, self.component)

    def __eq__(self, other):
        return self.__key() == other.__key()

    def __hash__(self):
        return hash(self.__key())

    def __repr__(self):
        return "%s/%s" % (self.component, self.time)



class RatesConfigError(Exception):
    '''Error for invalid constructor arguments to the ProperMotion class.'''
    pass
