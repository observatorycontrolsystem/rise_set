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


    def __init__(self, value, unit='milli-arcseconds/year'):
        if unit == 'milli-arcseconds/year':
            #TODO: Implement new Angle mode
            self.angle       = Angle(seconds_of_arc=value/1000)
            self.time_period = 'year'


    def in_radians_per_year(self):
        # Get the angle in radians

        # Convert the timescale if necessary

        # Return it

        pass
