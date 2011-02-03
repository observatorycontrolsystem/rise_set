#!/usr/bin/python

from __future__ import division

from astrometry import calc_rise_set


class Telescope(object):

    def __init__(self, latitude, longitude, horizon=0):

        self.latitude  = latitude
        self.longitude = longitude


    def set_date(self, date):

        self.date = date


    def calc_rise_set(self, target):

        return calc_rise_set(target, {'latitude':self.latitude, 'longitude':self.longitude}, self.date)
