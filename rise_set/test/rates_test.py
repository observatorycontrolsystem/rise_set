#!/usr/bin/env python

from __future__ import division

from nose.tools import eq_, assert_equal, assert_almost_equal, raises, nottest

#Import the module to test
from rise_set.rates import ProperMotion


class TestProperMotion(object):
    '''Unit tests for the ProperMotion object.'''

    def setup(self):
        self.pm = ProperMotion(30000, unit='milli-arcseconds/year')


    def teardown(self):
        pass


    def test_something(self):
        assert_equal(self.pm.in_radians_per_year(), 0.0021816615649929119)
