#!/usr/bin/env python

'''
moving_objects_test.py - Tests for moving object support

Author: Eric Saunders
December 2013
'''

from rise_set.angle import Angle
from rise_set.moving_objects import elem_to_topocentric_apparent

from datetime import datetime, timedelta
from nose.tools import assert_equal, assert_almost_equal


class TestMovingObjects(object):

    def setup(self):
        # Asteroid (470) Kilia
        # From MPC webpages: http://minorplanetcenter.net/iau/MPEph/MPEph.html
        # Epoch 2013 Nov. 4.0 TT = JDT 2456600.5                  MPC
        # M  54.47380              (2000.0)            P               Q
        # n   0.26424476     Peri.   47.60658     -0.75565416     +0.65480399             T = 2456394.35096 JDT
        # a   2.4050925      Node   173.25052     -0.63179958     -0.72278517             q =     2.1781735
        # e   0.0943494      Incl.    7.22565     -0.17267333     -0.22093739
        # P   3.73           H   10.07          G   0.15           U   0
        # From 1313 observations at 37 oppositions, 1901-2013, mean residual 0".53.

        self.elements = {
                          'epoch'          : datetime(2013, 11, 4),
                          'inclination'    : Angle(degrees=7.22565),
                          'long_node'      : Angle(degrees=173.25052),
                          'arg_perihelion' : Angle(degrees=47.60658),
                          'semi_axis'      : 2.4050925,
                          'eccentricity'   : 0.0943494,
                          'mean_anomaly'   : Angle(degrees=54.47380),
                        }

        self.elp = {
                     'latitude'  : Angle(degrees=30.6801),
                     'longitude' : Angle(degrees=-104.015194444),
                   }


    def test_elem_to_topocentric_apparent_time1(self):
        dt = datetime(2013, 12, 9)
        tdb = dt - timedelta(seconds=67.184)

        ra, dec = elem_to_topocentric_apparent(tdb, self.elements, self.elp)


        # Coordinates of (470) Kilia from JPL Horizons: 284.70994 -18.19420
        assert_almost_equal(ra.in_degrees(), 284.70994, places=2)
        assert_almost_equal(dec.in_degrees(), -18.19420, places=4)


    def test_elem_to_topocentric_apparent_time2(self):
        dt = datetime(2013, 6, 9)
        tdb = dt - timedelta(seconds=67.184)

        ra, dec = elem_to_topocentric_apparent(tdb, self.elements, self.elp)


        # Coordinates from JPL Horizons: 225.22256  -5.29589
        assert_almost_equal(ra.in_degrees(), 225.22256, places=2)
        assert_almost_equal(dec.in_degrees(), -5.29589, places=2)
