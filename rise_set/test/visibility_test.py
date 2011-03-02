#!/usr/bin/env python

from __future__ import division

from nose.tools import eq_, assert_equal, assert_almost_equal, raises, nottest
import datetime

#Import the module to test
from rise_set.visibility import Visibility

# Additional support modules
from rise_set.angle import Angle


class TestIntervals(object):

    def setUp(self):

        self.some_adjacent_intervals = [
               # Two contiguous intervals
               (datetime.datetime(year=2011, month=2, day=1, hour=12, minute=0),
                datetime.datetime(year=2011, month=2, day=1, hour=12, minute=30)),
               (datetime.datetime(year=2011, month=2, day=1, hour=12, minute=30),
                datetime.datetime(year=2011, month=2, day=1, hour=13, minute=0)),
               # Three contiguous intervals
               (datetime.datetime(year=2011, month=2, day=1, hour=15, minute=0),
                datetime.datetime(year=2011, month=2, day=1, hour=15, minute=30)),
               (datetime.datetime(year=2011, month=2, day=1, hour=15, minute=30),
                datetime.datetime(year=2011, month=2, day=1, hour=16, minute=0)),
               (datetime.datetime(year=2011, month=2, day=1, hour=16, minute=0),
                datetime.datetime(year=2011, month=2, day=1, hour=16, minute=30)),
               # Two non-contiguous intervals
               (datetime.datetime(year=2011, month=2, day=1, hour=17, minute=0),
                datetime.datetime(year=2011, month=2, day=1, hour=17, minute=30)),
               (datetime.datetime(year=2011, month=2, day=1, hour=18, minute=0),
                datetime.datetime(year=2011, month=2, day=1, hour=18, minute=30)),
              ]


        self.expected_coalescence = [
               # First coalescence
               (datetime.datetime(year=2011, month=2, day=1, hour=12, minute=0),
                datetime.datetime(year=2011, month=2, day=1, hour=13, minute=0)),
               # Second coalescence
               (datetime.datetime(year=2011, month=2, day=1, hour=15, minute=0),
                datetime.datetime(year=2011, month=2, day=1, hour=16, minute=30)),
               # Remaining, non-continous intervals
               (datetime.datetime(year=2011, month=2, day=1, hour=17, minute=0),
                datetime.datetime(year=2011, month=2, day=1, hour=17, minute=30)),
               (datetime.datetime(year=2011, month=2, day=1, hour=18, minute=0),
                datetime.datetime(year=2011, month=2, day=1, hour=18, minute=30)),
            ]

        self.sun        = 'sun'
        self.bpl        = dict(latitude = 34.4332222222, longitude = -119.863045833)
        self.dt         = datetime.datetime(year=2011, month=2, day=9)
        self.start_date = datetime.datetime(year=2011, month=2, day=9)
        self.end_date   = datetime.datetime(year=2011, month=2, day=11)
        self.twilight   = 'sunrise'


        self.capella    = {
                     'ra'                : Angle(),
                     'dec'               : Angle(),
                     'ra_proper_motion'  : Angle(),
                     'dec_proper_motion' : Angle(),
                     'parallax'          : 0.07729,   # Units: arcsec
                     'rad_vel'           : 30.2,      # Units: km/s (-ve approaches)
                     'epoch'             : 2000,
                   }

        self.capella['ra'].from_sexegesimal('05 16 41.36', ra=True)
        self.capella['dec'].from_sexegesimal('+45 59 52.8', dec=True)


        # Aladin units are mas/yr...
        self.capella['ra_proper_motion'].from_sexegesimal('00 00 00.07552', ra=True)
        self.capella['dec_proper_motion'].from_sexegesimal('-00 00 00.42711', dec=True)

        self.visibility = Visibility(self.bpl, self.start_date, self.end_date,
                                     self.twilight)


    def test_coalesce_adjacent_intervals(self):
        received = self.visibility.coalesce_adjacent_intervals(
                                                        self.some_adjacent_intervals)

        assert_equal(received, self.expected_coalescence)


    def test_can_get_sun_up_intervals(self):

        expected = [
                     (self.dt, self.dt.replace(hour=1, minute=36, second=2)),
                     (self.dt.replace(hour=14, minute=50, second=43),
                      self.dt.replace(day=10))
                    ]
        received = self.visibility.find_when_target_is_up(self.sun, self.dt)

        assert_equal(received, expected)


    def test_can_get_sun_down_intervals(self):

        expected = [
                    (self.dt.replace(hour=1, minute=36, second=2),
                     self.dt.replace(hour=14, minute=50, second=43))
                    ]
        received = self.visibility.find_when_target_is_down(self.sun, self.dt)

        assert_equal(received, expected)


    def test_can_get_dark_intervals(self):
        # Same as previous test, but with a simpler interface.

        expected = [
                     (
                       datetime.datetime(2011, 2, 9, 1, 36, 2),
                       datetime.datetime(2011, 2, 9, 14, 50, 43)
                     ),
                     (
                       datetime.datetime(2011, 2, 10, 1, 37),
                       datetime.datetime(2011, 2, 10, 14, 49, 47)
                     )
                   ]
        received = self.visibility.get_dark_intervals()

        assert_equal(received, expected)


    def test_get_target_up_intervals(self):

        expected = [
                     (
                       datetime.datetime(2011, 2, 9, 0, 0),
                       datetime.datetime(2011, 2, 9, 13, 6, 34)
                     ),
                     (
                       datetime.datetime(2011, 2, 9, 18, 52, 16),
                       datetime.datetime(2011, 2, 10, 13, 2, 38)
                     ),
                     (
                       datetime.datetime(2011, 2, 10, 18, 48, 20),
                       datetime.datetime(2011, 2, 11, 0, 0)
                     )
                   ]

        received = self.visibility.get_target_intervals(self.capella, up=True)

        assert_equal(received, expected)


    def test_time_tuple_to_datetime(self):

        time_tuple = (11.0, 45.0, 39.329936620018202)
        expected = datetime.datetime(year=2011, month=2, day=9,
                                     hour=11, minute=45, second=39)

        received = self.visibility.time_tuple_to_datetime(time_tuple, self.dt)

        assert_equal(received, expected)
