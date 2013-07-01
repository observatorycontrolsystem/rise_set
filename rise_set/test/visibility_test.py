#!/usr/bin/env python

from __future__ import division

from nose.tools import eq_, assert_equal, assert_almost_equal, raises, nottest
import datetime

#Import the module to test
from rise_set.visibility import Visibility, coalesce_adjacent_intervals

# Additional support modules
from rise_set.angle import Angle
from rise_set.sky_coordinates import RightAscension, Declination
from rise_set.rates import ProperMotion


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
        self.bpl        = {
                          'latitude'  : Angle(degrees = 34.4332222222),
                          'longitude' : Angle(degrees = -119.863045833)
                          }
        self.dt         = datetime.datetime(year=2011, month=2, day=9)
        self.start_date = datetime.datetime(year=2011, month=2, day=9)
        self.end_date   = datetime.datetime(year=2011, month=2, day=11)
        self.horizon    = 0
        self.twilight   = 'sunrise'


        self.capella    = {
                     'ra'                : RightAscension('05 16 41.36'),
                     'dec'               : Declination('+45 59 52.8'),
                     'ra_proper_motion'  : ProperMotion(RightAscension('00 00 00.07552')),
                     'dec_proper_motion' : ProperMotion(Declination('-00 00 00.42711')),
                     'parallax'          : 0.07729,   # Units: arcsec
                     'rad_vel'           : 30.2,      # Units: km/s (-ve approaches)
                     'epoch'             : 2000,
                   }


        self.visibility = Visibility(self.bpl, self.start_date, self.end_date,
                                     self.horizon, self.twilight)


    def test_coalesce_adjacent_intervals(self):
        received = coalesce_adjacent_intervals(self.some_adjacent_intervals)

        assert_equal(received, self.expected_coalescence)


    def test_can_get_sun_up_intervals(self):

        expected = [
                     (self.dt,
                      self.dt.replace(hour=1,    minute=36,
                                      second=1,  microsecond=743923)),
                     (self.dt.replace(hour=14,   minute=50,
                                      second=42, microsecond=912323),
                      self.dt.replace(day=10))
                    ]
        received = self.visibility.find_when_target_is_up(self.sun, self.dt)

        assert_equal(received, expected)


    def test_can_get_sun_down_intervals(self):

        expected = [
                    (self.dt.replace(hour=1,   minute=36,
                                     second=1, microsecond=743923),
                     self.dt.replace(hour=14,   minute=50,
                                     second=42, microsecond=912323))
                    ]
        received = self.visibility.find_when_target_is_down(self.sun, self.dt)

        assert_equal(received, expected)


    def test_can_get_dark_intervals(self):
        # Same as previous test, but with a simpler interface.

        expected = [
                     (
                       datetime.datetime(2011, 2, 9, 1, 36, 1, 743923),
                       datetime.datetime(2011, 2, 9, 14, 50, 42, 912323)
                     ),
                     (
                       datetime.datetime(2011, 2, 10, 1, 37, 0, 294325),
                       datetime.datetime(2011, 2, 10, 14, 49, 47, 105564)
                     )
                   ]
        received = self.visibility.get_dark_intervals()

        assert_equal(received, expected)


    def test_get_target_up_intervals(self):

        expected = [
                     (
                       datetime.datetime(2011, 2, 9, 0, 0),
                       datetime.datetime(2011, 2, 9, 13, 6, 34, 127859)
                     ),
                     (
                       datetime.datetime(2011, 2, 9, 18, 52, 15, 778256),
                       datetime.datetime(2011, 2, 10, 13, 2, 38, 213278)
                     ),
                     (
                       datetime.datetime(2011, 2, 10, 18, 48, 19, 841582),
                       datetime.datetime(2011, 2, 11, 0, 0)
                     )
                   ]

        received = self.visibility.get_target_intervals(self.capella, up=True)

        assert_equal(received, expected)


    def test_up_intervals_rise_set_transit_within_same_day(self):
        start = datetime.datetime(year=2011, month=10, day=13, hour=4, minute=30)
        end   = datetime.datetime(year=2011, month=10, day=13, hour=10, minute=30)

        observatory = {
                'latitude'  : Angle(degrees=34.4325),
                'longitude' : Angle('-119 51 46'),
                }

        rachels_target = {
                     'ra'                : RightAscension('00 00 34.23'),
                     'dec'               : Declination('-30 45 54.6'),
                     'ra_proper_motion'  : ProperMotion(RightAscension('00 00 00.0')),
                     'dec_proper_motion' : ProperMotion(Declination('-00 00 00.0')),
                     'parallax'          : 0.0,   # Units: arcsec
                     'rad_vel'           : 0.0,      # Units: km/s (-ve approaches)
                     'epoch'             : 2000,
                   }


        v = Visibility(observatory, start, end)

        received = v.get_target_intervals(rachels_target)
        expected = [
                     (
                       datetime.datetime(2011, 10, 13, 2, 8, 2, 406269),
                       datetime.datetime(2011, 10, 13, 11, 1, 25, 259816),
                     )
                    ]

        assert_equal(received, expected)
