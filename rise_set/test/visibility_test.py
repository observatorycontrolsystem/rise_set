#!/usr/bin/env python

from __future__ import division

from nose.tools import eq_, assert_equal, assert_almost_equal, raises, nottest
import datetime

#Import the module to test
from rise_set.visibility import Visibility



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

        self.target     = 'sun'
        self.bpl        = dict(latitude = 34.4332222222, longitude = -119.863045833)
        self.dt         = datetime.datetime(year=2011, month=2, day=9)
        self.start_date = datetime.datetime(year=2011, month=2, day=9)
        self.end_date   = datetime.datetime(year=2011, month=2, day=10)
        self.twilight   = 'sunrise'

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
        received = self.visibility.find_when_target_is_up(self.target, self.dt)

        assert_equal(received, expected)


    def test_can_get_sun_down_intervals(self):

        expected = [
                    (self.dt.replace(hour=1, minute=36, second=2),
                     self.dt.replace(hour=14, minute=50, second=43))
                    ]
        received = self.visibility.find_when_target_is_down(self.target, self.dt)

        assert_equal(received, expected)
