#!/usr/bin/env python

'''
test_utils.py - Unit tests for rise_set.utils

description

Author: Eric Saunders
December 2013
'''

from rise_set.utils import coalesce_adjacent_intervals

import datetime
from nose.tools import assert_equal


class TestUtils(object):

    def setup(self):
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


    def test_coalesce_adjacent_intervals(self):
        input_intervals = list(self.some_adjacent_intervals)
        received = coalesce_adjacent_intervals(self.some_adjacent_intervals)

        assert_equal(received, self.expected_coalescence)

        # Confirm the input list was not modified in place
        assert_equal(input_intervals, self.some_adjacent_intervals)


