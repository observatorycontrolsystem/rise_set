#!/usr/bin/env python

'''
test_utils.py - Unit tests for rise_set.utils

description

Author: Eric Saunders
December 2013
'''
from builtins import object

from rise_set.utils import (coalesce_adjacent_intervals, intersect_intervals,
                            intersect_many_intervals)

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


    def test_intersect_case01(self):
        # case 1      |......|
        #     |.....|
        expected = []

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15,  0, 0, 0), datetime.datetime(2013, 6, 15,  3, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert_equal(received, expected)


    def test_intersect_case02(self):
        # case 2      |......|
        #          |.....|
        expected = [(datetime.datetime(2013,6,15,12,0,0),datetime.datetime(2013,6,15,14,0,0))]

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15,  9, 0, 0), datetime.datetime(2013, 6, 15, 14, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert_equal(received, expected)


    def test_intersect_case03(self):
        # case 3      |......|
        #              |....|
        expected = [(datetime.datetime(2013,6,15,13,0,0),datetime.datetime(2013,6,15,14,0,0))]

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15, 13, 0, 0), datetime.datetime(2013, 6, 15, 14, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert_equal(received, expected)


    def test_intersect_case04(self):
        # case 4      |......|
        #                 |.....|
        expected = [(datetime.datetime(2013,6,15,13,0,0),datetime.datetime(2013,6,15,15,0,0))]

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15, 13, 0, 0), datetime.datetime(2013, 6, 15, 16, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert_equal(received, expected)


    def test_intersect_case05(self):
        # case 5      |......|
        #                      |......|
        expected = []

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15, 16, 0, 0), datetime.datetime(2013, 6, 15, 18, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert_equal(received, expected)


    def test_intersect_case06(self):
        # case 6      |......|
        #          |............|
        expected = [(datetime.datetime(2013,6,15,12,0,0),datetime.datetime(2013,6,15,15,0,0))]

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15, 11, 0, 0), datetime.datetime(2013, 6, 15, 16, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert_equal(received, expected)

## Edge cases ##

    def test_intersect_case07(self):
        # case 7      |......|
        #             |......|
        expected = [(datetime.datetime(2013,6,15,12,0,0),datetime.datetime(2013,6,15,15,0,0))]

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert_equal(received, expected)


    def test_intersect_case08(self):
        # case 8      |......|
        #      |......|
        expected = []

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15,  9, 0, 0), datetime.datetime(2013, 6, 15, 12, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert_equal(received, expected)


    def test_intersect_case09(self):
        # case 9      |......|
        #                    |......|
        expected = []

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15, 15, 0, 0), datetime.datetime(2013, 6, 15, 18, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert_equal(received, expected)


    def test_intersect_case10(self):
        # case 10     |......|
        #             |........|
        expected = [(datetime.datetime(2013,6,15,12,0,0),datetime.datetime(2013,6,15,15,0,0))]

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 16, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert_equal(received, expected)

    def test_intersect_case11(self):
        # case 11     |......|
        #           |........|
        expected = [(datetime.datetime(2013,6,15,12,0,0),datetime.datetime(2013,6,15,15,0,0))]

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15, 11, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert_equal(received, expected)

## Multiple intersections

    def test_intersect_case12(self):
        # case 12     |......|
        #           |...| |....|
        expected = [(datetime.datetime(2013,6,15,12,0,0),datetime.datetime(2013,6,15,13,0,0)),
                    (datetime.datetime(2013,6,15,14,0,0),datetime.datetime(2013,6,15,15,0,0))]

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15, 11, 0, 0), datetime.datetime(2013, 6, 15, 13, 0, 0)),
                (datetime.datetime(2013, 6, 15, 14, 0, 0), datetime.datetime(2013, 6, 15, 16, 0, 0))]

        received = intersect_intervals(int1, int2)

        assert_equal(received, expected)


    def test_intersect_many(self):
        # |......|
        #    |......|
        #      |......|
        expected = [(datetime.datetime(2013,6,15,14,0,0),datetime.datetime(2013,6,15,15,0,0))]

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15, 13, 0, 0), datetime.datetime(2013, 6, 15, 16, 0, 0))]
        int3 = [(datetime.datetime(2013, 6, 15, 14, 0, 0), datetime.datetime(2013, 6, 15, 17, 0, 0))]

        received = intersect_many_intervals(int1,int2,int3)
        assert_equal(received,expected)
