#!/usr/bin/env python

'''
test_utils.py - Unit tests for rise_set.utils

description

Author: Eric Saunders
December 2013
'''
from builtins import object

from rise_set.utils import (coalesce_adjacent_intervals, intersect_intervals,
                            intersect_many_intervals, inverse_intervals)

import datetime


class TestUtils(object):

    def setup_method(self):
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

        assert received == self.expected_coalescence

        # Confirm the input list was not modified in place
        assert input_intervals == self.some_adjacent_intervals


    def test_intersect_case01(self):
        # case 1      |......|
        #     |.....|
        expected = []

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15,  0, 0, 0), datetime.datetime(2013, 6, 15,  3, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert received == expected


    def test_intersect_case02(self):
        # case 2      |......|
        #          |.....|
        expected = [(datetime.datetime(2013,6,15,12,0,0),datetime.datetime(2013,6,15,14,0,0))]

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15,  9, 0, 0), datetime.datetime(2013, 6, 15, 14, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert received == expected


    def test_intersect_case03(self):
        # case 3      |......|
        #              |....|
        expected = [(datetime.datetime(2013,6,15,13,0,0),datetime.datetime(2013,6,15,14,0,0))]

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15, 13, 0, 0), datetime.datetime(2013, 6, 15, 14, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert received == expected


    def test_intersect_case04(self):
        # case 4      |......|
        #                 |.....|
        expected = [(datetime.datetime(2013,6,15,13,0,0),datetime.datetime(2013,6,15,15,0,0))]

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15, 13, 0, 0), datetime.datetime(2013, 6, 15, 16, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert received == expected


    def test_intersect_case05(self):
        # case 5      |......|
        #                      |......|
        expected = []

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15, 16, 0, 0), datetime.datetime(2013, 6, 15, 18, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert received == expected


    def test_intersect_case06(self):
        # case 6      |......|
        #          |............|
        expected = [(datetime.datetime(2013,6,15,12,0,0),datetime.datetime(2013,6,15,15,0,0))]

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15, 11, 0, 0), datetime.datetime(2013, 6, 15, 16, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert received == expected

## Edge cases ##

    def test_intersect_case07(self):
        # case 7      |......|
        #             |......|
        expected = [(datetime.datetime(2013,6,15,12,0,0),datetime.datetime(2013,6,15,15,0,0))]

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert received == expected


    def test_intersect_case08(self):
        # case 8      |......|
        #      |......|
        expected = []

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15,  9, 0, 0), datetime.datetime(2013, 6, 15, 12, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert received == expected


    def test_intersect_case09(self):
        # case 9      |......|
        #                    |......|
        expected = []

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15, 15, 0, 0), datetime.datetime(2013, 6, 15, 18, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert received == expected


    def test_intersect_case10(self):
        # case 10     |......|
        #             |........|
        expected = [(datetime.datetime(2013,6,15,12,0,0),datetime.datetime(2013,6,15,15,0,0))]

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 16, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert received == expected

    def test_intersect_case11(self):
        # case 11     |......|
        #           |........|
        expected = [(datetime.datetime(2013,6,15,12,0,0),datetime.datetime(2013,6,15,15,0,0))]

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15, 11, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        received = intersect_intervals(int1, int2)

        assert received == expected

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

        assert received == expected


    def test_intersect_many(self):
        # |......|
        #    |......|
        #      |......|
        expected = [(datetime.datetime(2013,6,15,14,0,0),datetime.datetime(2013,6,15,15,0,0))]

        int1 = [(datetime.datetime(2013, 6, 15, 12, 0, 0), datetime.datetime(2013, 6, 15, 15, 0, 0))]
        int2 = [(datetime.datetime(2013, 6, 15, 13, 0, 0), datetime.datetime(2013, 6, 15, 16, 0, 0))]
        int3 = [(datetime.datetime(2013, 6, 15, 14, 0, 0), datetime.datetime(2013, 6, 15, 17, 0, 0))]

        received = intersect_many_intervals(int1,int2,int3)
        assert received == expected


    def test_inverse_intervals_outside_bound(self):
        start_bound = datetime.datetime(2013, 6, 15)
        end_bound = datetime.datetime(2013, 6, 16)
        intervals = [(datetime.datetime(2013, 6, 15, 9), datetime.datetime(2013, 6, 15, 11)), (datetime.datetime(2013, 6, 15, 15, 30), datetime.datetime(2013, 6, 15, 17))]
        inverse = inverse_intervals(intervals, start_bound, end_bound)
        expected_intervals = [(start_bound, intervals[0][0]), (intervals[0][1], intervals[1][0]), (intervals[1][1], end_bound)]
        assert inverse == expected_intervals


    def test_inverse_intervals_inside_bound_between_intervals(self):
        start_bound = datetime.datetime(2013, 6, 15)
        end_bound = datetime.datetime(2013, 6, 16)
        intervals = [(datetime.datetime(2013, 6, 14, 9), datetime.datetime(2013, 6, 15, 11)), (datetime.datetime(2013, 6, 15, 15, 30), datetime.datetime(2013, 6, 16, 17))]
        inverse = inverse_intervals(intervals, start_bound, end_bound)
        expected_intervals = [(intervals[0][1], intervals[1][0])]
        assert inverse == expected_intervals


    def test_inverse_intervals_inside_bound_within_intervals(self):
        start_bound = datetime.datetime(2013, 6, 15)
        end_bound = datetime.datetime(2013, 6, 16)
        intervals = [
            (datetime.datetime(2013, 6, 14, 9), datetime.datetime(2013, 6, 14, 11)), (datetime.datetime(2013, 6, 15, 9), datetime.datetime(2013, 6, 15, 11)),
            (datetime.datetime(2013, 6, 15, 15, 30), datetime.datetime(2013, 6, 15, 17)), (datetime.datetime(2013, 6, 16, 15, 30), datetime.datetime(2013, 6, 16, 17))
        ]
        inverse = inverse_intervals(intervals, start_bound, end_bound)
        expected_intervals = [(start_bound, intervals[1][0]), (intervals[1][1], intervals[2][0]), (intervals[2][1], end_bound)]
        assert inverse == expected_intervals
