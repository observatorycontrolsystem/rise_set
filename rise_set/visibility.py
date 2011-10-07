#!/usr/bin/env python

'''
visibility.py - Visibility interval calculations.

TODO: description

Author: Eric Saunders (esaunders@lcogt.net)
February 2011
'''

# Required for true (non-integer) division
from __future__ import division

# Standard libary imports
import datetime

# Internal imports
from astrometry import calc_sunrise_set, calc_rise_set


# Set convenient constants
ONE_DAY = datetime.timedelta(days=1)
MIDNIGHT = datetime.time()

class Visibility(object):

    def __init__(self, site, start_date, end_date, horizon=0, twilight='sunrise'):
        self.site       = site
        self.start_date = start_date
        self.end_date   = end_date
        self.horizon    = horizon
        self.twilight   = twilight

        self.dark_intervals = []


    def get_dark_intervals(self):
        '''Returns a set of datetime 2-tuples, each of which represents an interval
           of uninterrupted darkness. The set of tuples gives the complete dark
           intervals between the Visibility object's start and end date.
        '''

        # Don't compute this again if we've already done it
        if self.dark_intervals:
            return self.dark_intervals

        target = 'sun'

        self.dark_intervals = self.get_target_intervals(target, up=False)

        return self.dark_intervals


    def get_target_intervals(self, target, up=True, horizon=None):
        '''Returns a set of datetime 2-tuples, each of which represents an interval
           of uninterrupted time when the target was below the horizon. The set of
           tuples gives the complete target down intervals between the Visibility
           object's start and end date.
        '''

        if up:
            day_interval_func = self.find_when_target_is_up
        else:
            day_interval_func = self.find_when_target_is_down

        # Find rise/set/transit for each day
        intervals = []
        current_date = self.start_date
        while current_date < self.end_date:
            one_day_intervals = day_interval_func(target, current_date, horizon)

            # Add today's intervals to the accumulating list of intervals
            intervals.extend(one_day_intervals)

            # Move on to tomorrow
            current_date += ONE_DAY

        # Collapse adjacent intervals into continuous larger intervals
        intervals = self.coalesce_adjacent_intervals(intervals)

        return intervals


    def find_when_target_is_down(self, target, dt, horizon=None):
        '''Returns a set of datetime 2-tuples, each of which represents an interval
           of uninterrupted time below the horizon at the specified site, for the
           requested date.

           Note: Even though this function currently ignores times, the dt object
           must be a datetime, *not* a date.
        '''

        # Ensure we only deal with dates, because our rise/set/transit tuple is
        # day specific.
        # TODO: Extend to arbitrary start/end times
        dt = dt.replace(hour=0, minute=0, second=0, microsecond=0)

        # We will calculate down intervals as the inverse of the up intervals
        print "dt:", dt
        up_intervals = self.find_when_target_is_up(target, dt, horizon)
        for i in up_intervals:
            print "up interval:", i

        if not up_intervals:
            print "Got no up intervals!"
            print "dt was:", dt
            print "target was:", target
            exit()

        down_intervals = []

        # If the first value has time 00:00:00, then the target starts up
        if up_intervals[0][0].time() == MIDNIGHT:
            pass

        # Otherwise the target starts down - so there's one extra interval at start
        else:
            down_start = dt
            down_end   = up_intervals[0][0]

            down_intervals.append((down_start, down_end))


        # Proceed through the intervals, extracting the gaps
        for i in range(len(up_intervals) - 1):
            down_start = up_intervals[i][1]
            down_end   = up_intervals[i+1][0]

            down_intervals.append((down_start, down_end))


        # If the target sets before the end of the day, grab that as an
        # extra down interval
        if up_intervals[-1][1].time() != MIDNIGHT:
            down_start = up_intervals[-1][1]
            down_end   = dt + ONE_DAY

            down_intervals.append((down_start, down_end))

        return down_intervals


    def find_when_target_is_up(self, target, dt, horizon=None):
        '''Returns a set of datetime 2-tuples, each of which represents an
           interval of uninterrupted time above the horizon at the specified
           site, for the requested date.

           Note: Even though this function currently ignores times, the dt
           object must be a datetime, *not* a date.
        '''

        # Remove any time component of the provided datetime object
        dt = dt.replace(hour=0, minute=0, second=0, microsecond=0)

        # Get the rise/set/transit times for the target, for this day
        if target == 'sun':
            (transit, rise, set) = calc_sunrise_set(self.site, dt, self.twilight)
        else:
            (transit, rise, set) = calc_rise_set(target, self.site, dt, horizon)

        intervals = []

        print "latitude", self.site['latitude'].in_degrees()
        print "longitude", self.site['longitude'].in_degrees()
        print "twilight", self.twilight
        print "dt", dt
        print "rise:", rise, dt + rise
        print "transit:", transit, dt + transit
        print "set:", set, dt + set

        # Case 1: Overlapping start of day boundary
        # Target rose yesterday, and sets today. Rises again later today.
        #         |       x                                     |
        #         |    x     x                                  |
        #         | x           x                               | x
        #        x|                x                           x|
        #     x   |                   x                     x   |
        #   Rise  0hr  Transit       Set                  Rise 24hr
        if (rise > transit) and (set > transit):

            # Store the first interval - start of day until target set
            absolute_set = set + dt
            intervals.append((dt, absolute_set))

            # Store the second interval - target rise until end of day
            absolute_rise = rise + dt
            intervals.append((absolute_rise, dt + ONE_DAY))


        # Case 2: Rise, set and transit all fall within the day, in order
        # Target rises today, transits, and sets before the day ends
        #         |                      x                     |
        #         |                   x     x                  |
        #         |                x           x               |
        #         |             x                 x            |
        #         |          x                       x         |
        #         |       x                             x      |
        #         0hr   Rise           Transit          Set   24hr
        elif (rise < transit) and (set > transit):
            # Only one interval - rise until target set
            absolute_rise = rise + dt
            absolute_set  = set + dt
            intervals.append((absolute_rise, absolute_set))


        # Case 3: Overlapping end of day boundary
        # Target rose yesterday, and sets today. Rises again later today.
        #                 x    |                                        x    |
        #              x     x |                                     x     x |
        #           x          |x                                 x          |x
        #        x             |   x                           x             |   x
        #     x                |      x                     x                |
        #   Rise       Transit 0hr   Set                  Rise      Transit 24hr
        elif (rise < transit) and (set < transit):
            # Same code as case 1!
            # Store the first interval - start of day until target set
            absolute_set = set + dt
            intervals.append((dt, absolute_set))

            # Store the second interval - target rise until end of day
            absolute_rise = rise + dt
            intervals.append((absolute_rise, dt + ONE_DAY))


        return intervals


    def coalesce_adjacent_intervals(self, intervals):
        '''Combine a set of datetime 2-tuples, coalescing adjacent intervals into
           larger intervals wherever possible.
        '''

        coalesced_intervals = [intervals[0]]
        for interval in intervals[1:]:

            # If the current interval end matches the next interval start...
            if coalesced_intervals[-1][1] == interval[0]:
                # ...the two intervals are contiguous - combine them
                coalesced_intervals[-1] = (coalesced_intervals[-1][0], interval[1])
            else:
                # ...the two intervals are not contiguous - store seperately
                coalesced_intervals.append(interval)

        return coalesced_intervals



    def dark_and_up(self):

        # Get the set of dark intervals
        dark_intervals = self.get_dark_intervals()

        # Get the set of target up intervals
        target_up_intervals = self.get_target_up_intervals()
        # Calculate the intersection between the two interval sets
        # Return the intersection
