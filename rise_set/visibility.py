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
import math

# Internal imports
from rise_set.astrometry import calc_sunrise_set, calc_rise_set, RiseSetError, Star
from rise_set.angle      import Angle

# Import logging modules
import logging

# Set the logger name
_log = logging.getLogger('rise_set.visibility')


# Set convenient constants
ONE_DAY  = datetime.timedelta(days=1)
MIDNIGHT = datetime.time()




def set_airmass_limit(airmass, horizon):
    ''' Compare the provided maximum airmass limit with the horizon of a telescope, and
        return the effective horizon for rise/set purposes (whichever is higher elevation). If
        no airmass is provided, we default to the horizon.

        airmass = 1 / cos(zenith)
        horizon = 90 - zenith'''

    if not airmass:
        return horizon

    # We convert the horizon to airmass, not vice versa, to avoid small angle problems if a
    # very large airmass is provided
    zenith_distance = 90 - horizon
    horizon_airmass = 1 / math.cos(math.radians(zenith_distance))

    effective_horizon = horizon
    if airmass < horizon_airmass:
        effective_horizon = 90 - math.degrees(math.acos(1 / airmass))


    return effective_horizon


def coalesce_adjacent_intervals(intervals):
    '''Combine a set of datetime 2-tuples, coalescing adjacent intervals into
       larger intervals wherever possible.
    '''

    # Catch the special case where the target never rose
    if len(intervals) == 0:
        return intervals

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


class Visibility(object):

    def __init__(self, site, start_date, end_date, horizon=0, twilight='sunrise'):
        self.site       = site
        self.start_date = start_date
        self.end_date   = end_date
        self.horizon    = Angle(degrees=horizon)
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


    def get_target_intervals(self, target, up=True, airmass=None):
        '''Returns a set of datetime 2-tuples, each of which represents an interval
           of uninterrupted time when the target was above the horizon (or below, if
           up=False). The set of tuples gives the complete target down intervals
           between the Visibility object's start and end date.
        '''
        effective_horizon = set_airmass_limit(airmass, self.horizon.in_degrees())

        star = Star(self.site['latitude'], target, effective_horizon)

        if up:
            day_interval_func = self.find_when_target_is_up
        else:
            day_interval_func = self.find_when_target_is_down

        # Find rise/set/transit for each day
        intervals = []
        current_date = self.start_date
        while current_date < self.end_date:
            one_day_intervals = day_interval_func(target, current_date, star, airmass)

            # Add today's intervals to the accumulating list of intervals
            intervals.extend(one_day_intervals)

            # Move on to tomorrow
            current_date += ONE_DAY

        # Collapse adjacent intervals into continuous larger intervals
        intervals = coalesce_adjacent_intervals(intervals)

        return intervals


    def find_when_target_is_down(self, target, dt, star=None, airmass=None):
        '''Returns a single datetime 2-tuple, representing an interval
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
        _log.debug("dt: %s", dt)
        up_intervals = self.find_when_target_is_up(target, dt, star, airmass)

        if not up_intervals:
            _log.warn("Got no up intervals!")
            _log.warn("dt was: %s", dt)
            _log.warn("target was: %s", target)
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


    def find_when_target_is_up(self, target, dt, star=None, airmass=None):
        '''Returns a single datetime 2-tuple, representing an interval
           of uninterrupted time above the horizon at the specified
           site, for the requested date.

           Note: Even though this function currently ignores times, the dt
           object must be a datetime, *not* a date.

           TODO: Clean up this unpleasant interface - no star need be passed if
           the sun is the target, but one *must* be passed otherwise. This isn't
           obvious from the method signature.
        '''
        effective_horizon = set_airmass_limit(airmass, self.horizon.in_degrees())

        # Remove any time component of the provided datetime object
        dt = dt.replace(hour=0, minute=0, second=0, microsecond=0)

        # Get the rise/set/transit times for the target, for this day
        intervals = []
        # TODO: Catch the RiseSetError for circumpolar stars
        # TODO: Return either a complete or empty interval, as appropriate
        # TODO: This requires either introspecting state in the error, or
        # TODO: calling an is_circumpolar method before this calculation

        if target == 'sun':
            transits, rises, sets = calc_sunrise_set(self.site, dt, self.twilight)
        else:
            # Test for circumpolarity
            if star.is_always_up(dt):
                # Return a full interval over the entire day
                intervals.append((dt, dt + ONE_DAY))
                return intervals

            # Catch target never rising
            try:
                effective_horizon_angle = Angle(degrees=effective_horizon)
                transits, rises, sets = calc_rise_set(target, self.site,
                                                      dt, effective_horizon_angle)
            except RiseSetError:
                return intervals


        _log.debug("latitude: %s",     self.site['latitude'].in_degrees())
        _log.debug("longitude: %s",    self.site['longitude'].in_degrees())
        _log.debug("twilight: %s",     self.twilight)
        _log.debug("dt: %s",           dt)
        _log.debug("rise: %s (%s)",    rises, dt + rises)
        _log.debug("transit: %s (%s)", transits, dt + transits)
        _log.debug("set: %s (%s)",     sets, dt + sets)

        # Case 1: Overlapping start of day boundary
        # Target rose yesterday, and sets today. Rises again later today.
        #         |       x                                     |
        #         |    x     x                                  |
        #         | x           x                               | x
        #        x|                x                           x|
        #     x   |                   x                     x   |
        #   Rise  0hr  Transit       Set                  Rise 24hr
        if (rises > transits) and (sets > transits):

            # Store the first interval - start of day until target set
            absolute_set = sets + dt
            intervals.append((dt, absolute_set))

            # Store the second interval - target rise until end of day
            absolute_rise = rises + dt
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
        elif (rises < transits) and (sets > transits):
            # Only one interval - rise until target set
            absolute_rise = rises + dt
            absolute_set  = sets + dt
            intervals.append((absolute_rise, absolute_set))


        # Case 3: Overlapping end of day boundary
        # Target rose yesterday, and sets today. Rises again later today.
        #                 x    |                                        x    |
        #              x     x |                                     x     x |
        #           x          |x                                 x          |x
        #        x             |   x                           x             |   x
        #     x                |      x                     x                |
        #   Rise       Transit 0hr   Set                  Rise      Transit 24hr
        elif (rises < transits) and (sets < transits):
            # Same code as case 1!
            # Store the first interval - start of day until target set
            absolute_set = sets + dt
            intervals.append((dt, absolute_set))

            # Store the second interval - target rise until end of day
            absolute_rise = rises + dt
            intervals.append((absolute_rise, dt + ONE_DAY))


        return intervals



    def __key(self):
        return (self.site, self.start_date, self.end_date, self.horizon, self.twilight)

    def __eq__(self, other):
        return self.__key() == other.__key()

    def __hash__(self):
        return hash(self.__key())
