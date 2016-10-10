#!/usr/bin/env python

'''
utils.py - Low level rise_set utility methods

description

Author: Eric Saunders
        Jason Eastman
December 2013
'''

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


def intersect_intervals(int1, int2):

    ''' Computes the intersections of two sets of datetime 2-tuples, each of which
    represents the overlap between the two tuples.
    '''
    intersect = []
    for start1, end1 in int1:
        for start2, end2 in int2:

            start0 = max(start1,start2)
            end0 = min(end1,end2)

            if start0 < end0:
                intersect.append((start0,end0))

    return intersect


def intersect_many_intervals(*args):

    ''' Generalizes intersect_intervals to an arbitrary number of lists
    e.g., when the object is up, sun is down, and hour angle within limits
    '''
    intersection = args[0]
    for interval in args:
        intersection = intersect_intervals(intersection,interval)

    return intersection


def is_moving_object(target):
    # If a type is not specified, default to sidereal objects
    if 'type' not in target:
        return False

    if target['type'].lower() in ('mpc_minor_planet', 'mpc_comet'):
        return True

    return False


class MovingViolation(Exception):
    '''Exception for moving object errors.'''
    pass
