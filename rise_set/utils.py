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
    if 'type' in target and target['type'].lower() in ('mpc_minor_planet', 'mpc_comet', 'jpl_major_planet'):
        return True

    return False


def is_static_target(target):
    # Static type targets are treated differently in that their exact intervals are not
    # calculated. They are assumed to be observable when it is nighttime.

    # We don't support explicit ephemeris calculations for satellite targets.
    if 'type' in target and target['type'].lower() == 'satellite':
        return True

    # Hour angle target types are already tied to a particular location and time.
    if 'type' not in target and 'hour_angle' in target:
        return True

    return False


def is_sidereal_target(target):
    if 'type' not in target and 'ra' in target:
        return True

    return False


def target_to_jform(target):
    target_type = target['type'].lower()
    if target_type == 'mpc_minor_planet':
        jform = 2
    elif target_type == 'mpc_comet':
        jform = 3
    elif target_type == 'jpl_major_planet':
        jform = 1
    else:
        raise MovingViolation("Unsupported target type: '{}'".format(target_type))

    return jform


class MovingViolation(Exception):
    '''Exception for moving object errors.'''
    pass
