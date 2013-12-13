#!/usr/bin/env python

'''
utils.py - Low level rise_set utitility methods

description

Author: Eric Saunders
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
