#!/usr/bin/env python

'''
darkness.py - When are the dark periods for a given site, over a given timeframe?

This code uses the rise_set library to determine the set of dark periods. It
doesn't do anything clever.

Author: Eric Saunders
February 2011
'''

import datetime
from astrometry import calc_sunrise_set
from visibility import Visibility
from angle      import Angle

# BPL East 0.4m - from http://tlister-linux/django/configdb/telescope/16/
sba = dict(
            latitude  = Angle(degrees=34.4332222222),
            longitude = Angle(degrees=-119.863045833),
           )

date = datetime.datetime(year=2011, month=2, day=9)

(transit, rise, set) = calc_sunrise_set(sba, date, 'sunrise')

print "Transit:", transit
print "Rise:", rise
print "Set:", set



# Question: What are the dark intervals between two given datetimes?




start = datetime.datetime(year=2011, month=2, day=9)
end   = datetime.datetime(year=2011, month=2, day=16)

visibility = Visibility(sba, start, end)
dark_intervals = visibility.get_dark_intervals()

for interval in dark_intervals:
    print "Start:", interval[0]
    print "Stop:",  interval[1]