#!/usr/bin/python

from __future__ import division     # floating point division (not integer)

from math import pi
from optparse import OptionParser

opt_parser = OptionParser("Usage: your face")
opt_parser.add_option("-r", "--ra", dest="ra", type="string",
                      help="RA in sexegesimal (e.g. 18 30 45.0)")

(options, args) = opt_parser.parse_args()
if len(args) <1:
    opt_parser.error("Your face is dumb.")

ra_string = options.ra

print "RA string:", ra_string

# Take the ra components provided and store them in a dictionary
if ra_string:
    ra = dict(zip(("hr", "min", "sec"), ra_string.split(None)))

print ra



