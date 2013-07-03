#!/usr/bin/python

from angle import Angle



ra = Angle(ra = '02 17 14.0343')
dec = Angle(dec = '+13 43 46.044')

print 'RA (deg):', ra.in_degrees()
print 'dec (deg):', dec.in_degrees()
