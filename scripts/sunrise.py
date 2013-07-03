#!/usr/bin/python

from __future__ import division
import datetime
from angle import Angle
import slalib as sla
from astrometry import gregorian_to_ut_mjd


# Date
date = datetime.date(year=2010, month=4, day=27)

# Site
st_andrews = {
    'latitude'  : Angle(degrees=56),
    'longitude' : Angle(degrees=3),
}


ut_mjd = gregorian_to_ut_mjd(date)    
tdb = ut_mjd + (sla.sla_dtt(ut_mjd)/86400.)


(app_ra_rads, app_dec_rads, ang_diameter) = sla.sla_rdplan(tdb, 0, 
                                  st_andrews['longitude'].in_radians(), 
                                  st_andrews['latitude'].in_radians())


print "app_ra_rads:", app_ra_rads
print "app_dec_rads:", app_dec_rads

app_ra  = Angle(radians=app_ra_rads)
app_dec = Angle(radians=app_dec_rads)

print "Apparent RA, Dec (degs):", app_ra.in_degrees(), app_dec.in_degrees()
