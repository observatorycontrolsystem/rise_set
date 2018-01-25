#!/usr/bin/env python

'''
example_rise_set.py - Example usage of rise_set functions.

This example gives some simple demonstrations of how to use the rise_set library to
calculate rise, set and transit times for a star and the sun.

Author: Eric Saunders
August 2011
'''


from rise_set.angle           import Angle
from rise_set.sky_coordinates import RightAscension, Declination
from rise_set.rates           import ProperMotion
from rise_set.astrometry      import calc_rise_set, calc_sunrise_set

import datetime
import logging


def rise_set_of_deneb_from_maui():
    # Target
    deneb = {
                   'ra'  : RightAscension('20 41 25.91'),
                   'dec' : Declination('+45 16 49.22'),
                   # Aladin proper motion units are mas/yr...
                   'ra_proper_motion'  : ProperMotion(RightAscension('00 00 00.00156')),
                   'dec_proper_motion' : ProperMotion(Declination('00 00 00.00155')),
                   'parallax'          : 0.00101,  # Units: arcsec
                   'rad_vel'           : -4.5,  # Units: km/s (-ve approaches)
                   'epoch'             : 2000,
                  }



    # Site (East +ve longitude)
    maui = {
       'latitude'  : Angle(degrees=20.7069444444),
       'longitude' : Angle(degrees=-156.258055556),
    }


    # Date
    date = datetime.datetime(year=2010, month=10, day=25)


    # Do the calculation
    (transit, rise, setting) = calc_rise_set(deneb, maui, date)

    return (transit, rise, setting)



def sunrise_sunset_from_st_andrews():

    # Date
    date = datetime.datetime(year=2010, month=4, day=27)

    # Site (East +ve longitude)
    # Very rough St. Andrews coords, for easy almanac comparison
    st_andrews = {
        'latitude'  : Angle(degrees = 56),
        'longitude' : Angle(degrees = 3),
    }

    # Do the calculation
    (transit, rise, setting) = calc_sunrise_set(st_andrews, date, 'sunrise')

    return (transit, rise, setting)



if __name__ == '__main__':

    # Configure logging
    logging.basicConfig(
        format = "%(levelname)s %(lineno)-6d %(message)s",
        level  = logging.INFO
    )



    (deneb_transit, deneb_rise, deneb_set) = rise_set_of_deneb_from_maui()
    (sun_transit, sun_rise, sun_set)       = sunrise_sunset_from_st_andrews()

    print '''Deneb from Maui:
           Rise:     %s
           Transit : %s
           Set :     %s
          ''' % (deneb_rise, deneb_transit, deneb_set)

    print '''Sun from St Andrews:
           Rise:     %s
           Transit : %s
           Set :     %s
          ''' % (sun_rise, sun_transit, sun_set)
