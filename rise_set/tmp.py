#!/usr/bin/python

from angle import Angle
from astrometry import *
import datetime
import logging

def test_deneb_from_maui():
    # Target
    deneb = {
                   'ra'  : Angle(),
                   'dec' : Angle(),
                   'ra_proper_motion'  : Angle(),
                   'dec_proper_motion' : Angle(),
                   'parallax'          : 0.00101,  # Units: arcsec
                   'rad_vel'           : -4.5,  # Units: km/s (-ve approaches)
                   'epoch'             : 2000,
                  }


    deneb['ra'].from_sexegesimal('20 41 25.91', ra=True)
    deneb['dec'].from_sexegesimal('+45 16 49.22', dec=True)

    # Aladin units are mas/yr...
    deneb['ra_proper_motion'].from_sexegesimal('00 00 00.00156', ra=True)
    deneb['dec_proper_motion'].from_sexegesimal('00 00 00.00155', dec=True)


    # Site (East +ve longitude)
    maui = {
       'latitude'  : 20.7069444444,
       'longitude' : -156.258055556
    }


    # Date
    date = datetime.date(year=2010, month=10, day=25)



    (transit, rise, setting) = calc_rise_set(deneb, maui, date)

    return (transit, rise, setting)


def test_canopus_from_siding_spring():

    # Target
    canopus = {
                 'ra'                : Angle(),
                 'dec'               : Angle(),
                 'ra_proper_motion'  : Angle(),
                 'dec_proper_motion' : Angle(),
                 'parallax'          : 0.01043,   # Units: arcsec
                 'rad_vel'           : 20.5,      # Units: km/s (-ve approaches)
                 'epoch'             : 2000,
               }

    canopus['ra'].from_sexegesimal('06 23 57.11', ra=True)
    canopus['dec'].from_sexegesimal('-52 40 03.5', dec=True)


    # Aladin units are mas/yr...
    canopus['ra_proper_motion'].from_sexegesimal('00 00 00.01999', ra=True)
    canopus['dec_proper_motion'].from_sexegesimal('00 00 00.02367', dec=True)


    # Site (East +ve longitude)
    siding_spring = {
        'latitude'  : -31.273,
        'longitude' : 149.070593
    }


    # Date
    date = date.datetime(year=2010, month=3, day=12)

    (transit, rise, setting) = calc_rise_set(canopus, siding_spring, date)

    return (transit, rise, setting)



def test_canopus_from_st_andrews():

    # Target
    canopus = {
                 'ra'                : Angle(),
                 'dec'               : Angle(),
                 'ra_proper_motion'  : Angle(),
                 'dec_proper_motion' : Angle(),
                 'parallax'          : 0.01043,   # Units: arcsec
                 'rad_vel'           : 20.5,      # Units: km/s (-ve approaches)
                 'epoch'             : 2000,
               }

    canopus['ra'].from_sexegesimal('06 23 57.11', ra=True)
    canopus['dec'].from_sexegesimal('-52 40 03.5', dec=True)


    # Aladin units are mas/yr...
    canopus['ra_proper_motion'].from_sexegesimal('00 00 00.01999', ra=True)
    canopus['dec_proper_motion'].from_sexegesimal('00 00 00.02367', dec=True)


    # Site (East +ve longitude)
    # Very rough St. Andrews coords, for easy almanac comparison
    st_andrews = {
        'latitude'  : 56,
        'longitude' : 3
    }


    # Date
    date = date.datetime(year=2010, month=3, day=12)

    (transit, rise, setting) = calc_rise_set(canopus, st_andrews, date)

    return (transit, rise, setting)




def test_capella_from_st_andrews():

    # Target
    capella = {
                 'ra'                : Angle(),
                 'dec'               : Angle(),
                 'ra_proper_motion'  : Angle(),
                 'dec_proper_motion' : Angle(),
                 'parallax'          : 0.07729,   # Units: arcsec
                 'rad_vel'           : 30.2,      # Units: km/s (-ve approaches)
                 'epoch'             : 2000,
               }

    capella['ra'].from_sexegesimal('05 16 41.36', ra=True)
    capella['dec'].from_sexegesimal('+45 59 52.8', dec=True)


    # Aladin units are mas/yr...
    capella['ra_proper_motion'].from_sexegesimal('00 00 00.07552', ra=True)
    capella['dec_proper_motion'].from_sexegesimal('-00 00 00.42711', dec=True)


    # Site (East +ve longitude)
    # Very rough St. Andrews coords, for easy almanac comparison
    st_andrews = {
        'latitude'  : 56,
        'longitude' : 3
    }


    # Date
    date = datetime.date(year=2010, month=3, day=12)

    (transit, rise, setting) = calc_rise_set(capella, st_andrews, date)

    return (transit, rise, setting)



def get_sunrise():
    
    # Date
    date = datetime.date(year=2010, month=4, day=27)

    # Site (East +ve longitude)
    # Very rough St. Andrews coords, for easy almanac comparison
    st_andrews = {
        'latitude'  : 56,
        'longitude' : 3
    }

    (transit, rise, setting) = calc_sunrise_set(st_andrews, date, 'sunrise')

    return (transit, rise, setting)


# Configure logging
logging.basicConfig(
    format = "%(levelname)s %(lineno)-6d %(message)s",
    level  = logging.DEBUG
)



#(transit, rise, setting) = test_deneb_from_maui()
#(transit, rise, setting) = test_canopus_from_st_andrews()
#(transit, rise, setting) = test_capella_from_st_andrews()
(transit, rise, setting) = get_sunrise()

print               
print __name__
print 'Rise:', rise
print 'Transit:', transit
print 'Set:', setting

