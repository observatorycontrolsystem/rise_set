#!/usr/bin/env python

'''
moving_objects.py - Routines for moving object (asteroids, comets)

description

Authors: Eric Saunders
         Tim Lister
December 2013
'''

from __future__ import division

from rise_set.astrometry import date_to_tdb
from rise_set.angle      import Angle
import slalib as sla

def elem_to_topocentric_apparent(dt, elements, site):
    '''Given a datetime, set of MPC orbital elements and a site, return the
       apparent topocentric RA/Dec. This is what you'd use for a rise/set
       calculation, for example.'''
    tdb = date_to_tdb(dt)

    # Slalib convention - force minor planets only
    MINOR_PLANET_JFORM = 2
    MDM_PLACEHOLDER    = 0.0  # Only used for major planets

    epoch_mjd = date_to_tdb(elements['epoch'])

    status = 0
    ra_app_rads, dec_app_rads, earth_obj_dist, status = sla.sla_plante(
                                                    tdb,
                                                    site['longitude'].in_radians(),
                                                    site['latitude'].in_radians(),
                                                    MINOR_PLANET_JFORM,
                                                    epoch_mjd,
                                                    elements['inclination'].in_radians(),
                                                    elements['long_node'].in_radians(),
                                                    elements['arg_perihelion'].in_radians(),
                                                    elements['semi_axis'],
                                                    elements['eccentricity'],
                                                    elements['mean_anomaly'].in_radians(),
                                                    MDM_PLACEHOLDER,
                                                  )

    error = {
              0 : 'OK',
              1 : 'illegal JFORM',
              2 : 'illegal eccentricity',
              3 : 'illegal mean distance',
              4 : 'illegal mean daily motion',
              5 : 'numerical error',
            }

    if (status != 0):
        raise MovingViolation('Error:' + error[status])


    return Angle(radians=ra_app_rads), Angle(radians=dec_app_rads)



class MovingViolation(object):
    '''Exception for moving object errors.'''
    pass

