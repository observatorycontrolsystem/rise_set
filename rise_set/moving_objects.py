#!/usr/bin/env python

'''
moving_objects.py - Routines for moving object (asteroids, comets) calculations

description

Authors: Eric Saunders
         Tim Lister
December 2013
'''

from __future__ import division

from rise_set.astrometry import (gregorian_to_ut_mjd, date_to_tdb,
                                 calc_local_hour_angle, calculate_altitude)
from rise_set.angle      import Angle
from rise_set.utils      import coalesce_adjacent_intervals

import slalib as sla

import ast
from datetime import datetime,timedelta

def initialelemdict():
    '''Create an inital empty orbital elements dictionary for use by read_neocp_orbit.'''
    keys  = "name H G epoch mean_anomaly long_node arg_perihelion inclination "
    keys += "eccentricity MDM semi_axis n_obs n_nights type"

    return {}.fromkeys(keys.split())


def extract_mpc_epoch(epochstring):
    '''Convert packed MPC epoch format (e.g. 'J974L') from NEOCP orbit files
    into a datetime.datetime epoch (e.g. '1997 4 21'). Returns -1 if invalid
    length (no other sanity checking is done).'''

    if len(epochstring) != 5: return -1
    year = 100 * (ord(epochstring[0]) - ord('A') + 10) + \
           10 * (ord(epochstring[1]) - ord('0')) +       \
           (ord(epochstring[2]) - ord('0'))

    arr = [0, 0]
    for i in range(0,2,+1):
        if (ord(epochstring[i + 3]) >= ord('A')):
            arr[i] = ord(epochstring[i + 3]) - ord('A') + 10
        else:
            arr[i] = ord(epochstring[i + 3]) - ord('0')

    return datetime(year, arr[0], arr[1], 0)


def read_neocp_orbit(orbfile):
    '''Read an NEOCP orbit file, storing the contents as a dict.'''

    elements = initialelemdict()

    orb_fh = open(orbfile, 'r')
    for line in orb_fh.readlines():
        if line.lstrip()[0] == '<': continue
        if line[0:7] == ' REMARK' : continue
        if not len(line.strip())  : continue

        line   = line.rstrip()
        chunks = line.split()

        elements['type']           = 'MPC_MINOR_PLANET'  # TODO: Hard-coded to asteroids for now
        elements['name']           = chunks[0]
        elements['H']              = float(chunks[1])
        elements['G']              = float(chunks[2])

        epoch_dt                   = extract_mpc_epoch(chunks[3])
        elements['epoch']          = gregorian_to_ut_mjd(epoch_dt)

        elements['mean_anomaly']   = Angle(degrees=float(chunks[4]))
        elements['arg_perihelion'] = Angle(degrees=float(chunks[5]))
        elements['long_node']      = Angle(degrees=float(chunks[6]))
        elements['inclination']    = Angle(degrees=float(chunks[7]))
        elements['eccentricity']   = float(chunks[8])
        elements['MDM']            = Angle(degrees=float(chunks[9]))
        elements['semi_axis']      = float(chunks[10])

        # If it's an NEOCP-style line...
        if len(chunks) == 17:
            elements['n_obs']    = int(chunks[11])
            elements['n_nights'] = int(chunks[13]) + 1

        # ...or if it's an MPCORB-style line...
        elif len(chunks) == 25 or len(chunks) == 23:
            nobs_field = 12
            if len(chunks) == 25:
                nobs_field = 13
            elements['n_obs']    = int(chunks[nobs_field])
            try:
                num_nights_or_years = int(chunks[14]) + 1
            except ValueError:
                num_nights_or_years = chunks[14]

            elements['n_nights'] = num_nights_or_years

        # Otherwise complain
        else:
            print "Unknown number of columns"
            elements['n_obs']    = 0
            elements['n_nights'] = 0

        # We only want the first non-header line
        break


    orb_fh.close()
    return elements



def file_to_dicts(filename):
    with open(filename, 'r') as fh:
        data = fh.read()
        return ast.literal_eval(data)


def initialise_sites(site_filename):
    '''Factory method to create a Sites object from a telescopes.dat file.'''
    site_dicts = file_to_dicts(site_filename)

    sites = Sites()
    sites.add_sites(site_dicts)

    return sites


class Sites(object):
    '''Collection for a group of site dictionaries. Only the first site with a
       given lat/long is stored. This allows us to conveniently read a
       file with a list of telescopes, which may share the same site.'''

    def __init__(self):
        self.sites = {}
        self.site_lat_longs = []


    def add_sites(self, site_dicts):
        for site_dict in site_dicts:
            self.add_site(site_dict)

        return


    def add_site(self, site_dict):
        lat_long = (site_dict['latitude'], site_dict['longitude'])
        HOURS_TO_DEGREES = 15.0

        if lat_long not in self.site_lat_longs:
            self.site_lat_longs.append(lat_long)
            site_dict['latitude']         = Angle(degrees=site_dict['latitude'])
            site_dict['longitude']        = Angle(degrees=site_dict['longitude'])
            site_dict['horizon']          = Angle(degrees=site_dict['horizon'])
            site_dict['ha_limit_neg']     = Angle(degrees=site_dict['ha_limit_neg'] * HOURS_TO_DEGREES)
            site_dict['ha_limit_pos']     = Angle(degrees=site_dict['ha_limit_pos'] * HOURS_TO_DEGREES)
            self.sites[site_dict['name']] = site_dict

        return


    def get(self, site_name):
        return self.sites.get(site_name)


    def __iter__(self):
        for site_dict in self.sites.values():
            yield site_dict



def elem_to_topocentric_apparent(dt, elements, site, JFORM=2):
    '''Given a datetime, set of MPC orbital elements and a site, return the
       apparent topocentric RA/Dec. This is what you'd use for a rise/set
       calculation, for example.
       JFORM should be set to 2 (default) for asteroids/minor planets and
       to 3 for comets'''
    tdb = date_to_tdb(dt)

    MINOR_PLANET_JFORM = 2
    COMET_JFORM = 3 
    MDM_PLACEHOLDER    = 0.0  # Only used for major planets
    MEANANOM_PLACEHOLDER    = 0.0  # Not applicable for comets

    status = 0
    if JFORM == MINOR_PLANET_JFORM:
# Minor planets (asteroids)
        ra_app_rads, dec_app_rads, earth_obj_dist, status = sla.sla_plante(
                                                    tdb,
                                                    site['longitude'].in_radians(),
                                                    site['latitude'].in_radians(),
                                                    JFORM,
                                                    elements['epoch'],
                                                    elements['inclination'].in_radians(),
                                                    elements['long_node'].in_radians(),
                                                    elements['arg_perihelion'].in_radians(),
                                                    elements['semi_axis'],
                                                    elements['eccentricity'],
                                                    elements['mean_anomaly'].in_radians(),
                                                    MDM_PLACEHOLDER,
                                                  )
    elif JFORM == COMET_JFORM:
# Comets
        ra_app_rads, dec_app_rads, earth_obj_dist, status = sla.sla_plante(
                                                    tdb,
                                                    site['longitude'].in_radians(),
                                                    site['latitude'].in_radians(),
                                                    JFORM,
                                                    elements['epochofperih'],
                                                    elements['inclination'].in_radians(),
                                                    elements['long_node'].in_radians(),
                                                    elements['arg_perihelion'].in_radians(),
                                                    elements['perihdist'],
                                                    elements['eccentricity'],
                                                    MEANANOM_PLACEHOLDER,
                                                    MDM_PLACEHOLDER,
                                                  )
    else:
        status = -1

    error = {
               0 : 'OK',
              -1 : 'illegal JFORM',
              -2 : 'illegal eccentricity',
              -3 : 'illegal mean distance',
              -4 : 'illegal mean daily motion',
              -5 : 'numerical error',
            }

    if (status != 0):
        elem_string = 'Bad Elements:\n'
        for key in elements.keys():
            elem_string += key + ' = ' + str(elements[key]) + '\n'
        print elem_string
        raise MovingViolation('Error: ' + str(status) + ' (' + error[status] + ')')

    return Angle(radians=ra_app_rads), Angle(radians=dec_app_rads)



def chunk_windows(window, chunksize):
    '''Given a user window, split that window into interval tuples (start, end)
       of size chunksize.'''
    chunk_start = window['start']
    chunk_end   = chunk_start + chunksize

    intervals = []
    while chunk_end <= window['end']:
        interval = (chunk_start, chunk_end)
        intervals.append(interval)

        chunk_start += chunksize
        chunk_end   += chunksize

    # Pick up any partial remaining interval
    if chunk_start < window['end']:
        interval = (chunk_start, window['end'])
        intervals.append(interval)

    return intervals


def find_moving_object_network_up_intervals(window, elements, site_filename, chunksize):
    '''Network wide wrapper for find_moving_object_up_intervals.'''

    sites = initialise_sites(site_filename)
    up_intervals_at = {}
    for site in sites:
        up_intervals, altitudes = find_moving_object_up_intervals(window, elements, site,
                                                                  chunksize)
        coalesced_intervals = coalesce_adjacent_intervals(up_intervals)
        up_intervals_at[site['name']] = coalesced_intervals

    return up_intervals_at


def hour_angle_within_limits(hour_angle, neg_limit, pos_limit):
   if ( hour_angle.in_degrees() > neg_limit.in_degrees()  and
        hour_angle.in_degrees() < pos_limit.in_degrees()):

       return True

   return False


def ephemeris_chunk_within_ha_limits(ha1, ha2, neg_limit, pos_limit):
    if ( hour_angle_within_limits(ha1, neg_limit, pos_limit) and
         hour_angle_within_limits(ha2, neg_limit, pos_limit) ):
           return True

    return False


def ephemeris_chunk_above_horizon(alt1, alt2, horizon):
    if (  alt1.in_degrees() > horizon.in_degrees()  and
          alt2.in_degrees() > horizon.in_degrees()  ):
        return True

    return False


def find_moving_object_up_intervals(window, elements, site, chunksize=timedelta(minutes=15)):
    '''Return only the (interval, altitude) pairs for which the moving object is
       above the horizon and within the hour angle limits at the provided site.
       If you just want rising and setting visible intevals, supply hour angle limits
       of -12 and 12 in your site dict.
    '''

    # Map the type of moving object to SLALIB numeric convention
    target_type = elements['type'].lower()
    if target_type == 'mpc_minor_planet':
        jform = 2
    elif target_type == 'mpc_comet':
        jform = 3
    else:
        raise MovingViolation("Unsupported target type: '%s'" % str(target_type))

    coords = calc_ephemerides(window, elements, site, chunksize, jform)

    intervals   = []
    altitudes   = []
    hour_angles = []
    for coord in coords:
        local_hour_angle = calc_local_hour_angle(coord['ra_app'],
                                                 site['longitude'],
                                                 coord['start'])
        altitude = calculate_altitude(site['latitude'].in_degrees(),
                                      coord['dec_app'].in_degrees(),
                                      local_hour_angle.in_degrees())

        interval = (coord['start'], coord['end'])
        altitudes.append(altitude)
        intervals.append(interval)
        hour_angles.append(local_hour_angle)

    up_intervals = []
    up_altitudes = []
    for i in xrange(len(intervals)-1):
        # Keep the interval only if both it and the next are up
        # (in other words, throw away partial intervals)
        # This works because we have an extra pseudo-interval with size 0
        if not ephemeris_chunk_above_horizon(altitudes[i], altitudes[i+1], site['horizon']):
            continue

        if not ephemeris_chunk_within_ha_limits(hour_angles[i], hour_angles[i+1],
                                                site['ha_limit_neg'],
                                                site['ha_limit_pos']):
            continue

        up_intervals.append(intervals[i])
        up_altitudes.append(altitudes[i])

    return up_intervals, up_altitudes


def calc_ephemerides(window, elements, site, chunksize=timedelta(minutes=15), JFORM=2):
    '''Return the apparent RA/Dec for the moving object specified by elements, for a given
       site, every chunksize minutes within the window.'''

    chunked_intervals = chunk_windows(window, chunksize)

    coords = []
    for start, end in chunked_intervals:
        ra_app, dec_app = elem_to_topocentric_apparent(start, elements, site, JFORM)
        ephem = {
                  'ra_app'  : ra_app,
                  'dec_app' : dec_app,
                  'start'   : start,
                  'end'     : end,
                }
        coords.append(ephem)

    start, end = chunked_intervals[-1]
    ra_app, dec_app = elem_to_topocentric_apparent(end, elements, site, JFORM)
    ephem = {
              'ra_app'  : ra_app,
              'dec_app' : dec_app,
              'start'   : end,
              'end'     : end,
            }
    coords.append(ephem)

    return coords


class MovingViolation(Exception):
    '''Exception for moving object errors.'''
    pass

