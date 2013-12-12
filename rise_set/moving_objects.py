#!/usr/bin/env python

'''
moving_objects.py - Routines for moving object (asteroids, comets)

description

Authors: Eric Saunders
         Tim Lister
December 2013
'''

from __future__ import division

from rise_set.astrometry import date_to_tdb, calc_local_hour_angle, calculate_altitude
from rise_set.angle      import Angle
from rise_set.visibility import coalesce_adjacent_intervals

import slalib as sla

import ast
from datetime import timedelta

def initialelemdict():
    '''Create an inital empty orbital elements dictionary'''
    keys  = "Name H G epoch mean_anomaly long_node arg_perihelion inclination "
    keys += "eccentricity MDM semi_axis n_obs n_nights"

    return {}.fromkeys(keys.split())


def extract_mpc_epoch(epochstring):
    '''Convert packed MPC epoch format (e.g. 'J974L') from NEOCP orbit files
    into a datetime.datetime epoch (e.g. '1997 4 21'). Returns -1 if invalid
    length (no other sanity checking is done)'''

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

        elements['Name']           = chunks[0]
        elements['H']              = float(chunks[1])
        elements['G']              = float(chunks[2])
        elements['epoch']          = extract_mpc_epoch(chunks[3])
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
            elements['n_obs']    = int(chunks[13])
            elements['n_nights'] = int(chunks[14]) + 1

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

    def __init__(self):
        self.sites = {}
        self.site_lat_longs = []


    def add_sites(self, site_dicts):
        for site_dict in site_dicts:
            self.add_site(site_dict)

        return


    def add_site(self, site_dict):
        lat_long = (site_dict['latitude'], site_dict['longitude'])

        if lat_long not in self.site_lat_longs:
            self.site_lat_longs.append(lat_long)
            site_dict['latitude']  = Angle(degrees=site_dict['latitude'])
            site_dict['longitude'] = Angle(degrees=site_dict['longitude'])
            site_dict['horizon']   = Angle(degrees=site_dict['horizon'])
            self.sites[site_dict['name']] = site_dict

        return


    def get(self, site_name):
        return self.sites.get(site_name)


    def __iter__(self):
        for site_dict in self.sites.values():
            yield site_dict



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



def chunk_windows(window, chunksize):
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
    sites = initialise_sites(site_filename)
    up_intervals_at = {}
    for site in sites:
        up_intervals, altitudes = find_moving_object_up_intervals(window, elements, site,
                                                                  chunksize)
        coalesced_intervals = coalesce_adjacent_intervals(up_intervals)
        up_intervals_at[site['name']] = coalesced_intervals

    return up_intervals_at


def find_moving_object_up_intervals(window, elements, site, chunksize):
    coords = calc_ephemerides(window, elements, site, chunksize)

    intervals = []
    altitudes = []
    for coord in coords:
        local_hour_angle = calc_local_hour_angle(coord['ra_app'],
                                                 site['longitude'],
                                                 coord['start'])
        altitude = calculate_altitude(site['latitude'].in_degrees(),
                                      coord['dec_app'].in_degrees(),
                                      local_hour_angle)

        interval = (coord['start'], coord['end'])
        altitudes.append(altitude)
        intervals.append(interval)


    up_intervals = []
    up_altitudes = []
    for i in xrange(len(intervals)-1):
        # Keep the interval only if both it and the next are up
        # (in other words, throw away partial intervals)
        # This works because we have an extra pseudo-interval with size 0
        if ( ( altitudes[i].in_degrees() > site['horizon'].in_degrees()   ) and
             ( altitudes[i+1].in_degrees() > site['horizon'].in_degrees() ) ):

           up_intervals.append(intervals[i])
           up_altitudes.append(altitudes[i])


    return up_intervals, up_altitudes


def calc_ephemerides(interval, elements, site, chunksize=timedelta(minutes=15)):
    chunked_intervals = chunk_windows(interval, chunksize)

    coords = []
    for start, end in chunked_intervals:
        ra_app, dec_app = elem_to_topocentric_apparent(start, elements, site)
        ephem = {
                  'ra_app'  : ra_app,
                  'dec_app' : dec_app,
                  'start'   : start,
                  'end'     : end,
                }
        coords.append(ephem)

    start, end = chunked_intervals[-1]
    ra_app, dec_app = elem_to_topocentric_apparent(end, elements, site)
    ephem = {
              'ra_app'  : ra_app,
              'dec_app' : dec_app,
              'start'   : end,
              'end'     : end,
            }
    coords.append(ephem)

    return coords


class MovingViolation(object):
    '''Exception for moving object errors.'''
    pass

