#!/usr/bin/python

'''rise_set/astrometry.py - astrometry calculations for telescope scheduling.

This package provides routines for finding the positions of astronomical bodies
to reasonable precision. It can be used to calculate target uptime and sunrise,
sunset and twilight times.

The code uses SLALIB for the heavy lifting. The rise/set/transit algorithms here
are implementations of Astronomical Algorithms, Ch. 14 (Jean Meeus).


Author: Eric Saunders (esaunders@lcogt.net)

May 2010
'''

# Required for true (non-integer) division
from __future__ import division

# Standard library imports
# Trigonometry functions
from math import sin, cos, asin, acos, radians, degrees, modf
import datetime

# API for accessing data files after deployment, within an egg
from pkg_resources import resource_stream

# Third party imports
import slalib as sla
from angle import Angle

# Import logging modules
import logging
import logging.config


# Configure logger from config file - '_' stops names being exported
_log_config_file = 'logging.conf'
_log_config_location = resource_stream(__name__, _log_config_file)
logging.config.fileConfig(_log_config_location)
_log = logging.getLogger('rise_set.astrometry')


# Set convenient constants
ONE_DAY = datetime.timedelta(days=1)
MIDNIGHT = datetime.time()

def gregorian_to_ut_mjd(date):
    '''Convert Gregorian calendar date to UT MJD.'''

    error = {
            0 : 'OK',
            1 : 'bad year (MJD not computed)',
            2 : 'bad month (MJD not computed)',
            3 : 'bad day (MJD not computed)',
    }
    status = 0

    (mjd, status) = sla.sla_caldj(date.year, date.month, date.day)

    if (status != 0):
        raise InvalidDateError('Error:' + error[status])

    return mjd



def ut_mjd_to_gmst(mjd):
    '''Convert UT MJD to Greenwich mean sidereal time (an Angle).'''

    gmst_in_radians = sla.sla_gmst(mjd)

    return Angle(radians=gmst_in_radians)



def calc_apparent_sidereal_time(date):
    '''Return the apparent sidereal time (an Angle) for a given date.'''

    # Convert Gregorian to UT MJD
    mjd = gregorian_to_ut_mjd(date)

    # Convert UT MJD to Greenwich mean sidereal time (in radians)
    gmst = ut_mjd_to_gmst(mjd)

    # Find the apparent sidereal time (GAST)
    # GAST = GMST + Eqn. of equinoxes
    gast_in_rads = gmst.in_radians() + sla.sla_eqeqx(mjd)

    return Angle(radians=gast_in_rads)



def mean_to_apparent(target, tdb):
    '''Given a target and TDB, return an apparent (RA, Dec) tuple.
       Thin wrapper for SLA_MAP.
    '''

    (ra_app_rads, dec_app_rads) = sla.sla_map(target['ra'].in_radians(),
                                              target['dec'].in_radians(),
                                              target['ra_proper_motion'].in_radians(),
                                              target['dec_proper_motion'].in_radians(),
                                              target['parallax'],
                                              target['rad_vel'],
                                              target['epoch'],
                                              tdb)

    ra_apparent  = Angle(radians=ra_app_rads)
    dec_apparent = Angle(radians=dec_app_rads)

    return (ra_apparent, dec_apparent)



def calc_rise_set_hour_angle(latitude_degs, dec_apparent, std_altitude):
    '''Find the hour angle H_0 corresponding to the time of rise or set of
       a  celestial body, according to

                  sin(h_0) - [ sin(phi) sin(delta_2) ]
       cos H_0 =  -----------------------------------
                         cos(phi) cos(delta_2)
    '''


    latitude = Angle(degrees=latitude_degs)

    # Evaluate the numerator
    top    = ( sin(std_altitude.in_radians())
               - sin(latitude.in_radians())
               * sin(dec_apparent.in_radians()) )

    # Evaluate the denominator
    bottom = cos(latitude.in_radians()) * cos(dec_apparent.in_radians())

    # Evaluate the fraction
    cos_h_0 = top/bottom

    # TODO: Return status object, with is_circumpolar() etc. methods
    # Throw exceptions if the rise/set accessor on this object is called

    # Sanity check the result
    if ( cos_h_0 > 1 ):
        return (None, 'Target never rises at this latitude')

    elif ( cos_h_0 < -1):
        return (None, 'Target never sets at this latitude')


    # Extract and return the hour angle
    hour_angle_in_rads = acos( cos_h_0 )
    hour_angle = Angle(radians=hour_angle_in_rads)

    return (hour_angle, 'OK')



def calc_transit_day_fraction(ra, longitude_in_degrees, app_sidereal_time):
    '''Find the time, expressed as a fraction of a day, when the transit occurs.

       Note: Sign of longitude is reversed relative to Astro Alg., p.98, because
       SLALIB uses East +ve in all functions except SLA_OBS, but Astro. Alg.
       uses West +ve. <Sigh>

              alpha_2 - L - theta_0
       m_0 = ------------------------
                      360
    '''

    m_0 = (ra.in_degrees() - longitude_in_degrees
            - app_sidereal_time.in_degrees()) / 360

    m_0 = normalise_day(m_0)

    return m_0



def calc_rising_day_fraction(m_0, hour_angle):
    '''Find the time, expressed as a fraction of a day, when the target rises.

                     H_0
       m_1 = m_0 - -------
                     360
    '''

    m_1 = m_0 - (hour_angle.in_degrees()/360)
    m_1 = normalise_day(m_1)

    return m_1



def calc_setting_day_fraction(m_0, hour_angle):
    '''Find the time, expressed as a fraction of a day, when the target sets.

                     H_0
       m_2 = m_0 + -------
                     360
    '''

    m_2 = m_0 + (hour_angle.in_degrees()/360)
    m_2 = normalise_day(m_2)

    return m_2



def normalise_day(day_frac):
    '''Adjust the day fraction to correspond to the current day.
      day_frac is a fractional day, so should be adjusted to fall in the range
      0-1 if necessary (Astro. Alg. p.98)
      TODO: This may not be what you want, if you are concerned about the next
      rise or set, as opposed to one that happened today.
    '''

    if ( day_frac < 0 ):
        day_frac += 1

    if ( day_frac > 1 ):
        day_frac -= 1

    return day_frac



def calc_rise_set(target, site, date):
    '''Return a string tuple (transit, rise, set) of times. Each time is a tuple
       of (hour, minute, second) values.
    '''

    std_alt_of_stars = Angle(degrees=-0.5667)

    ut_mjd = gregorian_to_ut_mjd(date)
    tdb = ut_mjd + (sla.sla_dtt(ut_mjd)/86400)


    app_ra, app_dec   = mean_to_apparent(target, tdb)
    app_sidereal_time = calc_apparent_sidereal_time(date)
    (hour_angle, msg) = calc_rise_set_hour_angle(site['latitude'], app_dec,
                                                 std_alt_of_stars)

    if ( not hour_angle ):
        raise RiseSetError(msg)

    m_0 = calc_transit_day_fraction(app_ra, site['longitude'], app_sidereal_time)
    m_1 = calc_rising_day_fraction(m_0, hour_angle)
    m_2 = calc_setting_day_fraction(m_0, hour_angle)

    transit = day_frac_to_hms(m_0)
    rise    = day_frac_to_hms(m_1)
    set     = day_frac_to_hms(m_2)

    _log.info('Rise time - unrefined (h, m, s): %s' % (rise,))
    _log.info('Transit time - unrefined (h, m, s): %s' % (transit,))
    _log.info('Set time - unrefined (h, m, s): %s' % (set,))

    (m_0, m_1, m_2) = refine_day_fraction(app_sidereal_time, m_0, m_1, m_2, tdb,
                                         target, site, std_alt_of_stars)

    transit = day_frac_to_hms(m_0)
    rise    = day_frac_to_hms(m_1)
    set     = day_frac_to_hms(m_2)


    return (transit, rise, set)



def calc_sunrise_set(site, date, twilight):
    '''Return a string tuple (transit, rise, set) of times. Each time is a tuple
       of (hour, minute, second) values.
    '''

    ut_mjd = gregorian_to_ut_mjd(date)
    tdb = ut_mjd + (sla.sla_dtt(ut_mjd)/86400)

    target = dict(planet = 'sun')
    (app_ra, app_dec) = apparent_planet_pos(target['planet'], tdb, site)

    _log.info("RA, Dec (apparent, degrees) for %s: (%s, %s)" 
                 % (target['planet'], app_ra.in_degrees(), app_dec.in_degrees()))

    app_sidereal_time = calc_apparent_sidereal_time(date)

    sun_std_alt = {
                     'sunrise'           : Angle(degrees=-5/6),
                     'sunset'            : Angle(degrees=-5/6),
                     'civil'             : Angle(degrees=-6),
                     'nautical'          : Angle(degrees=-12),
                     'astronomical'      : Angle(degrees=-18)
                    }


    (hour_angle, msg) = calc_rise_set_hour_angle(site['latitude'], app_dec,
                                                 sun_std_alt[twilight])

    if ( not hour_angle ):
        raise RiseSetError(msg)

    m_0 = calc_transit_day_fraction(app_ra, site['longitude'], app_sidereal_time)
    m_1 = calc_rising_day_fraction(m_0, hour_angle)
    m_2 = calc_setting_day_fraction(m_0, hour_angle)

    transit = day_frac_to_hms(m_0)
    rise    = day_frac_to_hms(m_1)
    set     = day_frac_to_hms(m_2)

    _log.info('Rise time - unrefined (h, m, s): %s' % (rise,))
    _log.info('Transit time - unrefined (h, m, s): %s' % (transit,))
    _log.info('Set time - unrefined (h, m, s): %s' % (set,))

    (m_0, m_1, m_2) = refine_day_fraction(app_sidereal_time, m_0, m_1, m_2, tdb,
                                          target, site, sun_std_alt[twilight])

    transit = day_frac_to_hms(m_0)
    rise    = day_frac_to_hms(m_1)
    set     = day_frac_to_hms(m_2)


    return (transit, rise, set)



def apparent_planet_pos(planet_name, tdb, site):
    '''Return the topocentric apparent position (ra, dec) tuple of a planet at
       a particular time, from a particular site.
       Thin wrapper for SLA_RDPLAN.
    '''

    latitude  = Angle(degrees=site['latitude'])
    longitude = Angle(degrees=site['longitude'])

    # TODO: Raise an error if an invalid body is provided.

    planet = dict(
                   mercury = 1,
                   venus   = 2,
                   moon    = 3,
                   mars    = 4,
                   jupiter = 5,
                   saturn  = 6,
                   uranus  = 7,
                   neptune = 8,
                   pluto   = 9,
                   sun     = 0
                  )


    (app_ra_rads, app_dec_rads, ang_diameter) = sla.sla_rdplan(tdb,
                                                        planet[planet_name],
                                                        longitude.in_radians(),
                                                         latitude.in_radians())

    app_ra  = Angle(radians=app_ra_rads)
    app_dec = Angle(radians=app_dec_rads)

    return (app_ra, app_dec)



def day_frac_to_hms(day_frac):
    '''Convert a fractional day into an (hr, min, sec) tuple.'''

    # Discard any surplus days
    (day_frac, junk) = modf(day_frac)

    # Convert the fractional day into hours
    hrs_frac = day_frac * 24

    # Extract the fractional and integer parts of the hours
    (mins_frac, hrs) = modf(hrs_frac)

    # Turn the remainder into fractional minutes
    mins_frac = mins_frac * 60

    # Extract the fractional and integer parts of the minutes
    (secs, mins) = modf(mins_frac)

    # Turn the remainder into seconds
    secs = secs * 60

    return (hrs, mins, secs)



def refine_day_fraction(app_sidereal_time, m_0, m_1, m_2, tdb, target, site,
                        std_altitude):
    '''Take an approximate value for transit, rise and set, and interpolate
       across the date boundary to obtain corrections to the values. The
       refined times are accurate to the nearest minute.
    '''


    # Find the sidereal time at Greenwich (in degrees)
    sidereal_time_transit = sidereal_time_at_greenwich(app_sidereal_time, m_0)
    sidereal_time_rise    = sidereal_time_at_greenwich(app_sidereal_time, m_1)
    sidereal_time_set     = sidereal_time_at_greenwich(app_sidereal_time, m_2)


    _log.debug('gwich sidereal_time (rise): %s'    % sidereal_time_rise)
    _log.debug('gwich sidereal_time (transit): %s' % sidereal_time_transit)
    _log.debug('gwich sidereal_time (set): %s'     % sidereal_time_set)


    # Calculate 'n' as per book instructions
    n_0 = calc_tabular_interval(m_0, tdb)
    n_1 = calc_tabular_interval(m_1, tdb)
    n_2 = calc_tabular_interval(m_2, tdb)


    # Calculate RA/Dec over 3 days for interpolation
    if ('planet' in target ):
        (alpha_1, delta_1) = apparent_planet_pos(target['planet'], tdb-1, site)
        (alpha_2, delta_2) = apparent_planet_pos(target['planet'], tdb, site)
        (alpha_3, delta_3) = apparent_planet_pos(target['planet'], tdb+1, site)
    else:
        (alpha_1, delta_1) = mean_to_apparent(target, tdb-1)
        (alpha_2, delta_2) = mean_to_apparent(target, tdb)
        (alpha_3, delta_3) = mean_to_apparent(target, tdb+1)

    _log.debug('alpha_1 (yesterday): %s' % alpha_1.in_degrees())
    _log.debug('alpha_2 (today): %s' % alpha_2.in_degrees())
    _log.debug('alpha_3 (tomorrow): %s' % alpha_3.in_degrees())


    # TODO: Handle wrapping across 24 hr boundary



    # Construct the first and second differences
    a = alpha_2.in_degrees() - alpha_1.in_degrees()
    b = alpha_3.in_degrees() - alpha_2.in_degrees()
    c = b - a

    _log.debug('a: %s' % a)
    _log.debug('b: %s' % b)
    _log.debug('c: %s' % c)

    interp_alpha_2_transit = interpolate(alpha_2.in_degrees(), n_0, a, b, c)
    interp_alpha_2_rise    = interpolate(alpha_2.in_degrees(), n_1, a, b, c)
    interp_alpha_2_set     = interpolate(alpha_2.in_degrees(), n_2, a, b, c)


    _log.debug('Original alpha_2: %s' % alpha_2.in_degrees())
    _log.debug('interpolated_alpha_2_rise: %s' % interp_alpha_2_rise)
    _log.debug('interpolated_alpha_2_transit: %s' % interp_alpha_2_transit)
    _log.debug('interpolated_alpha_2_set: %s' % interp_alpha_2_set)

    # Construct the first and second differences
    a = delta_2.in_degrees() - delta_1.in_degrees()
    b = delta_3.in_degrees() - delta_2.in_degrees()
    c = b - a

    _log.debug('a: %s' % a)
    _log.debug('b: %s' % b)
    _log.debug('c: %s' % c)


    interp_delta_2_transit = interpolate(delta_2.in_degrees(), n_0, a, b, c)
    interp_delta_2_rise    = interpolate(delta_2.in_degrees(), n_1, a, b, c)
    interp_delta_2_set     = interpolate(delta_2.in_degrees(), n_2, a, b, c)

    _log.debug('Original delta_2: %s' % delta_2.in_degrees())
    _log.debug('interpolated_delta_2_rise: %s' % interp_delta_2_rise)
    _log.debug('interpolated_delta_2_transit: %s' % interp_delta_2_transit)
    _log.debug('interpolated_delta_2_set: %s' % interp_delta_2_set)

    # Calculate the local hour angle (in degrees)
    local_hour_angle_transit = (sidereal_time_transit 
                                 + site['longitude']
                                 - interp_alpha_2_transit)

    local_hour_angle_rise = (sidereal_time_rise + site['longitude']
                                 - interp_alpha_2_rise)

    local_hour_angle_set = (sidereal_time_set + site['longitude']
                                 - interp_alpha_2_set)


    # TODO: Check local hour angle lies between -180 and +180
    _log.debug('local_hour_angle_rise: %s' % local_hour_angle_rise)
    _log.debug('local_hour_angle_transit: %s' % local_hour_angle_transit)
    _log.debug('local_hour_angle_set: %s' % local_hour_angle_set)

    refined_m_0 = correct_transit(m_0, local_hour_angle_transit)

    refined_m_1 = correct_rise_set(m_1, site['latitude'], interp_delta_2_rise,
                                   local_hour_angle_rise, std_altitude)

    refined_m_2 = correct_rise_set(m_2, site['latitude'], interp_delta_2_set,
                                   local_hour_angle_set, std_altitude)


    return (refined_m_0, refined_m_1, refined_m_2)



def sidereal_time_at_greenwich(app_sidereal_time, m):
    return app_sidereal_time.in_degrees() + (360.985647 * m)



def calc_tabular_interval(m, tdb):
    '''Find n, the tabular interval (Ast.Alg., p.99).'''
    return m + (sla.sla_dtt(tdb) / 86400)



def interpolate(y_2, n, a, b, c):
    '''See Eqn 3.3, p.25 Astro. Alg.'''
    y = y_2 + (n/2) * (a + b + n*c)

    return y



def correct_transit(m, local_hour_angle):
    '''Given an hour angle corrected for time of day, calculate and apply the
       correction to the transit time, as a fraction of a day.
       (Astro. Alg., p99)
                   H
       delta_m = -----
                  360
    '''

    # Local hour angle must be normalised to between -180 and +180
    if ( local_hour_angle < -180 ):
        local_hour_angle += 360
    elif ( local_hour_angle > 180 ):
        local_hour_angle -= 360


    # Find the correction to m_0 (transit)
    delta_m = -1 * (local_hour_angle / 360)

    _log.debug('delta_m (transit): %s' % delta_m)

    # Return the corrected value
    return m + delta_m



def correct_rise_set(m, latitude, dec, local_hour_angle, std_altitude):
    '''Given an hour angle and declination corrected for time of day, calculate
       and apply the correction to the rise or set time, as a fraction of a day.
       (Astro. Alg., p99)
                                h - h_0
       delta_m =     ------------------------------
                     360 cos(delta) cos(phi) sin(H)
    '''

    altitude = calculate_altitude(latitude, dec, local_hour_angle)

    _log.debug('altitude (rise/set): %s' % altitude.in_degrees())



    delta_m = ((altitude.in_degrees() - std_altitude.in_degrees())
              / ( 360 * cos(radians(dec)) * cos(radians(latitude))
                      * sin(radians(local_hour_angle)) ))

    _log.debug('delta_m (rise/set): %s' % delta_m)

    return m + delta_m



def calculate_altitude(latitude, dec, local_hour_angle):
    '''Find the altitude of the target.
       Eqn 12.6 Ast.Alg.
       sin h = sin(phi) sin(delta) + cos(phi)cos(delta)cos(H)
    '''

    altitude = asin( (sin(radians(latitude)) * sin(radians(dec)))
                 + (cos(radians(latitude)) * cos(radians(dec))
                    * cos(radians(local_hour_angle)))
               )

    return Angle(radians=altitude)


def get_dark_intervals(site, start_date, end_date):

    target = 'sun'

    # Find rise/set/transit for each day
    intervals = []
    current_date = start_date
    while current_date < end_date:
        # find_when_target_is_up()
        day_intervals = find_when_target_is_down(target, site, current_date)

        # Add today's intervals to the accumulating list of intervals
        intervals.extend(day_intervals)

        # Move on to tomorrow
        current_date += ONE_DAY

    # Collapse adjacent intervals into continuous larger intervals
    intervals = coalesce_adjacent_intervals(intervals)

    return intervals


def find_when_target_is_down(target, site, dt):

    dt = dt.replace(hour=0, minute=0, second=0, microsecond=0)

    up_intervals = find_when_target_is_up(target, site, current_date)

    # Treat the edges specially

    # If the first value has time 00:00:00, then the target starts up
    if up_intervals[0][0].time() == MIDNIGHT:
        pass

    # Otherwise the target starts down - there's one extra interval at start
    else:
        down_start = dt
        down_end   = up_intervals[0][0]

        down_intervals.append(down_start, down_end)


    # Proceed through the intervals, extracting the gaps
    for i in range(len(up_intervals) - 1):
        down_start = up_intervals[i][1]
        down_end   = up_intervals[i+1][0]

        down_intervals.append(down_start, down_end)


    # If the target sets before the end of the day, grab that as an
    # extra down interval
    if up_intervals[-1][1].time() != MIDNIGHT:
        down_start = up_intervals[-1][1]
        down_end   = dt + ONE_DAY

        down_intervals.append(down_start, down_end)

    return down_intervals


def find_when_target_is_up(target, site, dt):

    # Remove any time component of the provided datetime object
    dt = dt.replace(hour=0, minute=0, second=0, microsecond=0)

    # Get the rise/set/transit times for this day
    if target == 'sun':
        (transit, rise, set) = calc_sunrise_set(site, dt, 'sunrise')
    else:
        (transit, rise, set) = calc_rise_set(target, site, dt)

    intervals = []

    # Case 1: Overlapping start of day boundary
    # Target rose yesterday, and sets today. Rises again later today.
    # Rise  0hr    Transit      Set    Rise  24hr
    if (rise > transit) and (set > transit):
        # Store the first interval - start of day until target set
        absolute_set = dt.replace(hour=set[0], minute=set[1], second=set[2])
        intervals.append((dt, absolute_set))

        # Store the second interval - target rise until end of day
        absolute_rise = dt.replace(hour=rise[0], minute=rise[1], second=rise[2])
        intervals.append((absolute_rise, dt + ONE_DAY))


    # Case 2: Rise, set and transit all fall within the day, in order
    # Target rises today, transits, and sets before the day ends
    # 0hr    Rise      Transit      Set    24hr
    elif (rise < transit) and (set > transit):
        # Only one interval - rise until target set
        absolute_rise = dt.replace(hour=rise[0], minute=rise[1], second=rise[2])
        absolute_set  = dt.replace(hour=set[0], minute=set[1], second=set[2])
        intervals.append(absolute_rise, absolute_set)


    # Case 3: Overlapping end of day boundary
    # Target rose yesterday, and sets today. Rises again later today.
    # Rise  Transit 0hr    Set    Rise  Transit 24hr
    elif (rise < transit) and (set < transit):
        # Same code as case 1!
        # Store the first interval - start of day until target set
        absolute_set = dt.replace(hour=set[0], minute=set[1], second=set[2])
        intervals.append((dt, absolute_set))

        # Store the second interval - target rise until end of day
        absolute_rise = dt.replace(hour=rise[0], minute=rise[1], second=rise[2])
        intervals.append((absolute_rise, dt + ONE_DAY))


    return intervals


def coalesce_adjacent_intervals(intervals):
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



class InvalidDateError(Exception):
    '''Raised when an invalid date is encountered.'''

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value



class RiseSetError(Exception):
    '''Raised when a target either never rises or never sets.'''

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value
