#!/usr/bin/env python

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
from math import sin, cos, asin, acos, radians, modf
from datetime import timedelta

# Third party imports
import slalib as sla

# Internal imports
from rise_set.angle import Angle
from rise_set.sky_coordinates import RightAscension, Declination
from rise_set.rates import ProperMotion

# Import logging modules
import logging


_log = logging.getLogger('rise_set.astrometry')


class Star(object):
    # TODO: This is a crap name - change it

    def __init__(self, latitude, target, horizon=0.0):
        '''
             A star is circumpolar to a Northern hemisphere observer if the
             latitude plus the declination is greater than 90 degrees.

             A star is circumpolar to a Southern hemisphere observer if the
             latitude plus the declination is less than 90 degrees.

             If no horizon is provided, we approximate the effects of refraction.
             For positive horizons, the effects of refraction rapidly become
             negligible (and in any case require air temperature and pressure to
             calculate accurately).
        '''

        horizon_angle = Angle(degrees=horizon)

        # TODO: This function should be part of this object
        self.horizon  = apply_refraction_to_horizon(horizon_angle)
        self.target   = target
        self.latitude = latitude


    def is_always_up(self, date):

        lat, hor, dec = self._get_lat_hor_and_app_dec_in_degrees(date)

        # If the observer is in the Northern hemisphere...
        if lat > 0.0:
            if lat + dec - hor > 90.0:
                return True
            return False

        # Otherwise, if the observer is in the Southern hemisphere...
        elif lat < 0.0:
            if lat + dec + hor < -90.0:
                return True
            return False

        # We're on the exact equator - no circumpolar stars here!
        return False


    def is_always_down(self, date):

        lat, hor, dec = self._get_lat_hor_and_app_dec_in_degrees(date)

        # If the observer is in the Northern hemisphere...
        if lat > 0.0:
            if dec - lat - hor < -90.0:
                return True
            return False

        # Otherwise, if the observer is in the Southern hemisphere...
        elif lat < 0.0:
            if dec - lat + hor > 90:
                return True
            return False

        # We're on the exact equator - no circumpolar stars here!
        return False


    def _get_lat_hor_and_app_dec_in_degrees(self, date):
        tdb        = date_to_tdb(date)
        _, app_dec = mean_to_apparent(self.target, tdb)

        # TODO: Fix angles so they can be added together natively
        lat = self.latitude.in_degrees()
        dec = app_dec.in_degrees()
        hor = self.horizon.in_degrees()

        return (lat, hor, dec)

    def __repr__(self):
        return str(self.__dict__)



def gregorian_to_ut_mjd(date):
    '''Convert Gregorian calendar date to UTC MJD.'''

    # Do the date part
    caldj_error = {
                     0 : 'OK',
                     1 : 'bad year (MJD not computed)',
                     2 : 'bad month (MJD not computed)',
                     3 : 'bad day (MJD not computed)',
                  }
    caldj_status = 0

    mjd, caldj_status = sla.sla_caldj(date.year, date.month, date.day)

    if caldj_status != 0:
        raise InvalidDateTimeError('Error:' + caldj_error[caldj_status])


    # Do the time part
    dtf2d_error = {
                    0 : 'OK',
                    1 : 'IHOUR outside range 0-23',
                    2 : 'IMIN outside range 0-59',
                    3 : 'SEC outside range 0-59.999',
                  }
    dtf2d_status = 0

    days, dtf2d_status = sla.sla_dtf2d(date.hour, date.minute, date.second + (date.microsecond / 1e6))

    if dtf2d_status != 0:
        raise InvalidDateTimeError('Error:' + dtf2d_error[dtf2d_status])


    mjd += days

    return mjd



def ut_mjd_to_gmst(mjd):
    '''Convert UTC MJD to Greenwich mean sidereal time (an Angle).
       Note: We are assuming that UTC == UT1 here, which is what
       sla_gmst really expects. UT1 can't be easily determined, and
       we can only be out by less than 0.9s for as long as leap seconds
       persist.'''

    gmst_in_radians = sla.sla_gmst(mjd)

    return Angle(radians=gmst_in_radians)



def calc_apparent_sidereal_time(date):
    '''Return the apparent sidereal time (an Angle) for a given date.'''

    # Convert Gregorian to UT MJD
    mjd_utc = gregorian_to_ut_mjd(date)
    mjd_tdb = date_to_tdb(date)

    # Convert UT MJD to Greenwich mean sidereal time (in radians)
    gmst = ut_mjd_to_gmst(mjd_utc)

    # Find the apparent sidereal time (GAST)
    # GAST = GMST + Eqn. of equinoxes
    gast_in_rads = gmst.in_radians() + sla.sla_eqeqx(mjd_tdb)

    return Angle(radians=gast_in_rads)


def calc_local_sidereal_time(longitude, date):
    app_sidereal_time = calc_apparent_sidereal_time(date)

    local_sidereal_time = Angle(degrees=app_sidereal_time.in_degrees() + longitude.in_degrees())

    return local_sidereal_time


def calc_local_hour_angle(ra_app, longitude, date):
    ''' Ast. Algorithms p.92 (with reversed longitude convention)
        H = theta_0 + L - alpha

        where:
            theta_0 = GAST (Greenwich apparent sidereal time)
            L       = Site longitude (east +ve)
            alpha   = Apparent Right Ascension
     '''
    app_sidereal_time = calc_apparent_sidereal_time(date)

    local_hour_angle = calc_local_sidereal_time(longitude, date).in_degrees() - ra_app.in_degrees()

    if local_hour_angle < -180.0:
        local_hour_angle += 360.0
    elif local_hour_angle > 180.0:
        local_hour_angle -= 360.0

    return Angle(degrees=local_hour_angle)



def make_ra_dec_target(ra, dec, ra_proper_motion=None, dec_proper_motion=None, parallax=None,
                       rad_vel=None, epoch=None):

    target = {
                'ra'                : ra,
                'dec'               : dec,
                'ra_proper_motion'  : ra_proper_motion or ProperMotion(RightAscension(0), time='year'),
                'dec_proper_motion' : dec_proper_motion or ProperMotion(Declination(0), time='year'),
                'parallax'          : parallax or 0.0,
                'rad_vel'           : rad_vel or 0.0,
                'epoch'             : epoch or 2000,
             }

    return target


def make_minor_planet_target(target_type, epoch, inclination, long_node, arg_perihelion,
                              semi_axis, eccentricity, mean_anomaly):

    target = {
               'type'           : target_type,
               'epoch'          : epoch,
               'inclination'    : Angle(degrees=inclination),
               'long_node'      : Angle(degrees=long_node),
               'arg_perihelion' : Angle(degrees=arg_perihelion),
               'semi_axis'      : semi_axis,
               'eccentricity'   : eccentricity,
               'mean_anomaly'   : Angle(degrees=mean_anomaly),
             }

    return target


def make_comet_target(target_type, epoch, epochofperih, inclination, long_node, arg_perihelion,
                              perihdist, eccentricity):

    target = {
               'type'           : target_type,
               'epoch'          : epoch,
               'epochofperih'   : epochofperih,
               'inclination'    : Angle(degrees=inclination),
               'long_node'      : Angle(degrees=long_node),
               'arg_perihelion' : Angle(degrees=arg_perihelion),
               'perihdist'      : perihdist,
               'eccentricity'   : eccentricity,
             }

    return target


def mean_to_apparent(target, tdb):
    '''Given a target and TDB, return an apparent (RA, Dec) tuple.
       Thin wrapper for SLA_MAP.
    '''

    # Complain if the minimum target fields aren't present
    if not target.get('ra'):
        raise IncompleteTargetError("Missing RA in target definition")

    if not target.get('dec'):
        raise IncompleteTargetError("Missing Declination in target definition")

    target = make_ra_dec_target(target['ra'], target['dec'],
                                ra_proper_motion=target.get('ra_proper_motion'),
                                dec_proper_motion=target.get('dec_proper_motion'),
                                parallax=target.get('parallax'),
                                rad_vel=target.get('rad_vel'),
                                epoch=target.get('epoch'))



    (ra_app_rads, dec_app_rads) = sla.sla_map(
                                  target['ra'].in_radians(),
                                  target['dec'].in_radians(),
                                  target['ra_proper_motion'].in_radians_per_year(),
                                  target['dec_proper_motion'].in_radians_per_year(),
                                  target['parallax'],
                                  target['rad_vel'],
                                  target['epoch'],
                                  tdb)

    ra_apparent  = Angle(radians=ra_app_rads)
    dec_apparent = Angle(radians=dec_app_rads)

    return (ra_apparent, dec_apparent)



def calc_rise_set_hour_angle(latitude, dec_apparent, std_altitude):
    '''Find the hour angle H_0 corresponding to the time of rise or set of
       a  celestial body, according to

                  sin(h_0) - [ sin(phi) sin(delta_2) ]
       cos H_0 =  -----------------------------------
                         cos(phi) cos(delta_2)
    '''



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



def calc_transit_day_fraction(ra, longitude, app_sidereal_time):
    '''Find the time, expressed as a fraction of a day, when the transit occurs.

       Note: Sign of longitude is reversed relative to Astro Alg., p.98, because
       SLALIB uses East +ve in all functions except SLA_OBS, but Astro. Alg.
       uses West +ve. <Sigh>

              alpha_2 - L - theta_0
       m_0 = ------------------------
                      360
    '''

    m_0 = (ra.in_degrees() - longitude.in_degrees()
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



def apply_refraction_to_horizon(horizon):
    '''If using a horizon of zero, adds an average refraction term to improve
       the effective horizon. If no horizon is provided, a horizon of zero is
       assumed.

       Non-zero horizons fall through and are not modified.
    '''

    # For non-zero horizons, we shall ignore refraction effects
    std_alt_of_stars = Angle(degrees=0.0)

    if not horizon:
        # Default to the Earth's horizon
        horizon = Angle(degrees=0.0)

    # Approximate the effect of refraction if we are using the true horizon
    if horizon.in_degrees() == 0.0:
        std_alt_of_stars = Angle(degrees=-0.5667)

    effective_horizon_in_deg = std_alt_of_stars.in_degrees() + horizon.in_degrees()
    effective_horizon = Angle(degrees=effective_horizon_in_deg)

    return effective_horizon



def date_to_tdb(date):
    '''Converts a given UTC datetime (<date>; e.g. datetime(2015, 4, 8, 0, 0))
    to a *TT* Modified Julian Date (MJD; e.g. 57120.000777592591).
    TT is ahead of UTC by a fixed +32.184 seconds from TAI and a 
    variable (currently (2015-04-14) 35 seconds) offset, which includes 
    leapseconds, from UTC to TAI.
    Notes: 
    1) This code relies on SLALIB's sla_dat routine to be updated and recompiled
    when new leapseconds are announced,
    2) Despite the (terrible) name, this routine *does not* perform the 
    relativstic clock corrections (e.g. sla_rcc) to produce TDB - max error
    is ~2ms'''
    ut_mjd = gregorian_to_ut_mjd(date)
    tdb = ut_mjd + (sla.sla_dtt(ut_mjd)/86400)

    return tdb



def calc_rise_set(target, site, date, horizon=None):
    '''Return a tuple (transit, rise, set) of timedelta objects, describing the
       time offset for each event from the start of the provided date.
    '''

    # Remove any time component of the provided datetime object
    date = date.replace(hour=0, minute=0, second=0, microsecond=0)

    effective_horizon = apply_refraction_to_horizon(horizon)
    tdb = date_to_tdb(date)


    app_ra, app_dec   = mean_to_apparent(target, tdb)
    app_sidereal_time = calc_apparent_sidereal_time(date)
    (hour_angle, msg) = calc_rise_set_hour_angle(site['latitude'], app_dec,
                                                 effective_horizon)

    if ( not hour_angle ):
        msg += " (ra=%s, dec=%s, lat=%s)" % ( target.get('ra').in_sexegesimal(),
                                              target.get('dec').in_sexegesimal(),
                                              site['latitude'].in_degrees())
        raise RiseSetError(msg)

    m_0 = calc_transit_day_fraction(app_ra, site['longitude'], app_sidereal_time)
    m_1 = calc_rising_day_fraction(m_0, hour_angle)
    m_2 = calc_setting_day_fraction(m_0, hour_angle)

    transits = day_frac_to_hms(m_0)
    rises    = day_frac_to_hms(m_1)
    sets     = day_frac_to_hms(m_2)

    _log.info('Rise time - unrefined (h, m, s): %s', rises)
    _log.info('Transit time - unrefined (h, m, s): %s', transits)
    _log.info('Set time - unrefined (h, m, s): %s', sets)

    (m_0, m_1, m_2) = refine_day_fraction(app_sidereal_time, m_0, m_1, m_2, tdb,
                                         target, site, effective_horizon)

    transits = timedelta(days=m_0)
    rises    = timedelta(days=m_1)
    sets     = timedelta(days=m_2)


    return (transits, rises, sets)



def calc_sunrise_set(site, date, twilight):
    '''Return a tuple (transit, rise, set) of timedelta objects, describing the
       time offset for each event from the start of the provided date.
    '''

    # Remove any time component of the provided datetime object
    date = date.replace(hour=0, minute=0, second=0, microsecond=0)

    ut_mjd = gregorian_to_ut_mjd(date)
    tdb = ut_mjd + (sla.sla_dtt(ut_mjd)/86400)

    target = dict(planet = 'sun')
    (app_ra, app_dec) = apparent_planet_pos(target['planet'], tdb, site)

    _log.info("RA, Dec (apparent, degrees) for %s: (%s, %s)",
              target['planet'], app_ra.in_degrees(), app_dec.in_degrees())

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

    transits = day_frac_to_hms(m_0)
    rises    = day_frac_to_hms(m_1)
    sets     = day_frac_to_hms(m_2)

    _log.info('Rise time - unrefined (h, m, s): %s', rises)
    _log.info('Transit time - unrefined (h, m, s): %s', transits)
    _log.info('Set time - unrefined (h, m, s): %s', sets)

    (m_0, m_1, m_2) = refine_day_fraction(app_sidereal_time, m_0, m_1, m_2, tdb,
                                          target, site, sun_std_alt[twilight])

    transits = timedelta(days=m_0)
    rises    = timedelta(days=m_1)
    sets     = timedelta(days=m_2)

    return (transits, rises, sets)



def apparent_planet_pos(planet_name, tdb, site):
    '''Return the topocentric apparent position (ra, dec) tuple of a planet at
       a particular time, from a particular site.
       Thin wrapper for SLA_RDPLAN.
    '''

    latitude  = site['latitude']
    longitude = site['longitude']

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


    (app_ra_rads, app_dec_rads, _) = sla.sla_rdplan(tdb,
                                                    planet[planet_name],
                                                    longitude.in_radians(),
                                                    latitude.in_radians())

    app_ra  = Angle(radians=app_ra_rads)
    app_dec = Angle(radians=app_dec_rads)

    return (app_ra, app_dec)



def day_frac_to_hms(day_frac):
    '''Convert a fractional day into an (hr, min, sec) tuple.'''

    # Discard any surplus days
    day_frac, _ = modf(day_frac)

    # Convert the fractional day into hours
    hrs_frac = day_frac * 24

    # Extract the fractional and integer parts of the hours
    mins_frac, hrs = modf(hrs_frac)

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


    _log.debug('gwich sidereal_time (rise): %s', sidereal_time_rise)
    _log.debug('gwich sidereal_time (transit): %s', sidereal_time_transit)
    _log.debug('gwich sidereal_time (set): %s', sidereal_time_set)

    _log.debug('m_0: %s', m_0)
    _log.debug('m_1: %s', m_1)
    _log.debug('m_2: %s', m_2)

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

    _log.debug('alpha_1 (yesterday): %s', alpha_1.in_degrees())
    _log.debug('alpha_2 (today): %s', alpha_2.in_degrees())
    _log.debug('alpha_3 (tomorrow): %s', alpha_3.in_degrees())


    # Handle wrapping across 24 hr boundary
    # Why do we use 350 degrees? That's a very good question.
    # We need to determine when the alpha has wrapped, but we can't do
    # the obvious test, because of valid cases.
    # Instead, we use the fact that the differences between the alphas
    # need to be large, otherwise the object is moving extremely fast,
    # to determine the wrapping scenario.
    if alpha_2.in_degrees() < alpha_1.in_degrees():
        if alpha_2.in_degrees() - alpha_1.in_degrees() < -350:
            norm_alpha2 = alpha_2.in_degrees() + 360
            alpha_2 = Angle(degrees=norm_alpha2)

    if alpha_3.in_degrees() < alpha_1.in_degrees():
        if alpha_3.in_degrees() - alpha_1.in_degrees() < -350:
            norm_alpha3 = alpha_3.in_degrees() + 360
            alpha_3 = Angle(degrees=norm_alpha3)

    _log.debug('alpha_1 normalised (yesterday): %s', alpha_1.in_degrees())
    _log.debug('alpha_2 normalised (today): %s', alpha_2.in_degrees())
    _log.debug('alpha_3 normalised (tomorrow): %s', alpha_3.in_degrees())

    # Construct the first and second differences
    a = alpha_2.in_degrees() - alpha_1.in_degrees()
    b = alpha_3.in_degrees() - alpha_2.in_degrees()
    c = b - a


    interp_alpha_2_transit = interpolate(alpha_2.in_degrees(), n_0, a, b, c)
    interp_alpha_2_rise    = interpolate(alpha_2.in_degrees(), n_1, a, b, c)
    interp_alpha_2_set     = interpolate(alpha_2.in_degrees(), n_2, a, b, c)


    # Construct the first and second differences
    a = delta_2.in_degrees() - delta_1.in_degrees()
    b = delta_3.in_degrees() - delta_2.in_degrees()
    c = b - a


    interp_delta_2_rise    = interpolate(delta_2.in_degrees(), n_1, a, b, c)
    interp_delta_2_set     = interpolate(delta_2.in_degrees(), n_2, a, b, c)


    # Calculate the local hour angle (in degrees)
    local_hour_angle_transit = (sidereal_time_transit
                                 + site['longitude'].in_degrees()
                                 - interp_alpha_2_transit)

    local_hour_angle_rise = (sidereal_time_rise + site['longitude'].in_degrees()
                                 - interp_alpha_2_rise)

    local_hour_angle_set = (sidereal_time_set + site['longitude'].in_degrees()
                                 - interp_alpha_2_set)


    # TODO: Check local hour angle lies between -180 and +180
    _log.debug('local_hour_angle_rise: %s',    local_hour_angle_rise)
    _log.debug('local_hour_angle_transit: %s', local_hour_angle_transit)
    _log.debug('local_hour_angle_set: %s',     local_hour_angle_set)

    refined_m_0 = correct_transit(m_0, local_hour_angle_transit)

    refined_m_1 = correct_rise_set(m_1, site['latitude'].in_degrees(),
                                   interp_delta_2_rise, local_hour_angle_rise,
                                   std_altitude)

    refined_m_2 = correct_rise_set(m_2, site['latitude'].in_degrees(),
                                   interp_delta_2_set, local_hour_angle_set,
                                   std_altitude)

    refined_m_0 = normalise_day(refined_m_0)
    refined_m_1 = normalise_day(refined_m_1)
    refined_m_2 = normalise_day(refined_m_2)

    return (refined_m_0, refined_m_1, refined_m_2)



def sidereal_time_at_greenwich(app_sidereal_time, m):
    sidereal_m =  app_sidereal_time.in_degrees() + (360.985647 * m)

    while sidereal_m > 360.0:
        sidereal_m -= 360.0

    return sidereal_m



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

    _log.debug('delta_m (transit): %s', delta_m)

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

    _log.debug('altitude (rise/set): %s', altitude.in_degrees())



    delta_m = ((altitude.in_degrees() - std_altitude.in_degrees())
              / ( 360 * cos(radians(dec)) * cos(radians(latitude))
                      * sin(radians(local_hour_angle)) ))

    _log.debug('delta_m (rise/set): %s', delta_m)

    return m + delta_m



def calculate_altitude(latitude, dec, local_hour_angle):
    '''Find the altitude of the target.
       Eqn 13.6 Ast.Alg.
       sin h = sin(phi) sin(delta) + cos(phi)cos(delta)cos(H)
    '''

    altitude = asin( (sin(radians(latitude)) * sin(radians(dec)))
                 + (cos(radians(latitude)) * cos(radians(dec))
                    * cos(radians(local_hour_angle)))
               )

    return Angle(radians=altitude)


class InvalidDateTimeError(Exception):
    '''Raised when an invalid date is encountered.'''
    pass


class RiseSetError(Exception):
    '''Raised when a target either never rises or never sets.'''
    pass


class IncompleteTargetError(Exception):
    '''Raised when a target is missing a key value (RA, Dec).'''
    pass
