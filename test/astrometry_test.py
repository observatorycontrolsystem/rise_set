#!/usr/bin/python
from __future__ import division

from nose.tools import assert_equal, assert_almost_equal, assert_less, raises, assert_greater
from datetime import datetime, timedelta

#Import the module to test
from rise_set.astrometry import (InvalidDateTimeError, IncompleteTargetError, RiseSetError,
                                 RightAscension, Declination, Star, ProperMotion,
                                 gregorian_to_ut_mjd, mean_to_apparent,
                                 date_to_tdb, calc_sunrise_set, calculate_airmass_at_times,
                                 calc_rise_set, calc_setting_day_fraction,
                                 calc_rise_set_hour_angle, calc_rising_day_fraction,
                                 calc_transit_day_fraction, day_frac_to_hms, calc_planet_rise_set,
                                 calc_local_hour_angle, make_ra_dec_target, apparent_planet_pos)

from rise_set.angle import Angle

class YiannisIsTryingToBreakMyDateCalculator(object):
    def __init__(self):
        self.day   = 20
        self.month = 30
        self.year  = 1988


class TestLeapSeconds(object):
    '''Verification that leap seconds are up-to-date within SLALIB.'''

    def test_pre2009_jan1_leapsecond(self):
        date     = datetime(2008, 12, 31, 12, 0, 0)
        received = date_to_tdb(date)
        expected = 54831.5 + (65.184/86400)
        assert_equal(received, expected)



class TestAstrometry(object):
    '''Unit tests for the astrometry module.'''

    def setup(self):
        self.date      = datetime(year=1988, month=3, day=20)
        self.bad_month = YiannisIsTryingToBreakMyDateCalculator()
        self.mjd       = 47240.0


    def test_gregorian_to_ut_mjd(self):
        assert_equal(gregorian_to_ut_mjd(self.date), self.mjd)

    def test_gregorian_to_ut_mjd_microseconds(self):
        dt = datetime(2014, 6, 14)
        tdb = dt - timedelta(seconds=67.184)
        mjd_micros = 56821.99922240741
        assert_equal(gregorian_to_ut_mjd(tdb), mjd_micros)


    def test_date_to_tdb(self):
        date     = datetime(2013, 11, 4)
        received = date_to_tdb(date)
        expected = 56600.0 + (67.184/86400)
        assert_equal(received, expected)


    @raises(InvalidDateTimeError)
    def test_gregorian_to_ut_mjd_bad_month(self):
        gregorian_to_ut_mjd(self.bad_month)

#    def test_ut_mjd_to_gmst_in_radians(self):
#        assert_equal(ut_mjd_to_gmst_in_radians(self.mjd), 4)


    @raises(IncompleteTargetError)
    def test_missing_ra_raises_exception(self):
        tdb = self.mjd
        target_missing_ra = dict(dec='18 26 27.3')
        mean_to_apparent(target_missing_ra, tdb)


    @raises(IncompleteTargetError)
    def test_missing_dec_raises_exception(self):
        tdb = self.mjd
        target_missing_dec = dict(ra='02 46 55.51')
        mean_to_apparent(target_missing_dec, tdb)


    def test_calc_local_hour_angle(self):
        ra_app        = Angle(degrees=30)
        elp_longitude = Angle(degrees=-104.015194444)
        date          = datetime(2013, 12, 10)

        hour_angle = calc_local_hour_angle(ra_app, elp_longitude, date)

        assert_almost_equal(hour_angle.in_degrees(), -55.128564469690645, places=13)


    def test_calc_local_hour_angle_normalises_negative_overrun(self):
        ra_app        = Angle(degrees=288.75)
        elp_longitude = Angle(degrees=-104.015194444)
        date          = datetime(2014, 07, 18, 7)

        hour_angle = calc_local_hour_angle(ra_app, elp_longitude, date)

        assert_almost_equal(hour_angle.in_degrees(), 8.2509, places=4)


    def test_calc_local_hour_angle_normalises_positive_overrun(self):
        ra_app        = Angle(degrees=108.75)
        coj_longitude = Angle(degrees=149.070593)
        date          = datetime(2014, 1, 12, 14)

        hour_angle = calc_local_hour_angle(ra_app, coj_longitude, date)

        assert_almost_equal(hour_angle.in_degrees(), 2.3088, places=4)


class TestMoon(object):
    def setup(self):
        self.site = {
                        'name': 'test',
                        'latitude': Angle(degrees=-30.0),
                        'longitude': Angle(degrees=0.0)
                    }
        # 5 arcsecond tolerance
        self.tolerance = 5.0 / 3600.0

        # rise/set/transit time tolerance in seconds
        self.time_tolerance = 7.0 * 60.0

    def test_apparent_position(self):
        dt_tdb = gregorian_to_ut_mjd(datetime(2012, 1, 3))

        assert_equal(dt_tdb, 55929.0)

        (apparent_ra, apparent_dec) = apparent_planet_pos("moon", dt_tdb, self.site)

        # values from JPL Horizons
        expected_ra = Angle(degrees=26.63848)
        expected_dec = Angle(degrees=15.61778)

        assert_less(abs(apparent_ra.in_degrees() - expected_ra.in_degrees()), self.tolerance)
        assert_less(abs(apparent_dec.in_degrees() - expected_dec.in_degrees()), self.tolerance)

    def test_rise_set_12_31(self):
        dt_utc = datetime(2011, 12, 31)
        # parallax and semidiameter from astronomers almanac
        horizontal_parallax = Angle(degrees=0.9150694)
        semidiameter = Angle(degrees=0.2492556)
        # this equation is from astronomers almanac
        h_0 = Angle(degrees= 0.5666667 + semidiameter.in_degrees() - horizontal_parallax.in_degrees())
        # this equation is from the astronomical algorithms book
        h_0_2 = Angle(degrees= 0.7275*horizontal_parallax.in_degrees() - 0.5666667)

        (transits, rises, sets) = calc_planet_rise_set(self.site, dt_utc, h_0_2, 'moon')

        # values from JPL Horizons
        # There is no expected set time for the january 2nd, 2012 - it spans the date boundary
        expected_rise = timedelta(hours=11, minutes=30)
        expected_transit = timedelta(hours=17, minutes=27)
        expected_set = timedelta(hours=23, minutes=22)

        assert_less(abs(sets.total_seconds() - expected_set.total_seconds()), self.time_tolerance)
        assert_less(abs(rises.total_seconds() - expected_rise.total_seconds()), self.time_tolerance)
        assert_less(abs(transits.total_seconds() - expected_transit.total_seconds()), self.time_tolerance)

    def test_rise_set_1_1(self):
        dt_utc = datetime(2012, 1, 1)
        # parallax and semidiameter from astronomers almanac
        horizontal_parallax = Angle(degrees=0.9082611)
        semidiameter = Angle(degrees=0.2474)
        # this equation is from astronomers almanac
        h_0 = Angle(degrees= 0.5666667 + semidiameter.in_degrees() - horizontal_parallax.in_degrees())
        # this equation is from the astronomical algorithms book
        h_0_2 = Angle(degrees= 0.7275*horizontal_parallax.in_degrees() - 0.5666667)

        (transits, rises, sets) = calc_planet_rise_set(self.site, dt_utc, h_0_2, 'moon')

        # values from JPL Horizons
        # There is no expected set time for the january 2nd, 2012 - it spans the date boundary
        expected_rise = timedelta(hours=12, minutes=23)
        expected_transit = timedelta(hours=18, minutes=9)
        expected_set = timedelta(hours=23, minutes=54)

        assert_less(abs(sets.total_seconds() - expected_set.total_seconds()), self.time_tolerance)
        assert_less(abs(rises.total_seconds() - expected_rise.total_seconds()), self.time_tolerance)
        assert_less(abs(transits.total_seconds() - expected_transit.total_seconds()), self.time_tolerance)

    def test_rise_set_1_2(self):
        dt_utc = datetime(2012, 1, 2)
        # parallax and semidiameter from astronomers almanac
        horizontal_parallax = Angle(degrees=0.9043333)
        semidiameter = Angle(degrees=0.2463306)
        # this equation is from astronomers almanac
        h_0 = Angle(degrees= 0.5666667 + semidiameter.in_degrees() - horizontal_parallax.in_degrees())
        # this equation is from the astronomical algorithms book
        h_0_2 = Angle(degrees= 0.7275*horizontal_parallax.in_degrees() - 0.5666667)

        (transits, rises, sets) = calc_planet_rise_set(self.site, dt_utc, h_0_2, 'moon')

        # values from JPL Horizons
        # There is no expected set time for the january 2nd, 2012 - it spans the date boundary
        expected_rise = timedelta(hours=13, minutes=16)
        expected_transit = timedelta(hours=18, minutes=52)

        assert_less(abs(rises.total_seconds() - expected_rise.total_seconds()), self.time_tolerance)
        assert_less(abs(transits.total_seconds() - expected_transit.total_seconds()), self.time_tolerance)

    def test_rise_set_1_3(self):
        dt_utc = datetime(2012, 1, 3)
        horizontal_parallax = Angle(degrees=0.9033333)
        semidiameter = Angle(degrees=0.2460583)
        h_0 = Angle(degrees= 0.5666667 + semidiameter.in_degrees() - horizontal_parallax.in_degrees())
        h_0_2 = Angle(degrees= 0.7275*horizontal_parallax.in_degrees() - 0.5666667)

        (transits, rises, sets) = calc_planet_rise_set(self.site, dt_utc, h_0_2, 'moon')

        # values from JPL Horizons
        expected_set = timedelta(minutes=27)
        expected_rise = timedelta(hours=14, minutes=10)
        expected_transit = timedelta(hours=19, minutes=37)

        assert_less(abs(sets.total_seconds() - expected_set.total_seconds()), self.time_tolerance)
        assert_less(abs(rises.total_seconds() - expected_rise.total_seconds()), self.time_tolerance)
        assert_less(abs(transits.total_seconds() - expected_transit.total_seconds()), self.time_tolerance)

    def test_rise_set_1_4(self):
        dt_utc = datetime(2012, 1, 4)
        horizontal_parallax = Angle(degrees=0.905125)
        semidiameter = Angle(degrees=0.2465472)
        h_0 = Angle(degrees= 0.5666667 + semidiameter.in_degrees() - horizontal_parallax.in_degrees())
        h_0_2 = Angle(degrees= 0.7275*horizontal_parallax.in_degrees() - 0.5666667)

        (transits, rises, sets) = calc_planet_rise_set(self.site, dt_utc, h_0_2, 'moon')

        # values from JPL Horizons
        expected_set = timedelta(hours=1, minutes=4)
        expected_rise = timedelta(hours=15, minutes=4)
        expected_transit = timedelta(hours=20, minutes=24)

        assert_less(abs(sets.total_seconds() - expected_set.total_seconds()), self.time_tolerance)
        assert_less(abs(rises.total_seconds() - expected_rise.total_seconds()), self.time_tolerance)
        assert_less(abs(transits.total_seconds() - expected_transit.total_seconds()), self.time_tolerance)

    def test_rise_set_1_31(self):
        dt_utc = datetime(2012, 1, 31)
        horizontal_parallax = Angle(degrees=0.9039722)
        semidiameter = Angle(degrees=0.2462306)
        h_0 = Angle(degrees= 0.5666667 + semidiameter.in_degrees() - horizontal_parallax.in_degrees())
        h_0_2 = Angle(degrees= 0.7275*horizontal_parallax.in_degrees() - 0.5666667)

        (transits, rises, sets) = calc_planet_rise_set(self.site, dt_utc, h_0_2, 'moon')

        # values from JPL Horizons
        expected_set = timedelta(hours=23, minutes=40)
        expected_rise = timedelta(hours=12, minutes=53)
        expected_transit = timedelta(hours=18, minutes=17)

        assert_less(abs(sets.total_seconds() - expected_set.total_seconds()), self.time_tolerance)
        assert_less(abs(rises.total_seconds() - expected_rise.total_seconds()), self.time_tolerance)
        assert_less(abs(transits.total_seconds() - expected_transit.total_seconds()), self.time_tolerance)

    def test_rise_set_2_1(self):
        dt_utc = datetime(2012, 2, 1)
        horizontal_parallax = Angle(degrees=0.9062222)
        semidiameter = Angle(degrees=0.2468444)
        h_0 = Angle(degrees= 0.5666667 + semidiameter.in_degrees() - horizontal_parallax.in_degrees())
        h_0_2 = Angle(degrees= 0.7275*horizontal_parallax.in_degrees() - 0.5666667)

        (transits, rises, sets) = calc_planet_rise_set(self.site, dt_utc, h_0_2, 'moon')

        # values from JPL Horizons
        expected_set = timedelta(minutes=23)
        expected_rise = timedelta(hours=13, minutes=47)
        expected_transit = timedelta(hours=19, minutes=5)

        assert_less(abs(sets.total_seconds() - expected_set.total_seconds()), self.time_tolerance)
        assert_less(abs(rises.total_seconds() - expected_rise.total_seconds()), self.time_tolerance)
        assert_less(abs(transits.total_seconds() - expected_transit.total_seconds()), self.time_tolerance)

    def test_rise_set_2_2(self):
        dt_utc = datetime(2012, 2, 2)
        horizontal_parallax = Angle(degrees=0.9113389)
        semidiameter = Angle(degrees=0.2482389)
        h_0 = Angle(degrees= 0.5666667 + semidiameter.in_degrees() - horizontal_parallax.in_degrees())
        h_0_2 = Angle(degrees= 0.7275*horizontal_parallax.in_degrees() - 0.5666667)

        (transits, rises, sets) = calc_planet_rise_set(self.site, dt_utc, h_0_2, 'moon')

        # values from JPL Horizons
        expected_set = timedelta(minutes=23)
        expected_rise = timedelta(hours=14, minutes=40)
        expected_transit = timedelta(hours=19, minutes=55)

        assert_less(abs(sets.total_seconds() - expected_set.total_seconds()), self.time_tolerance)
        assert_less(abs(rises.total_seconds() - expected_rise.total_seconds()), self.time_tolerance)
        assert_less(abs(transits.total_seconds() - expected_transit.total_seconds()), self.time_tolerance)

    def test_rise_set_2_3(self):
        dt_utc = datetime(2012, 2, 3)
        horizontal_parallax = Angle(degrees=0.9190528)
        semidiameter = Angle(degrees=0.2503389)
        h_0 = Angle(degrees= 0.5666667 + semidiameter.in_degrees() - horizontal_parallax.in_degrees())
        h_0_2 = Angle(degrees= 0.7275*horizontal_parallax.in_degrees() - 0.5666667)

        (transits, rises, sets) = calc_planet_rise_set(self.site, dt_utc, h_0_2, 'moon')

        # values from JPL Horizons
        expected_set = timedelta(hours=1, minutes=11)
        expected_rise = timedelta(hours=15, minutes=31)
        expected_transit = timedelta(hours=20, minutes=47)

        assert_less(abs(sets.total_seconds() - expected_set.total_seconds()), self.time_tolerance)
        assert_less(abs(rises.total_seconds() - expected_rise.total_seconds()), self.time_tolerance)
        assert_less(abs(transits.total_seconds() - expected_transit.total_seconds()), self.time_tolerance)


class TestSunriseSunset(object):

    def setup(self):
        self.lsc =  {
                      'name'      : '1m0a.domb.lsc',
                      'latitude'  : Angle(degrees=-30.1673472222),
                      'longitude' : Angle(degrees=-70.8046722222),
                    }

        self.elp = {
                     'name'      : '1m0a.doma.elp',
                     'latitude'  : Angle(degrees=30.6801),
                     'longitude' : Angle(degrees=-104.015194444),
                   }



    def test_nautical_twilight_from_lsc_no_time(self):
        date = datetime(2013, 12, 10)
        twilight = 'nautical'

        expected = ('unknown', timedelta(hours=8, minutes=34),
                    timedelta(hours=0, minutes=38))

        received = calc_sunrise_set(self.lsc, date, twilight)

        # We only know these to the nearest 30s from USNO online calculator
        assert_less(abs(expected[1]-received[1]), timedelta(seconds=30))
        assert_less(abs(expected[2]-received[2]), timedelta(seconds=30))


    def test_nautical_twilight_from_lsc_with_time(self):
        date = datetime(2013, 12, 10, 12)
        twilight = 'nautical'

        expected = ('unknown', timedelta(hours=8, minutes=34),
                    timedelta(hours=0, minutes=38))

        received = calc_sunrise_set(self.lsc, date, twilight)

        # We only know these to the nearest 30s from USNO online calculator
        assert_less(abs(expected[1]-received[1]), timedelta(seconds=30))
        assert_less(abs(expected[2]-received[2]), timedelta(seconds=30))


    def test_sunrise_set_from_elp_no_time(self):
        date = datetime(2012, 5, 11)
        twilight = 'sunrise'

        expected = ('unknown', timedelta(hours=12, minutes=4),
                    timedelta(hours=1, minutes=41))

        received = calc_sunrise_set(self.elp, date, twilight)

        # We only know these to the nearest 30s from USNO online calculator
        assert_less(abs(expected[1]-received[1]), timedelta(seconds=30))
        assert_less(abs(expected[2]-received[2]), timedelta(seconds=30))


    def test_sunrise_set_from_elp_with_time(self):
        date = datetime(2012, 5, 11, 17, 30)
        twilight = 'sunrise'

        expected = ('unknown', timedelta(hours=12, minutes=4),
                    timedelta(hours=1, minutes=41))

        received = calc_sunrise_set(self.elp, date, twilight)

        # We only know these to the nearest 30s from USNO online calculator
        assert_less(abs(expected[1]-received[1]), timedelta(seconds=30))
        assert_less(abs(expected[2]-received[2]), timedelta(seconds=30))



class TestVenusRiseTransitSet(object):
    '''Implementation of Example 14.a from Astronomical Algorithms, p.99'''

    def setup(self):
        self.boston = {
           'latitude'  : Angle(degrees = 42.3333),
           'longitude' : Angle(degrees = -71.0833)
        }

        self.app_sidereal_time = Angle(degrees = '11 50 58.10', units = 'time')
        self.alpha_2 = Angle(degrees = '02 46 55.51', units = 'time')
        self.delta_2 = Angle(degrees = '18 26 27.3')

        self.std_alt = Angle(degrees=-0.5667)


    def test_calc_transit_day_fraction(self):
        m_0 = calc_transit_day_fraction(self.alpha_2, self.boston['longitude'],
                                        self.app_sidereal_time)
        assert_almost_equal(m_0, 0.81965, places=5)


        (hour_angle, _) = calc_rise_set_hour_angle(self.boston['latitude'],
                                                     self.delta_2,
                                                     self.std_alt)
        assert_almost_equal(hour_angle.in_degrees(), 108.5344, places=4)

        m_1 = calc_rising_day_fraction(m_0, hour_angle)
        assert_almost_equal(m_1, 0.51816, places=5)

        m_2 = calc_setting_day_fraction(m_0, hour_angle)
        assert_almost_equal(m_2, 0.12113, places=5)


    def test_day_frac_to_hms(self):
        (hrs, mins, secs) = day_frac_to_hms(0.531712963)

        assert_equal(hrs, 12)
        assert_equal(mins, 45)
        assert_almost_equal(secs, 40, places=5)



class TestDenebFromMaui(object):
    '''Integration test: rise/set/transit of a Northern star from a Northern
       observatory.'''

    # IMPORTANT NOTE: LATITUDES
    # IAU convention         = East is +ve
    # Config DB convention   = East is +ve,  e.g. Siding Spring = +149 04 14.13
    # Astro. Alg. convention = West is +ve
    # SLALIB convention (SLA_OBS) = West is +ve
    # SLALIB (rest of library) convention = East is +ve
    # Thanks. You bastards.


    def setup(self):

        # Target
        # Note: Aladin units are mas/yr...
        self.deneb = {
               'ra'                : RightAscension('20 41 25.91'),
               'dec'               : Declination('+45 16 49.22'),
               'ra_proper_motion'  : ProperMotion(RightAscension('00 00 00.00156')),
               'dec_proper_motion' : ProperMotion(Declination('00 00 00.00155')),
               'parallax'          : 0.00101,  # Units: arcsec
               'rad_vel'           : -4.5,  # Units: km/s (-ve approaches)
               'epoch'             : 2000,
              }

        # Site (East +ve longitude)
        self.maui = {
           'latitude'  : Angle(degrees = 20.7069444444),
           'longitude' : Angle(degrees = -156.258055556)
        }

        # Date
        # We provide an hour to prove that this is ignored by rise_set
        self.date = datetime(year=2010, month=10, day=25, hour=8)


    def test_rise_set(self):
        (transit, rise, sets) = calc_rise_set(self.deneb, self.maui, self.date)

        # Second accuracy isn't meaningful, so just check we would round up or
        # down as apppropriate.
        exp_secs = 30

        # Expected times taken from http://aa.usno.navy.mil/data/docs/mrst.php

        # Transit
        exp_t_hr  = 4
        exp_t_min = 52   # Not 53, because we expect to round up seconds

        the_date     = self.date.replace(hour=0, minute=0, second=0, microsecond=0)
        transit_time = the_date + transit

        assert_equal(transit_time.hour, exp_t_hr)
        assert_equal(transit_time.minute, exp_t_min)
        msg = '%r !> %r' % (transit_time.second, exp_secs)
        assert transit_time.second >= exp_secs, msg


        # Rise
        exp_r_hr  = 21
        exp_r_min = 16   # Not 17, because we expect to round up seconds

        rise_time = the_date + rise

        assert_equal(rise_time.hour, exp_r_hr)
        assert_equal(rise_time.minute, exp_r_min)
        assert rise_time.second >= exp_secs, '%r !> %r' % (rise_time.second, exp_secs)


        # Set
        exp_s_hr  = 12
        exp_s_min = 25   # Not 24, because we expect to round down seconds

        set_time = the_date + sets

        assert_equal(set_time.hour, exp_s_hr)
        assert_equal(set_time.minute, exp_s_min)
        assert set_time.second <= exp_secs, '%r !< %r' % (set_time.second, exp_secs)



class TestCanopusFromSidingSpring(object):
    '''Integration test: rise/set/transit of a Southern star from a Southern
       observatory.'''

    # IMPORTANT NOTE: LATITUDES
    # IAU convention         = East is +ve
    # Config DB convention   = East is +ve,  e.g. Siding Spring = +149 04 14.13
    # Astro. Alg. convention = West is +ve
    # SLALIB convention (SLA_OBS) = West is +ve
    # SLALIB (rest of library) convention = East is +ve


    def setup(self):

        # Target
        # Note: Aladin units are mas/yr...
        self.canopus = {
             'ra'                : RightAscension('06 23 57.11'),
             'dec'               : Declination('-52 40 03.5'),
             #'ra_proper_motion'  : ProperMotion(RightAscension('00 00 00.01999')),
             #'dec_proper_motion' : ProperMotion(Declination('00 00 00.02367')),
             'ra_proper_motion'  : ProperMotion(RightAscension('00 00 00.0')),
             'dec_proper_motion' : ProperMotion(Declination('00 00 00.0')),
             #'parallax'          : 0.01043,   # Units: arcsec
             'parallax'          : 0.0,   # Units: arcsec
             #'rad_vel'           : 20.5,      # Units: km/s (-ve approaches)
             'rad_vel'           : 0.0,      # Units: km/s (-ve approaches)
             'epoch'             : 2000,
           }

        # Site (East +ve longitude)
        self.siding_spring = {
            'latitude'  : Angle(degrees = -31.273),
            'longitude' : Angle(degrees = 149.070593)
        }

        # Date
        self.date = datetime(year = 2010, month = 3, day = 12)


    def test_rise_set(self):
        (transit, rise, sets) = calc_rise_set(self.canopus, self.siding_spring,
                                             self.date)

        # Second accuracy isn't meaningful, so just check we would round up or
        # down as apppropriate.
        exp_secs = 30

        # Expected times taken from http://aa.usno.navy.mil/data/docs/mrst.php

        # Transit
        exp_t_hr  = 9
        exp_t_min = 8

        transit_time = self.date + transit

        assert_equal(transit_time.hour, exp_t_hr)
        assert_equal(transit_time.minute, exp_t_min)
        msg = '%r !> %r' % (transit_time.second, exp_secs)
        assert transit_time.second <= exp_secs, msg


        # Rise
        exp_r_hr  = 23
        exp_r_min = 27

        rise_time = self.date + rise

        assert_equal(rise_time.hour, exp_r_hr)
        assert_equal(rise_time.minute, exp_r_min)
        assert rise_time.second <= exp_secs, '%r !> %r' % (rise_time.second, exp_secs)


        # Set
        exp_s_hr  = 18
        exp_s_min = 45    # USNO actually gives 18.46, but let's not worry...

        set_time = self.date + sets

        assert_equal(set_time.hour, exp_s_hr)
        assert_equal(set_time.minute, exp_s_min)
        assert set_time.second <= exp_secs, '%r !< %r' % (set_time.second, exp_secs)


class TestNGC2997FromCPT(object):
    '''Integration test: rise/set/transit of a Southern star with telescope horizon.
       Should rise and set, *not* be circumpolar. This test asserts #5969.'''

    def setup(self):

        self.ngc2997 = {
                    'ra' : RightAscension(degrees=146.4116375),
                    'dec': Declination(degrees=-31.19108888)
                    }


        self.cpt = {
                 'latitude' : Angle(degrees=-32.3805542),
                 'longitude': Angle(degrees=20.8101815),
               }


        self.date = datetime(2013, 3, 26)

        self.horizon = Angle(degrees=30)

        self.star = Star(self.cpt['latitude'], self.ngc2997, self.horizon.in_degrees())



    def test_not_circumpolar(self):
        assert(not self.star.is_always_up(self.date))
        assert(not self.star.is_always_down(self.date))





class TestCanopusFromStAndrews(object):
    '''Integration test: rise/set/transit of a Southern star from a Northern
       observatory (never rises).'''

    # IMPORTANT NOTE: LATITUDES
    # IAU convention         = East is +ve
    # Config DB convention   = East is +ve,  e.g. Siding Spring = +149 04 14.13
    # Astro. Alg. convention = West is +ve
    # SLALIB convention (SLA_OBS) = West is +ve
    # SLALIB (rest of library) convention = East is +ve


    def setup(self):

        # Target
        # Note: Aladin units are mas/yr...
        self.canopus = {
             'ra'                : RightAscension('06 23 57.11'),
             'dec'               : Declination('-52 40 03.5'),
             'ra_proper_motion'  : ProperMotion(RightAscension('00 00 00.01999')),
             'dec_proper_motion' : ProperMotion(Declination('00 00 00.02367')),
             'parallax'          : 0.01043,   # Units: arcsec
             'rad_vel'           : 20.5,      # Units: km/s (-ve approaches)
             'epoch'             : 2000,
           }

        # Site (East +ve longitude)
        # Very rough St. Andrews coords, for easy almanac comparison
        self.st_andrews = {
            'latitude'  : Angle(degrees = 56),
            'longitude' : Angle(degrees =  3)
        }

        # Date
        self.date = datetime(year=2010, month=3, day=12)


    @raises(RiseSetError)
    def test_rise_set(self):
        (transit, rise, sets) = calc_rise_set(self.canopus, self.st_andrews,
                                              self.date)


    def test_star_is_always_down(self):

        star = Star(self.st_andrews['latitude'],
                    self.canopus,
                    horizon = 0.0)

        assert(not star.is_always_up(self.date))
        assert(star.is_always_down(self.date))




class TestCapellaFromStAndrews(object):
    '''Integration test: rise/set/transit of a circumpolar star from a Northern
        observatory (never sets).'''

       # IMPORTANT NOTE: LATITUDES
    # IAU convention         = East is +ve
    # Config DB convention   = East is +ve,  e.g. Siding Spring = +149 04 14.13
    # Astro. Alg. convention = West is +ve
    # SLALIB convention (SLA_OBS) = West is +ve
    # SLALIB (rest of library) convention = East is +ve


    def setup(self):
        # Target
        # Note: Aladin units for proper motion are mas/yr...
        self.capella = {
                 'ra'                : RightAscension('05 16 41.36'),
                 'dec'               : Declination('+45 59 52.8'),
                 'ra_proper_motion'  : ProperMotion(RightAscension('00 00 00.07552')),
                 'dec_proper_motion' : ProperMotion(Declination('-00 00 00.42711')),
                 'parallax'          : 0.07729,   # Units: arcsec
                 'rad_vel'           : 30.2,      # Units: km/s (-ve approaches)
                 'epoch'             : 2000,
               }

        # Site (East +ve longitude)
        # Very rough St. Andrews coords, for easy almanac comparison
        self.st_andrews = {
            'latitude'  : Angle(degrees = 56),
            'longitude' : Angle(degrees =  3)
        }

        # Date
        self.date = datetime(year=2010, month=3, day=12)


    @raises(RiseSetError)
    def test_rise_set(self):
        (transit, rise, sets) = calc_rise_set(self.capella, self.st_andrews,
                                              self.date)


    def test_star_is_always_up(self):

        star = Star(self.st_andrews['latitude'],
                    self.capella,
                    horizon = 0.0)

        assert(star.is_always_up(self.date))
        assert(not star.is_always_down(self.date))


    def test_reducing_horizon_keeps_star_down(self):

        star = Star(self.st_andrews['latitude'],
                    self.capella,
                    horizon = 80)

        assert(not star.is_always_down(self.date))


class TestPolarisFromSidingSpring(object):
    '''Integration test: rise/set/transit of a Northern star from a Southern
       observatory (never rises).'''

    # IMPORTANT NOTE: LATITUDES
    # IAU convention         = East is +ve
    # Config DB convention   = East is +ve,  e.g. Siding Spring = +149 04 14.13
    # Astro. Alg. convention = West is +ve
    # SLALIB convention (SLA_OBS) = West is +ve
    # SLALIB (rest of library) convention = East is +ve


    def setup(self):

        # Target
        # Note: Aladin units for proper motion are mas/yr...
        self.polaris = {
                         'ra'                : RightAscension('02 31 49.09'),
                         'dec'               : Declination('+89 15 50.8'),
                       }


        # Site (East +ve longitude)
        self.siding_spring = {
            'latitude'  : Angle(degrees = -31.273),
            'longitude' : Angle(degrees = 149.070593)
        }

        # Date
        self.date = datetime(year = 2010, month = 3, day = 12)



    @raises(RiseSetError)
    def test_rise_set(self):
        (transit, rise, sets) = calc_rise_set(self.polaris, self.siding_spring,
                                              self.date)




class TestMimosaFromSidingSpring(object):
    '''Integration test: rise/set/transit of a circumpolar star from a Southern
        observatory (never sets).'''

       # IMPORTANT NOTE: LATITUDES
    # IAU convention         = East is +ve
    # Config DB convention   = East is +ve,  e.g. Siding Spring = +149 04 14.13
    # Astro. Alg. convention = West is +ve
    # SLALIB convention (SLA_OBS) = West is +ve
    # SLALIB (rest of library) convention = East is +ve


    def setup(self):
        # Target
        # Note: Aladin units for proper motion are mas/yr...
        self.mimosa = {
                         'ra'                : RightAscension('12 47 43.26'),
                         'dec'               : Declination('-59 41 19.549'),
                       }

        # Site (East +ve longitude)
        # Very rough St. Andrews coords, for easy almanac comparison
        self.siding_spring = {
            'latitude'  : Angle(degrees = -31.273),
            'longitude' : Angle(degrees = 149.070593)
        }

        # Date
        self.date = datetime(year=2010, month=3, day=12)


    @raises(RiseSetError)
    def test_rise_set(self):
        (transit, rise, sets) = calc_rise_set(self.mimosa, self.siding_spring,
                                              self.date)


class TestGetAirmassForTarget(object):
    '''Integration test: test getting the airmass for a target at a time. The values these results are compared to
        were obtained from testing the same source/target on PyEphem and Astropy. All have differing values by ~0.1
        but they are all close.
    '''
    def setup(self):
        self.target_1 = make_ra_dec_target(ra=Angle(degrees=148.925583), dec=Angle(degrees=69.673889))
        self.target_2 = make_ra_dec_target(ra=Angle(degrees=68.9791666667), dec=Angle(degrees=16.5))



        self.target_3 = make_ra_dec_target(ra=Angle(degrees=5.392944), dec=Angle(degrees=-69.756111),
                                           ra_proper_motion=ProperMotion(Angle(degrees=0.005519960155599643 / 3600.0, units='arc'), time='year'),
                                           dec_proper_motion=ProperMotion(Angle(degrees=0.000229 / 3600.0, units='arc'), time='year'))

        self.ogg_latitude = Angle(degrees=20.7069444444)
        self.ogg_longitude = Angle(degrees=-156.258055556)
        self.ogg_height = 3065.0 # meters

        self.cpt_latitude = Angle(degrees=-32.3805542)
        self.cpt_longitude = Angle(degrees=20.8101815)
        self.cpt_height = 1804.0 # meters

    def test_get_airmass_for_target_1_ogg(self):
        time = datetime(2016, 5, 20, 23, 12, 32, 247312)

        airmasses = calculate_airmass_at_times([time], self.target_1, self.ogg_latitude, self.ogg_longitude,
                                               self.ogg_height)
        assert_equal(1, len(airmasses))
        assert_almost_equal(2.50, airmasses[0], 2)

    def test_get_airmass_for_target_1_cpt_fail(self):
        time = datetime(2016, 5, 20, 23, 12, 32, 247312)

        airmasses = calculate_airmass_at_times([time], self.target_1, self.cpt_latitude, self.cpt_longitude,
                                               self.cpt_height)
        assert_equal(1, len(airmasses))
        # assert airmass is so high its not visible
        assert_greater(airmasses[0], 10)

    def test_get_airmass_for_target_2_ogg(self):
        time = datetime(2016, 5, 20, 20, 12, 32, 247312)

        airmasses = calculate_airmass_at_times([time], self.target_2, self.ogg_latitude, self.ogg_longitude,
                                               self.ogg_height)
        assert_equal(1, len(airmasses))
        assert_almost_equal(1.33, airmasses[0], 2)

    def test_get_airmass_for_target_2_ogg_fail(self):
        time = datetime(2016, 5, 20, 10, 12, 32, 247312)

        airmasses = calculate_airmass_at_times([time], self.target_2, self.ogg_latitude, self.ogg_longitude,
                                               self.ogg_height)
        assert_equal(1, len(airmasses))
        # assert airmass is so high its not visible
        assert_greater(airmasses[0], 10)

    def test_get_airmass_for_target_2_cpt(self):
        time = datetime(2016, 5, 20, 11, 12, 32, 247312)

        airmasses = calculate_airmass_at_times([time], self.target_2, self.cpt_latitude, self.cpt_longitude,
                                               self.cpt_height)
        assert_equal(1, len(airmasses))
        assert_almost_equal(1.52, airmasses[0], 2)

    def test_get_airmass_for_target_3_cpt(self):
        time = datetime(2016, 5, 20, 3, 12, 32, 247312)

        airmasses = calculate_airmass_at_times([time], self.target_3, self.cpt_latitude, self.cpt_longitude,
                                               self.cpt_height)
        assert_equal(1, len(airmasses))
        assert_almost_equal(1.52, airmasses[0], 2)

