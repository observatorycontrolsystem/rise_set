#!/usr/bin/python
from __future__ import division

from nose.tools import eq_, assert_equal, assert_almost_equal, raises, nottest
import datetime

#Import the module to test
from rise_set.astrometry import *
from rise_set.angle import Angle

class YiannisIsTryingToBreakMyDateCalculator(object):
    def __init__(self):
        self.day   = 20
        self.month = 30
        self.year  = 1988


class TestAstrometry(object):
    '''Unit tests for the astrometry module.'''

    def setUp(self):
        self.date      = datetime.date(year=1988, month=3, day=20)
        self.bad_month = YiannisIsTryingToBreakMyDateCalculator()
        self.mjd       = 47240.0


    def tearDown(self):
        pass


    def test_gregorian_to_ut_mjd(self):
        assert_equal(gregorian_to_ut_mjd(self.date), self.mjd)


    @raises(InvalidDateError)
    def test_gregorian_to_ut_mjd_bad_month(self):
        gregorian_to_ut_mjd(self.bad_month)

#    def test_ut_mjd_to_gmst_in_radians(self):
#        assert_equal(ut_mjd_to_gmst_in_radians(self.mjd), 4)



class TestVenusRiseTransitSet(object):
    '''Implementation of Example 14.a from Astronomical Algorithms, p.99'''

    def setUp(self):
        self.boston = {
           'latitude'  : 42.3333,
           'longitude' : -71.0833
        }

        self.app_sidereal_time = Angle(ra='11 50 58.10')
        self.alpha_2 = Angle(ra='02 46 55.51')
        self.delta_2 = Angle(dec='18 26 27.3')

        self.std_alt = Angle(degrees=-0.5667)


    def test_calc_transit_day_fraction(self):
        m_0 = calc_transit_day_fraction(self.alpha_2, self.boston['longitude'],
                                        self.app_sidereal_time)
        assert_almost_equal(m_0, 0.81965, places=5)


        (hour_angle, msg) = calc_rise_set_hour_angle(self.boston['latitude'],
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


    def setUp(self):

        # Target
        # Note: Aladin units are mas/yr...
        self.deneb = {
                       'ra'                : Angle(ra='20 41 25.91'),
                       'dec'               : Angle(dec='+45 16 49.22'),
                       'ra_proper_motion'  : Angle(ra='00 00 00.00156'),
                       'dec_proper_motion' : Angle(dec='00 00 00.00155'),
                       'parallax'          : 0.00101,  # Units: arcsec
                       'rad_vel'           : -4.5,  # Units: km/s (-ve approaches)
                       'epoch'             : 2000,
                      }

        # Site (East +ve longitude)
        self.maui = {
           'latitude'  : 20.7069444444,
           'longitude' : -156.258055556
        }

        # Date
        self.date = datetime.date(year=2010, month=10, day=25)


    def test_rise_set(self):
        (transit, rise, set) = calc_rise_set(self.deneb, self.maui, self.date)

        # Second accuracy isn't meaningful, so just check we would round up or
        # down as apppropriate.
        exp_secs = 30

        # Expected times taken from http://aa.usno.navy.mil/data/docs/mrst.php

        # Transit
        exp_t_hr  = 4
        exp_t_min = 52   # Not 53, because we expect to round up seconds

        assert_equal(transit[0], exp_t_hr)
        assert_equal(transit[1], exp_t_min)
        assert transit[2] >= exp_secs, '%r !> %r' % (transit[2], exp_secs)


        # Rise
        exp_r_hr  = 21
        exp_r_min = 16   # Not 17, because we expect to round up seconds

        assert_equal(rise[0], exp_r_hr)
        assert_equal(rise[1], exp_r_min)
        assert rise[2] >= exp_secs, '%r !> %r' % (rise[2], exp_secs)


        # Set
        exp_s_hr  = 12
        exp_s_min = 25   # Not 24, because we expect to round down seconds

        assert_equal(set[0], exp_s_hr)
        assert_equal(set[1], exp_s_min)
        assert set[2] <= exp_secs, '%r !< %r' % (set[2], exp_secs)



class TestCanopusFromSidingSpring(object):
    '''Integration test: rise/set/transit of a Southern star from a Southern
       observatory.'''

    # IMPORTANT NOTE: LATITUDES
    # IAU convention         = East is +ve
    # Config DB convention   = East is +ve,  e.g. Siding Spring = +149 04 14.13
    # Astro. Alg. convention = West is +ve
    # SLALIB convention (SLA_OBS) = West is +ve
    # SLALIB (rest of library) convention = East is +ve


    def setUp(self):

        # Target
        # Note: Aladin units are mas/yr...
        self.canopus = {
                         'ra'                : Angle(ra='06 23 57.11'),
                         'dec'               : Angle(dec='-52 40 03.5'),
                         'ra_proper_motion'  : Angle(ra='00 00 00.01999'),
                         'dec_proper_motion' : Angle(dec='00 00 00.02367'),
                         'parallax'          : 0.01043,   # Units: arcsec
                         'rad_vel'           : 20.5,      # Units: km/s (-ve approaches)
                         'epoch'             : 2000,
                       }

        # Site (East +ve longitude)
        self.siding_spring = {
            'latitude'  : -31.273,
            'longitude' : 149.070593
        }

        # Date
        self.date = datetime.date(year = 2010, month = 3, day = 12)


    def test_rise_set(self):
        (transit, rise, set) = calc_rise_set(self.canopus, self.siding_spring,
                                             self.date)

        # Second accuracy isn't meaningful, so just check we would round up or
        # down as apppropriate.
        exp_secs = 30

        # Expected times taken from http://aa.usno.navy.mil/data/docs/mrst.php

        # Transit
        exp_t_hr  = 9
        exp_t_min = 8

        assert_equal(transit[0], exp_t_hr)
        assert_equal(transit[1], exp_t_min)
        assert transit[2] <= exp_secs, '%r !> %r' % (transit[2], exp_secs)


        # Rise
        exp_r_hr  = 23
        exp_r_min = 27

        assert_equal(rise[0], exp_r_hr)
        assert_equal(rise[1], exp_r_min)
        assert rise[2] <= exp_secs, '%r !> %r' % (rise[2], exp_secs)


        # Set
        exp_s_hr  = 18
        exp_s_min = 45    # USNO actually gives 18.46, but let's not worry...

        assert_equal(set[0], exp_s_hr)
        assert_equal(set[1], exp_s_min)
        assert set[2] <= exp_secs, '%r !< %r' % (set[2], exp_secs)



class TestCanopusFromStAndrews(object):
    '''Integration test: rise/set/transit of a Southern star from a Northern
       observatory (never rises).'''

    # IMPORTANT NOTE: LATITUDES
    # IAU convention         = East is +ve
    # Config DB convention   = East is +ve,  e.g. Siding Spring = +149 04 14.13
    # Astro. Alg. convention = West is +ve
    # SLALIB convention (SLA_OBS) = West is +ve
    # SLALIB (rest of library) convention = East is +ve


    def setUp(self):

        # Target
        # Note: Aladin units are mas/yr...
        self.canopus = {
                         'ra'                : Angle(ra='06 23 57.11'),
                         'dec'               : Angle(dec='-52 40 03.5'),
                         'ra_proper_motion'  : Angle(ra='00 00 00.01999'),
                         'dec_proper_motion' : Angle(dec='00 00 00.02367'),
                         'parallax'          : 0.01043,   # Units: arcsec
                         'rad_vel'           : 20.5,      # Units: km/s (-ve approaches)
                         'epoch'             : 2000,
                       }

        # Site (East +ve longitude)
        # Very rough St. Andrews coords, for easy almanac comparison
        self.st_andrews = {
            'latitude'  : 56,
            'longitude' : 3
        }

        # Date
        self.date = datetime.date(year=2010, month=3, day=12)


    @raises(RiseSetError)
    def test_rise_set(self):
        (transit, rise, set) = calc_rise_set(self.canopus, self.st_andrews,
                                             self.date)


class TestCapellaFromStAndrews(object):
    '''Integration test: rise/set/transit of a circumpolar star from a Northern
        observatory (never sets).'''

       # IMPORTANT NOTE: LATITUDES
    # IAU convention         = East is +ve
    # Config DB convention   = East is +ve,  e.g. Siding Spring = +149 04 14.13
    # Astro. Alg. convention = West is +ve
    # SLALIB convention (SLA_OBS) = West is +ve
    # SLALIB (rest of library) convention = East is +ve


    def setUp(self):
        # Target
        # Note: Aladin units for proper motion are mas/yr...
        self.capella = {
                         'ra'                : Angle(ra='05 16 41.36'),
                         'dec'               : Angle(dec='+45 59 52.8'),
                         'ra_proper_motion'  : Angle(ra='00 00 00.07552'),
                         'dec_proper_motion' : Angle(dec='-00 00 00.42711'),
                         'parallax'          : 0.07729,   # Units: arcsec
                         'rad_vel'           : 30.2,      # Units: km/s (-ve approaches)
                         'epoch'             : 2000,
                       }

        # Site (East +ve longitude)
        # Very rough St. Andrews coords, for easy almanac comparison
        self.st_andrews = {
            'latitude'  : 56,
            'longitude' : 3
        }

        # Date
        self.date = datetime.date(year=2010, month=3, day=12)


    @raises(RiseSetError)
    def test_rise_set(self):
        (transit, rise, set) = calc_rise_set(self.capella, self.st_andrews,
                                     self.date)


class TestIntervals(object):

    def setUp(self):

        self.some_adjacent_intervals = [
               # Two contiguous intervals
               (datetime.datetime(year=2011, month=2, day=1, hour=12, minute=0),
                datetime.datetime(year=2011, month=2, day=1, hour=12, minute=30)),
               (datetime.datetime(year=2011, month=2, day=1, hour=12, minute=30),
                datetime.datetime(year=2011, month=2, day=1, hour=13, minute=0)),
               # Three contiguous intervals
               (datetime.datetime(year=2011, month=2, day=1, hour=15, minute=0),
                datetime.datetime(year=2011, month=2, day=1, hour=15, minute=30)),
               (datetime.datetime(year=2011, month=2, day=1, hour=15, minute=30),
                datetime.datetime(year=2011, month=2, day=1, hour=16, minute=0)),
               (datetime.datetime(year=2011, month=2, day=1, hour=16, minute=0),
                datetime.datetime(year=2011, month=2, day=1, hour=16, minute=30)),
               # Two non-contiguous intervals
               (datetime.datetime(year=2011, month=2, day=1, hour=17, minute=0),
                datetime.datetime(year=2011, month=2, day=1, hour=17, minute=30)),
               (datetime.datetime(year=2011, month=2, day=1, hour=18, minute=0),
                datetime.datetime(year=2011, month=2, day=1, hour=18, minute=30)),
              ]


        self.expected_coalescence = [
               # First coalescence
               (datetime.datetime(year=2011, month=2, day=1, hour=12, minute=0),
                datetime.datetime(year=2011, month=2, day=1, hour=13, minute=0)),
               # Second coalescence
               (datetime.datetime(year=2011, month=2, day=1, hour=15, minute=0),
                datetime.datetime(year=2011, month=2, day=1, hour=16, minute=30)),
               # Remaining, non-continous intervals
               (datetime.datetime(year=2011, month=2, day=1, hour=17, minute=0),
                datetime.datetime(year=2011, month=2, day=1, hour=17, minute=30)),
               (datetime.datetime(year=2011, month=2, day=1, hour=18, minute=0),
                datetime.datetime(year=2011, month=2, day=1, hour=18, minute=30)),
            ]

        self.target     = 'sun'
        self.bpl        = dict(latitude = 34.4332222222, longitude = -119.863045833)
        self.dt         = datetime.datetime(year=2011, month=2, day=9)
        self.start_date = datetime.datetime(year=2011, month=2, day=9)
        self.end_date   = datetime.datetime(year=2011, month=2, day=10)
        self.twilight   = 'sunrise'

        self.visibility = Visibility(self.bpl, self.start_date, self.end_date,
                                     self.twilight)


    def test_coalesce_adjacent_intervals(self):
        received = self.visibility.coalesce_adjacent_intervals(
                                                        self.some_adjacent_intervals)

        assert_equal(received, self.expected_coalescence)


    def test_can_get_sun_up_intervals(self):

        expected = [
                     (self.dt, self.dt.replace(hour=1, minute=36, second=1)),
                     (self.dt.replace(hour=14, minute=50, second=42),
                      self.dt.replace(day=10))
                    ]
        received = self.visibility.find_when_target_is_up(self.target, self.dt)

        assert_equal(received, expected)


    def test_can_get_sun_down_intervals(self):

        expected = [
                    (self.dt.replace(hour=1, minute=36, second=1),
                     self.dt.replace(hour=14, minute=50, second=42))
                    ]
        received = self.visibility.find_when_target_is_down(self.target, self.dt)

        assert_equal(received, expected)
