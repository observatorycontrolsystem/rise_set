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



class TestVenusRiseTransitSet(object):
    '''Implementation of Example 14.a from Astronomical Algorithms, p.99'''

    def setUp(self):
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
        self.date = datetime.datetime(year=2010, month=10, day=25)


    def test_rise_set(self):
        (transit, rise, set) = calc_rise_set(self.deneb, self.maui, self.date)

        # Second accuracy isn't meaningful, so just check we would round up or
        # down as apppropriate.
        exp_secs = 30

        # Expected times taken from http://aa.usno.navy.mil/data/docs/mrst.php

        # Transit
        exp_t_hr  = 4
        exp_t_min = 52   # Not 53, because we expect to round up seconds

        transit_time = self.date + transit

        assert_equal(transit_time.hour, exp_t_hr)
        assert_equal(transit_time.minute, exp_t_min)
        assert transit_time.second >= exp_secs, '%r !> %r' % (transit_time.second, exp_secs)


        # Rise
        exp_r_hr  = 21
        exp_r_min = 16   # Not 17, because we expect to round up seconds

        rise_time = self.date + rise

        assert_equal(rise_time.hour, exp_r_hr)
        assert_equal(rise_time.minute, exp_r_min)
        assert rise_time.second >= exp_secs, '%r !> %r' % (rise_time.second, exp_secs)


        # Set
        exp_s_hr  = 12
        exp_s_min = 25   # Not 24, because we expect to round down seconds

        set_time = self.date + set

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


    def setUp(self):

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
        self.date = datetime.datetime(year = 2010, month = 3, day = 12)


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

        transit_time = self.date + transit

        assert_equal(transit_time.hour, exp_t_hr)
        assert_equal(transit_time.minute, exp_t_min)
        assert transit_time.second <= exp_secs, '%r !> %r' % (transit_time.second, exp_secs)


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

        set_time = self.date + set

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


        self.date = datetime.date(2013, 3, 26)

        self.horizon=Angle(degrees=30)

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


    def setUp(self):

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
        self.date = datetime.date(year=2010, month=3, day=12)


    @raises(RiseSetError)
    def test_rise_set(self):
        (transit, rise, set) = calc_rise_set(self.canopus, self.st_andrews,
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


    def setUp(self):
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
        self.date = datetime.date(year=2010, month=3, day=12)


    @raises(RiseSetError)
    def test_rise_set(self):
        (transit, rise, set) = calc_rise_set(self.capella, self.st_andrews,
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


    def setUp(self):

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
        self.date = datetime.date(year = 2010, month = 3, day = 12)



    @raises(RiseSetError)
    def test_rise_set(self):
        (transit, rise, set) = calc_rise_set(self.polaris, self.siding_spring,
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


    def setUp(self):
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
        self.date = datetime.date(year=2010, month=3, day=12)


    @raises(RiseSetError)
    def test_rise_set(self):
        (transit, rise, set) = calc_rise_set(self.mimosa, self.siding_spring,
                                             self.date)
