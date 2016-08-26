#!/usr/bin/env python

from __future__ import division

from nose.tools import assert_equal, assert_almost_equals, assert_less, raises
from nose import SkipTest
from datetime import datetime, timedelta

#Import the module to test
from rise_set.visibility import Visibility, set_airmass_limit, InvalidHourAngleLimit

# Additional support modules
from rise_set.angle import Angle
from rise_set.sky_coordinates import RightAscension, Declination
from rise_set.astrometry import make_satellite_target
from rise_set.rates import ProperMotion
from rise_set.moving_objects import initialise_sites
from rise_set.utils          import intersect_many_intervals
from mock import patch


def intervals_almost_equal(received, expected, tolerance=1e-5):
    """
    Check if two interval lists are the same to some tolerance.

    Args:
        received: a list of datetimes
        expected: a list of datetimes
        tolerance: tolerance to say two intervals are close enough (seconds)

    """
    assert len(received) == len(expected)
    for i in range(len(received)):
        assert abs((received[i][0] - expected[i][0]).total_seconds()) < tolerance
        assert abs((received[i][1] - expected[i][1]).total_seconds()) < tolerance


def zero_out_microseconds(received):
    zeroed = []
    for start, end in received:
        z_start = start.replace(microsecond=0)
        z_end   = end.replace(microsecond=0)
        zeroed.append((z_start, z_end))


    return zeroed


class TestIntervals(object):

    def setup(self):

        self.sun        = 'sun'
        self.bpl        = {
                          'latitude'  : Angle(degrees = 34.4332222222),
                          'longitude' : Angle(degrees = -119.863045833)
                          }
        self.dt         = datetime(year=2011, month=2, day=9)
        self.start_date = datetime(year=2011, month=2, day=9)
        self.end_date   = datetime(year=2011, month=2, day=11)
        self.horizon    = 0
        self.twilight   = 'sunrise'


        self.capella    = {
                     'ra'                : RightAscension('05 16 41.36'),
                     'dec'               : Declination('+45 59 52.8'),
                     'ra_proper_motion'  : ProperMotion(RightAscension('00 00 00.07552')),
                     'dec_proper_motion' : ProperMotion(Declination('-00 00 00.42711')),
                     'parallax'          : 0.07729,   # Units: arcsec
                     'rad_vel'           : 30.2,      # Units: km/s (-ve approaches)
                     'epoch'             : 2000,
                   }

        self.rachels_target = {
            'ra'                : RightAscension('00 00 34.23'),
            'dec'               : Declination('-30 45 54.6'),
            'ra_proper_motion'  : ProperMotion(RightAscension('00 00 00.0')),
            'dec_proper_motion' : ProperMotion(Declination('-00 00 00.0')),
            'parallax'          : 0.0,   # Units: arcsec
            'rad_vel'           : 0.0,      # Units: km/s (-ve approaches)
            'epoch'             : 2000,
            }



        self.visibility = Visibility(self.bpl, self.start_date, self.end_date,
                                     self.horizon, self.twilight)


    def test_dark_intervals_have_positive_duration(self):
        site        = {
                        'latitude': Angle(degrees=-30.1673305556),
                        'longitude': Angle(degrees=-70.8046611111),
                      }
        visibility = Visibility(
                                 site=site,
                                 start_date=datetime(2013, 10, 1, 0, 0),
                                 end_date=datetime(2014, 4, 1, 0, 0),
                                 horizon=30,
                                 twilight='nautical',
                               )

        dark_intervals = visibility.get_dark_intervals()
#        dark_intervals = visibility.get_target_intervals('sun', up=False)
#        dark_intervals = visibility.find_when_target_is_up('sun', dt=datetime(2014, 3, 20))

        for start, end in dark_intervals:
            assert_less(start, end)


    def test_can_get_sun_up_intervals(self):

        expected = [
                     (
                       self.dt,
                       self.dt.replace(
                                        hour=1,
                                        minute=36,
                                        second=1,
                                      )
                     ),
                     (
                       self.dt.replace(
                                        hour=14,
                                        minute=50,
                                        second=42,
                                      ),
                        self.dt.replace(day=10)
                     )
                    ]
        received = self.visibility.find_when_target_is_up(self.sun, self.dt)


        # Ignore microseconds for these tests
        received = zero_out_microseconds(received)

        intervals_almost_equal(received, expected, 1e-6)


    def test_can_get_sun_down_intervals(self):

        expected = [
                    (self.dt.replace(hour=1,   minute=36,
                                     second=1,             ),
                     self.dt.replace(hour=14,   minute=50,
                                     second=42,            ))
                    ]
        received = self.visibility.find_when_target_is_down(self.sun, self.dt)

        # Ignore microseconds for these tests
        received = zero_out_microseconds(received)

        assert_equal(received, expected)


    def test_sun_down_intervals_have_positive_duration(self):
        dt = datetime(2013, 10, 31, 0, 0)
        site        = {
                        'latitude': Angle(degrees=-30.1673305556),
                        'longitude': Angle(degrees=-70.8046611111),
                      }
        visibility = Visibility(
                                 site=site,
                                 start_date=datetime(2013, 10, 1, 0, 0),
                                 end_date=datetime(2014, 4, 1, 0, 0),
                                 horizon=30,
                                 twilight='nautical',
                               )

        received = visibility.find_when_target_is_down(self.sun, dt, star=None,
                                                                     airmass=None)

        assert_less(received[0][0], received[0][1])



    def test_can_get_dark_intervals(self):
        # Same as previous test, but with a simpler interface.

        expected = [
                     (
                       datetime(2011, 2, 9, 1, 36, 1),
                       datetime(2011, 2, 9, 14, 50, 42)
                     ),
                     (
                       datetime(2011, 2, 10, 1, 37, 0),
                       datetime(2011, 2, 10, 14, 49, 47)
                     )
                   ]
        received = self.visibility.get_dark_intervals()

        # Ignore microseconds for these tests
        received = zero_out_microseconds(received)

        assert_equal(received, expected)


    def test_get_target_intervals_ra_dec(self):

        expected = [
                     (
                       datetime(2011, 2, 9, 0, 0),
                       datetime(2011, 2, 9, 13, 6, 34)
                     ),
                     (
                       datetime(2011, 2, 9, 18, 52, 15),
                       datetime(2011, 2, 10, 13, 2, 38)
                     ),
                     (
                       datetime(2011, 2, 10, 18, 48, 19),
                       datetime(2011, 2, 11, 0, 0)
                     )
                   ]

        received = self.visibility.get_target_intervals(self.capella, up=True)

        # Ignore microseconds for these tests
        received = zero_out_microseconds(received)

        assert_equal(received, expected)


    @patch('rise_set.visibility.Visibility.get_ra_target_intervals')
    @patch('rise_set.visibility.Visibility.get_moving_object_target_intervals')
    def test_get_target_intervals_no_target_specified(self, mock_func1, mock_func2):
        target = {}

        received = self.visibility.get_target_intervals(target)

        assert_equal(mock_func1.call_count, 0)
        assert_equal(mock_func2.call_count, 1)


    def test_get_target_intervals_satellite_target(self):
        target = make_satellite_target(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        target_intervals = self.visibility.get_target_intervals(target)
        dark_intervals = self.visibility.get_dark_intervals()

        assert_equal(target_intervals, dark_intervals)


    @patch('rise_set.visibility.Visibility.get_ra_target_intervals')
    @patch('rise_set.visibility.Visibility.get_moving_object_target_intervals')
    def test_get_target_intervals_mpc_comet_type_is_a_moving_object(self, moving_obj_func,
                                                                    ra_dec_func):
        target = {'type' : 'mpc_comet'}

        received = self.visibility.get_target_intervals(target)

        assert_equal(moving_obj_func.call_count, 1)
        assert_equal(ra_dec_func.call_count, 0)


    @patch('rise_set.visibility.find_moving_object_up_intervals')
    def test_get_target_intervals_mpc_minor_planet_type_gets_correct_intervals(self, mock_func):
        target = {
                   'type'           : 'mpc_minor_planet',
                   'epoch'          : datetime(2013, 11, 4),
                   'inclination'    : Angle(degrees=7.22565),
                   'long_node'      : Angle(degrees=173.25052),
                   'arg_perihelion' : Angle(degrees=47.60658),
                   'semi_axis'      : 2.4050925,
                   'eccentricity'   : 0.0943494,
                   'mean_anomaly'   : Angle(degrees=54.47380),
                 }

        mock_func.return_value = 1,2

        received = self.visibility.get_target_intervals(target)

        assert_equal(mock_func.call_count, 1)

        args, kwargs = mock_func.call_args

        expected_window = {
                            'start' : datetime(2011, 2, 9, 0, 0),
                            'end'   : datetime(2011, 2, 11, 0, 0)
                          }
        expected_target = target
        expected_site   = {
                            'latitude' : Angle(degrees=34.4332222222),
                            'horizon'  : Angle(degrees=0.0),
                            'longitude': Angle(degrees=-119.863045833)
                          }

        assert_equal(args[0], expected_window)
        assert_equal(args[1], expected_target)
        assert_equal(args[2], expected_site)



    def test_up_intervals_rise_set_transit_within_same_day(self):
        start = datetime(year=2011, month=10, day=13, hour=0, minute=30)
        end   = datetime(year=2011, month=10, day=13, hour=23, minute=30)

        observatory = {
                'latitude'  : Angle(degrees=34.4325),
                'longitude' : Angle('-119 51 46'),
                }

        v = Visibility(observatory, start, end)

        received = v.get_target_intervals(self.rachels_target)
        expected = [
                     (
                       datetime(2011, 10, 13, 2, 8, 2),
                       datetime(2011, 10, 13, 11, 1, 25),
                     )
                    ]

        # Ignore microseconds for these tests
        received = zero_out_microseconds(received)
        intervals_almost_equal(received, expected, tolerance=1e-5)


    def test_edge_interval(self):
        start = datetime(year=2011, month=10, day=13, hour=4, minute=30)
        end   = datetime(year=2011, month=10, day=13, hour=8, minute=30)

        observatory = {
            'latitude'  : Angle(degrees=34.4325),
            'longitude' : Angle('-119 51 46'),
            }

        v = Visibility(observatory, start, end)

        received = v.get_target_intervals(self.rachels_target)
        expected = [
                     (
                       datetime(2011, 10, 13, 4, 30, 0, 0),
                       datetime(2011, 10, 13, 8, 30, 0, 0),
                     )
                    ]

        intervals_almost_equal(received, expected, tolerance=1e-5)


    def test_visibility_end_loop(self):
        start = datetime(year=2011, month=10, day=13, hour=4, minute=30)
        end   = datetime(year=2011, month=10, day=14, hour=0, minute=30)

        observatory = {
            'latitude'  : Angle(degrees=34.4325),
            'longitude' : Angle('-119 51 46'),
            }

        jasons_target = {
            'ra'                : RightAscension('00 00 34.23'),
            'dec'               : Declination('89 45 54.6'),
            'ra_proper_motion'  : ProperMotion(RightAscension('00 00 00.0')),
            'dec_proper_motion' : ProperMotion(Declination('-00 00 00.0')),
            'parallax'          : 0.0,   # Units: arcsec
            'rad_vel'           : 0.0,      # Units: km/s (-ve approaches)
            'epoch'             : 2000,
            }

        v = Visibility(observatory, start, end)

        received = v.get_target_intervals(jasons_target)
        expected = [
                     (
                       datetime(2011, 10, 13, 4, 30, 0, 0),
                       datetime(2011, 10, 14, 0, 30, 0, 0),
                     )
                    ]

        intervals_almost_equal(received, expected, tolerance=1e-5)


    def test_hour_angle(self):
        expected = [(datetime(2013, 3, 22, 8, 24, 37, 986301),
                     datetime(2013, 3, 22, 18, 22, 59, 690688))]
        target = {
            'ra'                : RightAscension('20 41 25.91'),
            'dec'               : Declination('+20 00 00.00'),
            'epoch'             : 2000,
            }

        site  = {
                  'latitude':  Angle(degrees=-30.1673305556),
                  'longitude': Angle(degrees=-70.8046611111),
                }

        start_date = datetime(year=2013, month=3, day=22)
        end_date   = datetime(year=2013, month=3, day=23)

        v = Visibility(site, start_date, end_date,
                       ha_limit_neg=-5.0, ha_limit_pos=5.0)

        received = v.get_ha_intervals(target)

        intervals_almost_equal(received, expected, 1e-5)

    def test_ha_wrong_day(self):
        # for some windows/limits, the HA block did not start at the beginning of the window.
        # This test fails prior to 2013-02-20
        expected = [(datetime(2011, 11, 01, 06, 00, 00, 000000),datetime(2011, 11, 01, 07, 52, 00, 564199)),
                    (datetime(2011, 11, 02, 02, 01, 50, 423880),datetime(2011, 11, 02, 06, 00, 00, 000000))]

        target = {
            'ra' : RightAscension(degrees=310.35795833333333),
            'dec': Declination(degrees=45.280338888888885),
            'epoch': 2000,
            }

        site  = {
            'latitude':  Angle(degrees=34.433157),
            'longitude': Angle(degrees=-119.86308),
            }

        start_date = datetime(year=2011, month=11, day=1, hour=6)
        end_date   = datetime(year=2011, month=11, day=2, hour=6)

        v = Visibility(site, start_date, end_date, twilight='nautical', horizon=25,
                       ha_limit_neg=-12, ha_limit_pos=12)


#        # get the intervals of each separately
#        dark               = v.get_dark_intervals()
#        above_horizon      = v.get_target_intervals(target)
#        within_hour_angle  = v.get_ha_intervals(target)
#        # find the overlapping intervals between them
#        received = intersect_many_intervals(dark, above_horizon, within_hour_angle)

        received = v.get_observable_intervals(target, moon_distance=Angle(degrees=0))

        intervals_almost_equal(received, expected, tolerance=1e-5)


    def test_ha_gets_all_intervals(self):
        # For targets that cross 00 UTC, there was a missed interval
        # This failed prior to 2015-03-19
        expected = [(datetime(2015, 3, 16, 0, 5, 58, 794332),datetime(2015, 3, 16, 9, 16, 28, 362369)),
                    (datetime(2015, 3, 17, 0, 2, 2, 884862),datetime(2015, 3, 17, 9, 12, 32, 452899)),
                    (datetime(2015, 3, 17, 23, 58, 6, 975392),datetime(2015, 3, 18, 9, 8, 36,543429)),
                    (datetime(2015, 3, 18, 23, 54, 11, 65922),datetime(2015, 3, 19, 9, 4, 40, 633959)),
                    (datetime(2015, 3, 19, 23, 50, 15, 156452),datetime(2015, 3, 20, 9, 0, 44, 724489)),
                    (datetime(2015, 3, 20, 23, 46, 19, 246982),datetime(2015, 3, 21, 8, 56, 48, 815019))]

        target = {
            'ra' : RightAscension(degrees=172.962037037),
            'dec': Declination(degrees=-12.5255555556),
            'epoch': 2000,
            }

        site  = {
            'latitude':  Angle(degrees=-30.1673306),
            'longitude': Angle(degrees=-70.8046611),
            }

        start_date = datetime(year=2015, month=3, day=15, hour=18)
        end_date   = datetime(year=2015, month=3, day=21, hour=18)

        v = Visibility(site, start_date, end_date, twilight='nautical', horizon=30,
                       ha_limit_neg=-4.6, ha_limit_pos=4.6)

        received = v.get_ha_intervals(target)
        intervals_almost_equal(received, expected, tolerance=1e-5)


    def test_short_window_has_all_intervals(self):

        window_start = datetime(2016, 3, 14, 21, 35, 7, 488985)

        window_end = datetime(2016, 3, 15, 21, 35, 7, 488985)

        sitecoords_coj = {'latitude': Angle(-31.2733), 'longitude': Angle(149.439)}

        eff_horizon = 29.999999999999993

        target = {'rad_vel': 0.0, 'ra': RightAscension(degrees=91.784167),
                  'dec_proper_motion': 0.0, 'ra_proper_motion': 0.0,
                  'dec': Declination(degrees=-45.181167), 'type': 'SIDEREAL', 'epoch': 2000, 'parallax': 0}

        v1 = Visibility(site=sitecoords_coj, start_date=window_start, end_date=window_end,
                        twilight='nautical', horizon=eff_horizon)
        observable_intervals1 = v1.get_observable_intervals(target, moon_distance=Angle(degrees=0))

        v2 = Visibility(site=sitecoords_coj, start_date=window_start,
                        end_date=window_end + timedelta(days=1),
                        twilight='nautical', horizon=eff_horizon)

        observable_intervals2 = v2.get_observable_intervals(target, moon_distance=Angle(degrees=0))

        intervals_almost_equal(observable_intervals1, observable_intervals2[:-1], tolerance=1e-5)


    def test_short_window_has_all_intervals_lsc(self):
        window_start = datetime(2016, 3, 14, 21, 35, 7, 488985)

        window_end = datetime(2016, 3, 15, 21, 35, 7, 488985)

        sitecoords_lsc = {'latitude': Angle(degrees=-30.167367),
                          'longitude': Angle(degrees=-70.8049)}

        eff_horizon = 29.999999999999993
        target = {'rad_vel': 0.0, 'ra': RightAscension(degrees=91.784167),
                  'dec_proper_motion': 0.0, 'ra_proper_motion': 0.0,
                  'dec': Declination(degrees=-45.181167), 'type': 'SIDEREAL', 'epoch': 2000,
                  'parallax': 0}

        v1 = Visibility(site=sitecoords_lsc, start_date=window_start, end_date=window_end,
                        twilight='nautical', horizon=eff_horizon)

        observable_intervals1 = v1.get_observable_intervals(target, moon_distance=Angle(degrees=0))

        v2 = Visibility(site=sitecoords_lsc, start_date=window_start,
                        end_date=window_end + timedelta(days=1),
                        twilight='nautical', horizon=eff_horizon)

        observable_intervals2 = v2.get_observable_intervals(target, moon_distance=Angle(degrees=0))

        assert len(observable_intervals2) == 2
        intervals_almost_equal(observable_intervals1, observable_intervals2[:-1], tolerance=1e-5)

        v3 = Visibility(site=sitecoords_lsc, start_date=window_start,
                        end_date=window_end + timedelta(days=2),
                        twilight='nautical', horizon=eff_horizon)
        observable_intervals3 = v3.get_observable_intervals(target, moon_distance=Angle(degrees=0))
        assert len(observable_intervals3) == 3
        intervals_almost_equal(observable_intervals2, observable_intervals3[:-1])


    def test_ha_where_interval_close_to_a_day(self):
        site = {'latitude': Angle(degrees=-31.273),
                'longitude': Angle(degrees=149.070593)}

        target = {'rad_vel': 0.0, 'ra': RightAscension(degrees=251.3),
                  'dec_proper_motion': 0.0, 'ra_proper_motion': 0.0,
                  'dec': Declination(degrees=27.874), 'type': 'SIDEREAL', 'epoch': 2000,
                  'parallax': 0}
        start_date = datetime(2016, 6, 18, 22, 9, 6)
        end_date = datetime(2016, 6, 19, 22, 9, 6)

        v1 = Visibility(site, start_date, end_date, horizon=15, ha_limit_neg=-4.6, ha_limit_pos=4.6)

        hi1 = v1.get_ha_intervals(target)
        # returns no intervals before 2016-06-21

        end_date = datetime(2016, 6, 19, 23, 9, 6)
        v2 = Visibility(site, start_date, end_date, horizon=15, ha_limit_neg=-4.6, ha_limit_pos=4.6)

        hi2 = v2.get_ha_intervals(target)
        # returns an interval that should be within the previous set as well
        intervals_almost_equal(hi1, hi2, tolerance=1e-5)


    def test_ha_ut_mjd_is_truncated(self):
        # If a user window specifies a time, then hour angle calculations are entirely
        # wrong unless that time is discarded. This unit test protects against regressions
        # of this bug fix.
        target = {
            'ra'    : RightAscension(degrees=59.890410785),
            'dec'   : Declination(degrees=-89.2641974447),
            'epoch' : 2000,
            }

        site  = {
                  'name'         : 'coj',
                  'latitude'     : Angle(degrees=-31.273),
                  'longitude'    : Angle(degrees=149.070593),
                  'horizon'      : 30,
                  'ha_limit_neg' : -4.6,
                  'ha_limit_pos' : 4.6,
                }

        start_date = datetime(year=2014, month=6, day=1, hour=12, minute=10)
        end_date   = datetime(year=2014, month=6, day=1, hour=12, minute=15)

        v = Visibility(site, start_date, end_date,
                       twilight='nautical',
                       horizon=site['horizon'],
                       ha_limit_neg=site['ha_limit_neg'],
                       ha_limit_pos=site['ha_limit_pos'],
                       )

        assert_equal(v.get_observable_intervals(target, moon_distance=Angle(degrees=0)), [])



    @raises(InvalidHourAngleLimit)
    def test_negative_hour_angle_too_big_is_rejected(self):
        visibility = Visibility(None, None, None, ha_limit_neg=-13.0)

    @raises(InvalidHourAngleLimit)
    def test_negative_hour_angle_too_small_is_rejected(self):
        visibility = Visibility(None, None, None, ha_limit_neg=1.0)

    @raises(InvalidHourAngleLimit)
    def test_positive_hour_angle_too_big_is_rejected(self):
        visibility = Visibility(None, None, None, ha_limit_pos=13.0)

    @raises(InvalidHourAngleLimit)
    def test_positive_hour_angle_too_small_is_rejected(self):
        visibility = Visibility(None, None, None, ha_limit_pos=-1.0)


    def test_sites(self):
        site_filename="test/telescopes.dat"

        # expected for dec = -60
        expected = {
                     '1m0a.doma.elp' : [],
                     '1m0a.doma.coj' : [(datetime(2013, 3, 22, 17, 49, 35, 118039),
                                         datetime(2013, 3, 22, 20,  8, 49, 303466))],
                     '1m0a.doma.cpt' : [(datetime(2013, 3, 22,  2, 25,  9, 494210),
                                         datetime(2013, 3, 22,  4, 41, 27, 874760))],
                     '1m0a.doma.lsc' : [(datetime(2013, 3, 22,  8, 30, 37,   6004),
                                         datetime(2013, 3, 22, 10, 48,  0, 967403))]
                   }

        target = {
                   'ra'                : RightAscension('20 41 25.91'),
                   'dec'               : Declination('-60 00 00.00'),
                   'epoch'             : 2000,
                 }

        start_date = datetime(year=2013, month=3, day=22)
        end_date   = datetime(year=2013, month=3, day=23)
        sites      = initialise_sites(site_filename)
        received   = []

        for site in sites:
            v = Visibility(site, start_date, end_date, ha_limit_neg=-4.9, ha_limit_pos=4.9)
            intervals = v.get_observable_intervals(target, moon_distance=Angle(degrees=0))

            intervals_almost_equal(intervals, expected[site['name']], tolerance=1e-5)



class TestAirmassCalculation(object):

    def setup(self):
        pass

    def interval_for_airmass(self, airmass):
        horizon = 25

        # Site (East +ve longitude)
        site = {
           'latitude'  : Angle(degrees = 20.0),
           'longitude' : Angle(degrees = -150.0)
        }

        start = datetime(year=2010, month=10, day=25)
        end   = datetime(year=2010, month=10, day=26)

        # Target
        # Note: Aladin units are mas/yr...
        target = {
               'ra'                : RightAscension('20 41 25.91'),
               'dec'               : Declination('+20 00 00.00'),
               'epoch'             : 2000,
              }

        v = Visibility(site, start, end, horizon)
        intervals   = v.get_target_intervals(target=target, airmass=airmass)
        up_interval = intervals[0][1] - intervals[0][0]

        return up_interval

    def test_set_airmass_limit_no_airmass(self):
        horizon  = 30
        expected = horizon
        assert_equal(set_airmass_limit(None, horizon), expected)


    def test_set_airmass_limit_airmass_worse_than_horizon(self):
        horizon  = 30
        airmass  = 3
        expected = horizon
        assert_equal(set_airmass_limit(airmass, horizon), expected)


    def test_set_airmass_limit_airmass_better_than_horizon(self):
        horizon  = 30
        airmass  = 1.2
        expected = 56.44
        assert_almost_equals(set_airmass_limit(airmass, horizon), expected, places=2)


    def test_airmass_is_applied_if_above_horizon_and_not_otherwise(self):
        airmass_above_horizon = 1.082392200292394   # 22.5 degrees from the zenith
        airmass_below_horizon = 3
        interval_with_airmass    = self.interval_for_airmass(airmass_above_horizon)
        interval_without_airmass = self.interval_for_airmass(airmass_below_horizon)

        assert_less(interval_with_airmass, interval_without_airmass)



