#!/usr/bin/env python

from __future__ import division
from builtins import range
from builtins import object

from nose.tools import assert_equal, assert_almost_equals, assert_less, raises
from nose import SkipTest
from datetime import datetime, timedelta

#Import the module to test
from rise_set.visibility import Visibility, set_airmass_limit, InvalidHourAngleLimit

# Additional support modules
from rise_set.angle import Angle
from rise_set.sky_coordinates import RightAscension, Declination
from rise_set.astrometry import (make_satellite_target, make_minor_planet_target, make_comet_target,
                                 make_major_planet_target, make_hour_angle_target)
from rise_set.rates import ProperMotion
from rise_set.moving_objects import initialise_sites
from rise_set.utils          import intersect_many_intervals, coalesce_adjacent_intervals
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

        self.site       = {
            'latitude': Angle(degrees=34.4332222222),
            'longitude': Angle(degrees=-119.863045833),
            'ha_limit_neg': Angle(degrees=-4.6 * 15.0),
            'ha_limit_pos': Angle(degrees=4.6 * 15.0),
            'horizon': Angle(degrees=15.0)
        }


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

        # Details from JPL Horizons for Jupiter
        self.major_planet_target = make_major_planet_target('JPL_MAJOR_PLANET',
                                                            epochofel=55959.0,
                                                            inclination=1.303884172546506,
                                                            long_node=100.5093329813755,
                                                            arg_perihelion=274.0516181838379,
                                                            semi_axis=5.204023536751508,
                                                            eccentricity=0.04910768996084790,
                                                            mean_anomaly=26.60707699766562,
                                                            dailymot=0.08306200006467207)

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

    def test_get_target_intervals_hour_angle_target(self):
        target = make_hour_angle_target(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

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
        expected = [(datetime(2013, 3, 22, 8, 25, 13, 564442),
                     datetime(2013, 3, 22, 18, 23, 35, 268830))]
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


    def test_hour_angle_near_pole(self):
        expected = [(datetime(2018, 7, 4, 8, 33, 18, 130142),
                     datetime(2018, 7, 4, 17, 31, 49, 663014))]
        target = {
            'ra': RightAscension(degrees=260.1633416667),
            'dec': Declination(degrees=-89.0275836111),
            'epoch': 2000,
        }

        site = {
            'latitude': Angle(degrees=-31.272932),
            'longitude': Angle(degrees=149.070648),
        }

        start_date = datetime(year=2018, month=7, day=4)
        end_date = datetime(year=2018, month=7, day=5)

        v = Visibility(site, start_date, end_date,
                       ha_limit_neg=-4.533333, ha_limit_pos=4.4666667)

        received = v.get_ha_intervals(target)

        intervals_almost_equal(received, expected, 1e-5)


    def test_ha_wrong_day(self):
        # for some windows/limits, the HA block did not start at the beginning of the window.
        # This test fails prior to 2013-02-20

        expected = [(datetime(2011, 11, 1, 6, 0, 0, 0),datetime(2011, 11, 1, 7, 52, 0, 564199)),
                    (datetime(2011, 11, 2, 2, 1, 50, 423880),datetime(2011, 11, 2, 6, 0, 0, 0))]

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
        expected = [(datetime(2015, 3, 16, 0, 6, 46, 500233),datetime(2015, 3, 16, 9, 17, 16, 68270)),
                    (datetime(2015, 3, 17, 0, 2, 50, 590764),datetime(2015, 3, 17, 9, 13, 20, 158801)),
                    (datetime(2015, 3, 17, 23, 58, 54, 681295),datetime(2015, 3, 18, 9, 9, 24, 249332)),
                    (datetime(2015, 3, 18, 23, 54, 58, 771826),datetime(2015, 3, 19, 9, 5, 28, 339863)),
                    (datetime(2015, 3, 19, 23, 51, 2, 862357),datetime(2015, 3, 20, 9, 1, 32, 430394)),
                    (datetime(2015, 3, 20, 23, 47, 6, 952888),datetime(2015, 3, 21, 8, 57, 36, 520925))]

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
        print(received)

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
                     '1m0a.doma.coj' : [(datetime(2013, 3, 22, 17, 50, 38, 57999),
                                         datetime(2013, 3, 22, 20, 8, 49, 303466))],
                     '1m0a.doma.cpt' : [(datetime(2013, 3, 22,  2, 26,  12, 434169),
                                         datetime(2013, 3, 22,  4, 41, 27, 874760))],
                     '1m0a.doma.lsc' : [(datetime(2013, 3, 22, 8, 31, 39, 945963),
                                         datetime(2013, 3, 22, 10, 48, 0, 967403))]
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
            print(intervals)

            intervals_almost_equal(intervals, expected[site['name']], tolerance=1e-5)

    def test_target_up_intervals_jpl_major_planet(self):
        v = Visibility(self.site, datetime(2012, 2, 1), datetime(2012, 2, 2), ha_limit_neg=-4.6, ha_limit_pos=4.6)
        observable_intervals = v.get_observable_intervals(self.major_planet_target, moon_distance=Angle(degrees=0))

        # Jupiter is above the horizon and within ha_limits from rise time 1:27 until 5:55
        observable_intervals = coalesce_adjacent_intervals(observable_intervals)
        rise_time = datetime(2012, 2, 1, 1, 27, 51, 357328)
        assert_equal(observable_intervals[0][0], rise_time)
        assert_equal(observable_intervals[0][1], datetime(2012, 2, 1, 5, 45))

    def test_target_up_intervals_jpl_major_planet_high_ha(self):
        site = self.site.copy()
        site['ha_limit_neg'] = Angle(degrees=-360.0)
        site['ha_limit_pos'] = Angle(degrees=360.0)
        site['horizon'] = Angle(degrees=0.0)
        v = Visibility(site, datetime(2012, 2, 1), datetime(2012, 2, 2), ha_limit_neg=-12, ha_limit_pos=12,
                       horizon=0)
        observable_intervals = v.get_observable_intervals(self.major_planet_target, moon_distance=Angle(degrees=0))

        # Jupiter is above the horizon of 0 degrees from rise time 1:27 until 7:45
        observable_intervals = coalesce_adjacent_intervals(observable_intervals)
        rise_time = datetime(2012, 2, 1, 1, 27, 51, 357328)

        assert_equal(observable_intervals[0][0], rise_time)
        assert_equal(observable_intervals[0][1], datetime(2012, 2, 1, 7, 45))

    def test_target_up_intervals_jpl_major_planet_high_ha_and_no_horizon(self):
        site = self.site.copy()
        site['ha_limit_neg'] = Angle(degrees=-360.0)
        site['ha_limit_pos'] = Angle(degrees=360.0)
        site['horizon'] = Angle(degrees=-360.0)
        v = Visibility(site, datetime(2012, 2, 1), datetime(2012, 2, 2), ha_limit_neg=-12, ha_limit_pos=12,
                       horizon=-360)
        observable_intervals = v.get_observable_intervals(self.major_planet_target, moon_distance=Angle(degrees=0))
        dark_intervals = v.get_dark_intervals()
        # with maxed horizon and ha_limits, the observable intervals basically become the dark intervals
        observable_intervals = coalesce_adjacent_intervals(observable_intervals)

        assert_equal(observable_intervals[0][0], dark_intervals[0][0])
        assert_equal(observable_intervals[0][1], dark_intervals[0][1])


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


class TestMoonDistanceCalculation(object):
    ''' All time intervals to test against are obtained from JPL Horizons
    '''
    def setup(self):
        self.site = {
            'latitude'  : Angle(degrees = 20.0),
            'longitude' : Angle(degrees = -150.0),
            'ha_limit_neg': Angle(degrees=-4.6*12.0),
            'ha_limit_pos': Angle(degrees=4.6*12.0),
            'horizon': Angle(degrees=15.0)
        }

        self.horizon = 15.0
        self.sidereal_target = {
            'ra'                : RightAscension('20 41 25.91'),
            'dec'               : Declination('+20 00 00.00'),
            'epoch'             : 2000,
        }

        # Details from JPL Horizons for Ceres
        self.minor_planet_target = make_minor_planet_target('MPC_MINOR_PLANET',
                                                            epoch=49731,
                                                            inclination=10.60069567603618,
                                                            long_node=80.65851514365535,
                                                            arg_perihelion=71.44921526124109,
                                                            semi_axis=2.767218108003098,
                                                            eccentricity=0.07610292126891821,
                                                            mean_anomaly=340.389084821267)

        # Details from JPL Horizons for Comet 27P
        self.comet_target = make_comet_target('MPC_COMET',
                                              epoch=56364,
                                              epochofperih=55776.8910902,
                                              inclination=28.96687278723059,
                                              long_node=250.6264098390235,
                                              arg_perihelion=196.0253968913816,
                                              perihdist=0.748287467144728,
                                              eccentricity=0.9189810923126022)

        # Details from JPL Horizons for Jupiter
        self.major_planet_target = make_major_planet_target('JPL_MAJOR_PLANET',
                                                            epochofel=55959.0,
                                                            inclination=1.303884172546506,
                                                            long_node=100.5093329813755,
                                                            arg_perihelion=274.0516181838379,
                                                            semi_axis=5.204023536751508,
                                                            eccentricity=0.04910768996084790,
                                                            mean_anomaly=26.60707699766562,
                                                            dailymot=0.08306200006467207)


    def test_moon_distance_no_angle(self):
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)
        moon_distance_constraint = Angle(degrees=0)

        v = Visibility(self.site, start, end, self.horizon)
        target_intervals   = v.get_target_intervals(target=self.sidereal_target)
        moon_distance_intervals = v.get_moon_distance_intervals(self.sidereal_target, target_intervals, moon_distance_constraint)

        assert_equal(moon_distance_intervals, target_intervals)

    def test_moon_distance_none_removed(self):
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)
        # low angle given, the distance is always greater for this target so it should allow all intervals
        moon_distance_constraint = Angle(degrees=30)

        v = Visibility(self.site, start, end, self.horizon)
        target_intervals   = v.get_target_intervals(target=self.sidereal_target)
        moon_distance_intervals = v.get_moon_distance_intervals(self.sidereal_target, target_intervals, moon_distance_constraint)

        assert_equal(moon_distance_intervals, target_intervals)

    def test_moon_distance_half_removed(self):
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)
        # this target/site/date has a morning and evening target interval on this day
        # the morning interval has a moon/target distance < 70, and then evening interval is >70
        moon_distance_constraint = Angle(degrees=70)

        v = Visibility(self.site, start, end, self.horizon)
        target_intervals   = v.get_target_intervals(target=self.sidereal_target)
        moon_distance_intervals = v.get_moon_distance_intervals(self.sidereal_target, target_intervals, moon_distance_constraint)

        # test that just the evening interval is returned by the moon_distance_intervals
        assert_equal(moon_distance_intervals, [target_intervals[1]])

    def test_moon_distance_all_removed(self):
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)
        # this target/site/date has no moon distances less than 75 degrees
        moon_distance_constraint = Angle(degrees=75)

        v = Visibility(self.site, start, end, self.horizon)
        target_intervals   = v.get_target_intervals(target=self.sidereal_target)
        moon_distance_intervals = v.get_moon_distance_intervals(self.sidereal_target, target_intervals, moon_distance_constraint)

        # test that no moon distance intervals are returned due to the constraint of > 75 degrees distance
        assert_equal(moon_distance_intervals, [])

    def test_moon_distance_ignored_for_satellite_target(self):
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)
        v = Visibility(self.site, start, end, self.horizon)

        target = make_satellite_target(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        observable_intervals = v.get_observable_intervals(target)
        # check that this doesn't crash, and that intervals match night intervals
        night_intervals = v.get_dark_intervals()
        assert_equal(night_intervals, observable_intervals)

    def test_moon_distance_ignored_for_hour_angle_target(self):
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)
        v = Visibility(self.site, start, end, self.horizon)

        target = make_hour_angle_target(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        observable_intervals = v.get_observable_intervals(target)
        # check that this doesn't crash, and that intervals match night intervals
        night_intervals = v.get_dark_intervals()
        assert_equal(night_intervals, observable_intervals)

    def test_moon_distance_minor_planet_half_removed(self):
        start = datetime(2012, 8, 2)
        end = datetime(2012, 8, 3)
        v = Visibility(self.site, start, end, self.horizon)

        target = self.minor_planet_target.copy()
        # According to JPL horizons, this target has a moon distance > 115 only until 21:00
        moon_distance_constraint = Angle(degrees=115)

        # The target intervals start at 14:30 and go until 21:45
        target_intervals = v.get_target_intervals(target=target)
        moon_distance_intervals = v.get_moon_distance_intervals(target, target_intervals, moon_distance_constraint)

        # Verify that moon distance intervals only goes until 21:00
        assert_equal(moon_distance_intervals[0][0], target_intervals[0][0])
        assert_equal(moon_distance_intervals[0][1], datetime(2012, 8, 2, 21, 0))
        assert_equal(len(moon_distance_intervals), 1)

    def test_moon_distance_minor_planet_half_removed_2(self):
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)
        v = Visibility(self.site, start, end, self.horizon)

        target = self.minor_planet_target.copy()
        # According to JPL horizons, this target has a moon distance > 30 after ~6:15
        moon_distance_constraint = Angle(degrees=30)

        # The target intervals start at 00:00 and go until 6:45 and then start again at 23:30 until 24:00
        target_intervals = v.get_target_intervals(target=target)
        moon_distance_intervals = v.get_moon_distance_intervals(target, target_intervals, moon_distance_constraint)

        # Verify that moon distance intervals only appear from 6:00 to 6:45 and from 23:30 to 24:00
        assert_equal(len(moon_distance_intervals), 2)
        assert_equal(moon_distance_intervals[0][0], datetime(2012, 1, 2, 6, 15))
        assert_equal(moon_distance_intervals[0][1], datetime(2012, 1, 2, 6, 45))
        assert_equal(moon_distance_intervals[1][0], datetime(2012, 1, 2, 23, 30))
        assert_equal(moon_distance_intervals[1][1], datetime(2012, 1, 3, 0, 0))

    def test_moon_distance_minor_planet_none_removed(self):
        start = datetime(2012, 2, 2)
        end = datetime(2012, 2, 3)
        v = Visibility(self.site, start, end, self.horizon)

        target = self.minor_planet_target.copy()
        # According to JPL horizons, this target has a moon distance > 50 during the whole day
        moon_distance_constraint = Angle(degrees=50)

        # The target intervals start at 00:00 and go until 5:15 and then start again at 22:00 until 24:00
        target_intervals = v.get_target_intervals(target=target)
        moon_distance_intervals = v.get_moon_distance_intervals(target, target_intervals, moon_distance_constraint)

        # Verify that moon distance intervals only appear during the full target interval
        assert_equal(len(moon_distance_intervals), 2)
        assert_equal(moon_distance_intervals[0][0], datetime(2012, 2, 2, 0, 0))
        assert_equal(moon_distance_intervals[0][1], datetime(2012, 2, 2, 5, 15))
        assert_equal(moon_distance_intervals[1][0], datetime(2012, 2, 2, 22, 00))
        assert_equal(moon_distance_intervals[1][1], datetime(2012, 2, 3, 0, 0))

    def test_moon_distance_minor_planet_all_removed(self):
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)
        v = Visibility(self.site, start, end, self.horizon)

        target = self.minor_planet_target.copy()
        # According to JPL horizons, this target has a moon distance <40 all the time
        moon_distance_constraint = Angle(degrees=40)

        # The target intervals start at 00:00 and go until 6:45 and then start again at 23:30 until 24:00
        target_intervals = v.get_target_intervals(target=target)
        moon_distance_intervals = v.get_moon_distance_intervals(target, target_intervals, moon_distance_constraint)

        # Verify that moon distance intervals only appear from 6:00 to 6:45 and from 23:30 to 24:00
        assert_equal(len(moon_distance_intervals), 0)

    def test_moon_distance_comet_non_removed(self):
        start = datetime(2012, 7, 22)
        end = datetime(2012, 7, 23)
        v = Visibility(self.site, start, end, self.horizon)

        target = self.comet_target.copy()
        # According to JPL horizons, this target has a moon distance >45 all the time
        moon_distance_constraint = Angle(degrees=45)

        # The target intervals start at 00:15 to 07:15
        target_intervals = v.get_target_intervals(target=target)
        moon_distance_intervals = v.get_moon_distance_intervals(target, target_intervals, moon_distance_constraint)

        print(target_intervals)
        print(moon_distance_intervals)

        # Verify that moon distance intervals is the entire target interval
        assert_equal(len(moon_distance_intervals), 1)
        assert_equal(moon_distance_intervals[0][0], target_intervals[0][0])
        assert_equal(moon_distance_intervals[0][1], target_intervals[-1][1])

    def test_moon_distance_comet_all_removed(self):
        start = datetime(2012, 7, 22)
        end = datetime(2012, 7, 23)
        v = Visibility(self.site, start, end, self.horizon)

        target = self.comet_target.copy()
        # According to JPL horizons, this target has a moon distance <65 all the time
        moon_distance_constraint = Angle(degrees=65)

        # The target intervals start at 00:15 to 07:15
        target_intervals = v.get_target_intervals(target=target)
        moon_distance_intervals = v.get_moon_distance_intervals(target, target_intervals, moon_distance_constraint)

        # Verify that there are no moon distance intervals
        assert_equal(len(moon_distance_intervals), 0)

    def test_moon_distance_major_planet_none_removed(self):
        start = datetime(2012, 2, 1)
        end = datetime(2012, 2, 2)
        v = Visibility(self.site, start, end, self.horizon)

        target = self.major_planet_target.copy()
        # According to JPL horizons, this target has a moon distance angle < 30.5 and > 18.5 all the time
        target_intervals = v.get_target_intervals(target=target)
        moon_distance_intervals = v.get_moon_distance_intervals(target, target_intervals, Angle(degrees=18))
        target_intervals = coalesce_adjacent_intervals(target_intervals)

        # Verify that moon distance intervals is the entire target interval
        assert_equal(len(moon_distance_intervals), len(target_intervals))
        assert_equal(moon_distance_intervals[0][0], target_intervals[0][0])
        assert_equal(moon_distance_intervals[0][1], target_intervals[0][1])

    def test_moon_distance_major_planet_all_removed(self):
        start = datetime(2012, 2, 1)
        end = datetime(2012, 2, 2)
        v = Visibility(self.site, start, end, self.horizon)

        target = self.major_planet_target.copy()
        # According to JPL horizons, this target has a moon distance angle < 30.5 all the time
        target_intervals = v.get_target_intervals(target=target)
        moon_distance_intervals = v.get_moon_distance_intervals(target, target_intervals, Angle(degrees=31))

        # Verify that there are no moon distance intervals
        assert_equal(len(moon_distance_intervals), 0)

    def test_moon_distance_major_planet_some_removed(self):
        start = datetime(2012, 2, 1)
        end = datetime(2012, 2, 2)
        v = Visibility(self.site, start, end, self.horizon)

        target = self.major_planet_target.copy()
        # According to JPL horizons, this target has a moon distance angle > 20.0 from 03:45 on
        target_intervals = v.get_target_intervals(target=target)
        moon_distance_intervals = v.get_moon_distance_intervals(target, target_intervals, Angle(degrees=20))
        target_intervals = coalesce_adjacent_intervals(target_intervals)

        # Verify that the first moon distance interval truncated to start at 3:45 but otherwise matches the target int
        assert_equal(len(moon_distance_intervals), 2)
        assert_equal(moon_distance_intervals[0][0], datetime(2012, 2, 1, 3, 45))
        assert_equal(moon_distance_intervals[0][1], target_intervals[0][1])
        assert_equal(moon_distance_intervals[1][0], target_intervals[1][0])
        assert_equal(moon_distance_intervals[1][1], target_intervals[1][1])


class TestZenithDistanceCalculation(object):
    """ All time intervals to test against are obtained from JPL Horizons
    """
    def setup(self):
        self.site = {
            'latitude': Angle(degrees=20.0),
            'longitude': Angle(degrees=-150.0),
            'ha_limit_neg': Angle(degrees=-4.6*12.0),
            'ha_limit_pos': Angle(degrees=4.6*12.0),
            'horizon': Angle(degrees=15.0)
        }

        self.horizon = 15.0
        self.sidereal_target = {
            'ra': RightAscension('20 41 25.91'),
            'dec': Declination('+20 00 00.00'),
            'epoch': 2000,
        }

        # Details from JPL Horizons for Ceres
        self.minor_planet_target = make_minor_planet_target('MPC_MINOR_PLANET',
                                                            epoch=49731,
                                                            inclination=10.60069567603618,
                                                            long_node=80.65851514365535,
                                                            arg_perihelion=71.44921526124109,
                                                            semi_axis=2.767218108003098,
                                                            eccentricity=0.07610292126891821,
                                                            mean_anomaly=340.389084821267)

        # Details from JPL Horizons for Comet 27P
        self.comet_target = make_comet_target('MPC_COMET',
                                              epoch=56364,
                                              epochofperih=55776.8910902,
                                              inclination=28.96687278723059,
                                              long_node=250.6264098390235,
                                              arg_perihelion=196.0253968913816,
                                              perihdist=0.748287467144728,
                                              eccentricity=0.9189810923126022)

        # Details from JPL Horizons for Jupiter
        self.major_planet_target = make_major_planet_target('JPL_MAJOR_PLANET',
                                                            epochofel=55959.0,
                                                            inclination=1.303884172546506,
                                                            long_node=100.5093329813755,
                                                            arg_perihelion=274.0516181838379,
                                                            semi_axis=5.204023536751508,
                                                            eccentricity=0.04910768996084790,
                                                            mean_anomaly=26.60707699766562,
                                                            dailymot=0.08306200006467207)

    def test_zenith_distance_no_angle(self):
        """
        test that that when the zenith distance is zero, no intervals are removed from the
        original target intervals
        """
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)
        zenith_distance_constraint = Angle(degrees=0)

        v = Visibility(self.site, start, end, self.horizon)
        target_intervals = v.get_target_intervals(target=self.sidereal_target)
        zenith_distance_intervals = v.get_zenith_distance_intervals(self.sidereal_target, target_intervals,
                                                                    zenith_distance_constraint)

        assert_equal(zenith_distance_intervals, target_intervals)

    def test_zenith_distance_all_removed_180(self):
        """
        test that with a zenith distance of 180 degrees, all intervals are removed
        """
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)
        # giant angle given, no target could further away so all interval should be disallowed
        zenith_distance_constraint = Angle(degrees=180)

        v = Visibility(self.site, start, end, self.horizon)
        target_intervals = v.get_target_intervals(target=self.sidereal_target)
        zenith_distance_intervals = v.get_zenith_distance_intervals(self.sidereal_target, target_intervals,
                                                                    zenith_distance_constraint)

        assert_equal(len(zenith_distance_intervals), 0)

    # TODO: finish implementation
    def test_zenith_distance_none_removed(self):
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)
        # low angle given, the distance is always greater for this target so it should allow all intervals
        zenith_distance_constraint = Angle(degrees=1.23)

        v = Visibility(self.site, start, end, self.horizon)
        target_intervals = v.get_target_intervals(target=self.sidereal_target)
        zenith_distance_intervals = v.get_zenith_distance_intervals(self.sidereal_target, target_intervals,
                                                                    zenith_distance_constraint)

        assert_equal(zenith_distance_intervals, target_intervals)

    # TODO: finish implementation
    def test_zenith_distance_below_horizon(self):
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)
        # low angle given, the distance is always greater for this target so it should allow all intervals
        zenith_distance_constraint = Angle(degrees=self.horizon)

        v = Visibility(self.site, start, end, self.horizon)
        target_intervals = v.get_target_intervals(target=self.sidereal_target)
        zenith_distance_intervals = v.get_zenith_distance_intervals(self.sidereal_target, target_intervals,
                                                                    zenith_distance_constraint)

        assert_equal(len(zenith_distance_intervals), 0)


