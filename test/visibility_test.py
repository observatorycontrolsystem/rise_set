#!/usr/bin/env python

from __future__ import division
from builtins import range
from builtins import object

import pytest
from datetime import datetime, timedelta

#Import the module to test
from rise_set.visibility import Visibility, set_airmass_limit, InvalidHourAngleLimit

# Additional support modules
from rise_set.angle import Angle
from rise_set.sky_coordinates import RightAscension, Declination
from rise_set.astrometry import (make_satellite_target, make_minor_planet_target, make_comet_target,
                                 make_major_planet_target, make_hour_angle_target, calculate_altitude,
                                 calculate_zenith_distance, mean_to_apparent, date_to_tdb,
                                 calc_local_hour_angle, elem_to_topocentric_apparent, target_to_jform)
from rise_set.rates import ProperMotion
from rise_set.moving_objects import initialise_sites
from rise_set.utils          import coalesce_adjacent_intervals, intersect_intervals, intersect_many_intervals
from mock import patch
import unittest

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


class TestIntervals(unittest.TestCase):

    def setup_method(self, method):

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
#        dark_intervals = visibility._find_when_target_is_up('sun', dt=datetime(2014, 3, 20))

        for start, end in dark_intervals:
            assert start < end


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
        received = self.visibility._find_when_target_is_up(self.sun, self.dt)


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
        received = self.visibility._find_when_target_is_down(self.sun, self.dt)

        # Ignore microseconds for these tests
        received = zero_out_microseconds(received)

        assert received == expected


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

        received = visibility._find_when_target_is_down(self.sun, dt, star=None,
                                                                     airmass=None)

        assert received[0][0] < received[0][1]



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

        assert received == expected


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

        assert received == expected

    @patch('rise_set.visibility.Visibility._get_ra_target_intervals')
    @patch('rise_set.visibility.Visibility._get_moving_object_target_intervals')
    def test_get_target_intervals_no_target_specified(self, mock_func1, mock_func2):
        target = {}

        received = self.visibility.get_target_intervals(target)

        assert mock_func1.call_count == 0
        assert mock_func2.call_count == 1

    def test_get_target_intervals_satellite_target(self):
        target = make_satellite_target(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        target_intervals = self.visibility.get_target_intervals(target)
        dark_intervals = self.visibility.get_dark_intervals()

        assert target_intervals == dark_intervals

    def test_get_target_intervals_hour_angle_target(self):
        target = make_hour_angle_target(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        target_intervals = self.visibility.get_target_intervals(target)
        dark_intervals = self.visibility.get_dark_intervals()

        assert target_intervals == dark_intervals

    @patch('rise_set.visibility.Visibility._get_ra_target_intervals')
    @patch('rise_set.visibility.Visibility._get_moving_object_target_intervals')
    def test_get_target_intervals_mpc_comet_type_is_a_moving_object(self, moving_obj_func,
                                                                    ra_dec_func):
        target = {'type' : 'mpc_comet'}

        received = self.visibility.get_target_intervals(target)

        assert moving_obj_func.call_count == 1
        assert ra_dec_func.call_count == 0

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

        assert mock_func.call_count == 1

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

        assert args[0] == expected_window
        assert args[1] == expected_target
        assert args[2] == expected_site

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

        assert v.get_observable_intervals(target, moon_distance=Angle(degrees=0)) == []

    def test_negative_hour_angle_too_big_is_rejected(self):
        with pytest.raises(InvalidHourAngleLimit):
            visibility = Visibility(None, None, None, ha_limit_neg=-13.0)

    def test_negative_hour_angle_too_small_is_rejected(self):
        with pytest.raises(InvalidHourAngleLimit):
            visibility = Visibility(None, None, None, ha_limit_neg=1.0)

    def test_positive_hour_angle_too_big_is_rejected(self):
        with pytest.raises(InvalidHourAngleLimit):
            visibility = Visibility(None, None, None, ha_limit_pos=13.0)

    def test_positive_hour_angle_too_small_is_rejected(self):
        with pytest.raises(InvalidHourAngleLimit):
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
        assert observable_intervals[0][0] == rise_time
        assert observable_intervals[0][1] == datetime(2012, 2, 1, 5, 45)

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

        assert observable_intervals[0][0] == rise_time
        assert observable_intervals[0][1] == datetime(2012, 2, 1, 7, 45)

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

        assert observable_intervals[0][0] == dark_intervals[0][0]
        assert observable_intervals[0][1] == dark_intervals[0][1]


class TestAirmassCalculation(unittest.TestCase):

    def setup_method(self, method):
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
        assert set_airmass_limit(None, horizon) == expected


    def test_set_airmass_limit_airmass_worse_than_horizon(self):
        horizon  = 30
        airmass  = 3
        expected = horizon
        assert set_airmass_limit(airmass, horizon) == expected


    def test_set_airmass_limit_airmass_better_than_horizon(self):
        horizon  = 30
        airmass  = 1.2
        expected = 56.44
        self.assertAlmostEqual(set_airmass_limit(airmass, horizon), expected, places=2)


    def test_airmass_is_applied_if_above_horizon_and_not_otherwise(self):
        airmass_above_horizon = 1.082392200292394   # 22.5 degrees from the zenith
        airmass_below_horizon = 3
        interval_with_airmass    = self.interval_for_airmass(airmass_above_horizon)
        interval_without_airmass = self.interval_for_airmass(airmass_below_horizon)

        assert interval_with_airmass < interval_without_airmass


class TestMoonDistanceCalculation(object):
    ''' All time intervals to test against are obtained from JPL Horizons
    '''
    def setup_method(self):
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

        assert moon_distance_intervals == target_intervals

    def test_moon_distance_none_removed(self):
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)
        # low angle given, the distance is always greater for this target so it should allow all intervals
        moon_distance_constraint = Angle(degrees=30)

        v = Visibility(self.site, start, end, self.horizon)
        target_intervals   = v.get_target_intervals(target=self.sidereal_target)
        moon_distance_intervals = v.get_moon_distance_intervals(self.sidereal_target, target_intervals, moon_distance_constraint)

        assert moon_distance_intervals == target_intervals

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
        assert moon_distance_intervals == [target_intervals[1]]

    def test_moon_distance_all_removed(self):
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)
        # this target/site/date has no moon distances less than 75 degrees
        moon_distance_constraint = Angle(degrees=75)

        v = Visibility(self.site, start, end, self.horizon)
        target_intervals   = v.get_target_intervals(target=self.sidereal_target)
        moon_distance_intervals = v.get_moon_distance_intervals(self.sidereal_target, target_intervals, moon_distance_constraint)

        # test that no moon distance intervals are returned due to the constraint of > 75 degrees distance
        assert moon_distance_intervals == []

    def test_moon_distance_ignored_for_satellite_target(self):
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)
        v = Visibility(self.site, start, end, self.horizon)

        target = make_satellite_target(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        observable_intervals = v.get_observable_intervals(target)
        # check that this doesn't crash, and that intervals match night intervals
        night_intervals = v.get_dark_intervals()
        assert night_intervals == observable_intervals

    def test_moon_distance_ignored_for_hour_angle_target(self):
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)
        v = Visibility(self.site, start, end, self.horizon)

        target = make_hour_angle_target(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        observable_intervals = v.get_observable_intervals(target)
        # check that this doesn't crash, and that intervals match night intervals
        night_intervals = v.get_dark_intervals()
        assert night_intervals == observable_intervals

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
        assert moon_distance_intervals[0][0] == target_intervals[0][0]
        assert moon_distance_intervals[0][1] == datetime(2012, 8, 2, 21, 0)
        assert len(moon_distance_intervals) == 1

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
        assert len(moon_distance_intervals) == 2
        assert moon_distance_intervals[0][0] == datetime(2012, 1, 2, 6, 15)
        assert moon_distance_intervals[0][1] == datetime(2012, 1, 2, 6, 45)
        assert moon_distance_intervals[1][0] == datetime(2012, 1, 2, 23, 30)
        assert moon_distance_intervals[1][1] == datetime(2012, 1, 3, 0, 0)

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
        assert len(moon_distance_intervals) == 2
        assert moon_distance_intervals[0][0] == datetime(2012, 2, 2, 0, 0)
        assert moon_distance_intervals[0][1] == datetime(2012, 2, 2, 5, 15)
        assert moon_distance_intervals[1][0] == datetime(2012, 2, 2, 22, 00)
        assert moon_distance_intervals[1][1] == datetime(2012, 2, 3, 0, 0)

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
        assert len(moon_distance_intervals) == 0

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

        # Verify that moon distance intervals is the entire target interval
        assert len(moon_distance_intervals) == 1
        assert moon_distance_intervals[0][0] == target_intervals[0][0]
        assert moon_distance_intervals[0][1] == target_intervals[-1][1]

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
        assert len(moon_distance_intervals) == 0

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
        assert len(moon_distance_intervals) == len(target_intervals)
        assert moon_distance_intervals[0][0] == target_intervals[0][0]
        assert moon_distance_intervals[0][1] == target_intervals[0][1]

    def test_moon_distance_major_planet_all_removed(self):
        start = datetime(2012, 2, 1)
        end = datetime(2012, 2, 2)
        v = Visibility(self.site, start, end, self.horizon)

        target = self.major_planet_target.copy()
        # According to JPL horizons, this target has a moon distance angle < 30.5 all the time
        target_intervals = v.get_target_intervals(target=target)
        moon_distance_intervals = v.get_moon_distance_intervals(target, target_intervals, Angle(degrees=31))

        # Verify that there are no moon distance intervals
        assert len(moon_distance_intervals) == 0

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
        assert len(moon_distance_intervals) == 2
        assert moon_distance_intervals[0][0] == datetime(2012, 2, 1, 3, 45)
        assert moon_distance_intervals[0][1] == target_intervals[0][1]
        assert moon_distance_intervals[1][0] == target_intervals[1][0]
        assert moon_distance_intervals[1][1] == target_intervals[1][1]


class TestMoonPhaseCalculation(object):
    ''' All time intervals to test against are obtained from JPL Horizons
    '''
    def setup_method(self):
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

    def test_moon_phase_none_removed(self):
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)
        max_moon_phase = 1.0

        v = Visibility(self.site, start, end, self.horizon)
        target_intervals   = v.get_target_intervals(target=self.sidereal_target)
        moon_phase_intervals = v.get_moon_phase_intervals(target_intervals, max_moon_phase)

        assert moon_phase_intervals == target_intervals

        # Test the observable intervals as well, must combine moon phase with dark and ha intervals
        observable_intervals = v.get_observable_intervals(target=self.sidereal_target, moon_distance=Angle(degrees=0), moon_phase=max_moon_phase)
        ha_intervals = v.get_ha_intervals(self.sidereal_target)
        dark_intervals = v.get_dark_intervals()
        combined_moon_phase_intervals = intersect_many_intervals(moon_phase_intervals, ha_intervals, dark_intervals)
        assert combined_moon_phase_intervals == observable_intervals

    def test_moon_phase_all_removed_except_when_moon_down(self):
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)
        # Moon phase is always above 0.5 in this time range, so only moon dark periods will go through
        max_moon_phase = 0.5

        v = Visibility(self.site, start, end, self.horizon)
        target_intervals   = v.get_target_intervals(target=self.sidereal_target)
        moon_phase_intervals = v.get_moon_phase_intervals(target_intervals, max_moon_phase)
        moon_down_intervals = v.get_moon_dark_intervals()
        moon_down_intervals = intersect_intervals(target_intervals, moon_down_intervals)

        assert moon_phase_intervals == moon_down_intervals

        # Test the observable intervals as well, must combine moon phase with dark and ha intervals
        observable_intervals = v.get_observable_intervals(target=self.sidereal_target, moon_distance=Angle(degrees=0), moon_phase=max_moon_phase)
        ha_intervals = v.get_ha_intervals(self.sidereal_target)
        dark_intervals = v.get_dark_intervals()
        combined_moon_phase_intervals = intersect_many_intervals(moon_phase_intervals, ha_intervals, dark_intervals)
        assert combined_moon_phase_intervals == observable_intervals

    def test_moon_phase_some_removed(self):
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)
        # Moon phase goes above 0.6 at ~09:30
        max_moon_phase = 0.6

        v = Visibility(self.site, start, end, self.horizon)
        target_intervals   = v.get_target_intervals(target=self.sidereal_target)
        moon_phase_intervals = v.get_moon_phase_intervals(target_intervals, max_moon_phase)
        moon_down_intervals = v.get_moon_dark_intervals()
        # The first target interval is below the moon phase constraint, but the second is above it so its only good while the moon is down
        expected_intervals = [target_intervals[0], (target_intervals[1][0], moon_down_intervals[0][1])]
        assert moon_phase_intervals == expected_intervals

        # Test the observable intervals as well, must combine moon phase with dark and ha intervals
        observable_intervals = v.get_observable_intervals(target=self.sidereal_target, moon_distance=Angle(degrees=0), moon_phase=max_moon_phase)
        ha_intervals = v.get_ha_intervals(self.sidereal_target)
        dark_intervals = v.get_dark_intervals()
        combined_moon_phase_intervals = intersect_many_intervals(moon_phase_intervals, ha_intervals, dark_intervals)
        assert combined_moon_phase_intervals == observable_intervals

    def test_moon_phase_longer(self):
        start = datetime(2012, 1, 1)
        end = datetime(2012, 2, 1)
        # Somewhat low max moon phase, should only allow intervals between 1/18 - 1/29
        max_moon_phase = 0.3

        v = Visibility(self.site, start, end, self.horizon)
        target_intervals   = v.get_target_intervals(target=self.sidereal_target)
        moon_phase_intervals = v.get_moon_phase_intervals(target_intervals, max_moon_phase)
        moon_down_intervals = v.get_moon_dark_intervals()
        moon_down_intervals = intersect_intervals(target_intervals, moon_down_intervals)

        expected_intervals = []
        # Add in dark moon intervals while the constraint is not met
        for interval in moon_down_intervals:
            if interval[0] < datetime(2012, 1, 18, 5):
                expected_intervals.append(interval)

        # Add in target intervals while the constraint is met
        for interval in target_intervals:
            if interval[0] > datetime(2012, 1, 18, 5) and interval[0] < datetime(2012, 1, 28, 16):
                expected_intervals.append(interval)

        # Add in one interval to bridge the gap - this one ends at 1-28-22-57 is about when moon phase goes above the constraint again.
        expected_intervals.append((datetime(2012, 1, 28, 16, 49, 18, 193648), datetime(2012, 1, 28, 22, 57, 23, 36345)))
        # Add dark moon intervals when the moon phase constraint is no longer met
        for interval in moon_down_intervals:
            if interval[0] > datetime(2012, 1, 28, 23):
                expected_intervals.append(interval)

        assert moon_phase_intervals == expected_intervals


class TestZenithDistanceCalculation(unittest.TestCase):
    """ All time intervals to test against are obtained from JPL Horizons
    """
    def setup_method(self, method):
        self.site = {
            'latitude': Angle(degrees=20.0),
            'longitude': Angle(degrees=-150.0),
            # TODO: sort out if there's a problem with ha_limits in Visibility.get_target_intervals
            #'ha_limit_neg': Angle(degrees=-4.6*12.0),
            #'ha_limit_pos': Angle(degrees=4.6*12.0),
            'ha_limit_neg': Angle(degrees=-36*12.0),
            'ha_limit_pos': Angle(degrees=36*12.0),
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

        # for datetime(2012, 07, 22) to datetime(2012, 07, 23)
        # Details from JPL Horizons for Comet 27P
        self.comet_target = make_comet_target('MPC_COMET',
                                              epoch=56364,
                                              epochofperih=55776.8910902,
                                              inclination=28.96687278723059,
                                              long_node=250.6264098390235,
                                              arg_perihelion=196.0253968913816,
                                              perihdist=0.748287467144728,
                                              eccentricity=0.9189810923126022)

        # 1. this Epoch of Elements (epochofel) is for a modified Julian date (MJD) corresponding
        #    to datetime(2012, 2, 2)
        # 2. datetime(2012, 2, 1) is the start time used in tests using self.major_planet_target
        # 3. datetime(2012, 2, 2) is the end time used in tests using self.major_planet_target
        # 4. JD = MDJ + 2400000.5
        # 5. see https://aa.usno.navy.mil/jdconverter?ID=AA&jd=2455959.5
        #
        # Orbital Elements from JPL Horizons for Jupiter
        self.major_planet_target = make_major_planet_target('JPL_MAJOR_PLANET',
                                                            epochofel=55959.0,                 # date at top
                                                            inclination=1.303884172546506,     # IN
                                                            long_node=100.5093329813755,       # OM
                                                            arg_perihelion=274.0516181838379,  # W
                                                            semi_axis=5.204023536751508,       # A
                                                            eccentricity=0.04910768996084790,  # EC
                                                            mean_anomaly=26.60707699766562,    # MA
                                                            dailymot=0.08306200006467207)      # N


class TestObservableIntervalsZDIgnoredForStaticTargets(TestZenithDistanceCalculation):
    """
    The zenith distance calculation should be ignored for satellite and hour angle targets.
    """
    def setup_method(self, method):
        super(TestObservableIntervalsZDIgnoredForStaticTargets, self).setup_method(method)
        start = datetime(2012, 2, 1)
        end = datetime(2012, 2, 2)

        # giant zenith_blind_spot would normally limit intervals
        self.v = Visibility(self.site, start, end, self.horizon, zenith_blind_spot=180)

    def test_moon_distance_ignored_for_satellite_target(self):
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)

        target = make_satellite_target(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        observable_intervals = self.v.get_observable_intervals(target)
        # check that this doesn't crash, and that intervals match night intervals
        night_intervals = self.v.get_dark_intervals()
        assert night_intervals == observable_intervals

    def test_moon_distance_ignored_for_hour_angle_target(self):
        start = datetime(2012, 1, 2)
        end = datetime(2012, 1, 3)

        target = make_hour_angle_target(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        observable_intervals = self.v.get_observable_intervals(target)
        # check that this doesn't crash, and that intervals match night intervals
        night_intervals = self.v.get_dark_intervals()
        assert night_intervals == observable_intervals


class TestZDIntervalsZeroZenithDistance(TestZenithDistanceCalculation):
    """Test that a zero zenith distance leaves the target intervals untouched.
    """

    def setup_method(self, method):
        super(TestZDIntervalsZeroZenithDistance, self).setup_method(method)
        start = datetime(2012, 2, 1)
        end = datetime(2012, 2, 2)
        self.v = Visibility(self.site, start, end, self.horizon)  # zenith_blind_spot defaults to 0

    def test_zd_intervals_zero_zd_sidereal_target(self):
        target = self.sidereal_target

        target_intervals = self.v.get_target_intervals(target=target)
        zenith_distance_intervals = self.v.get_zenith_distance_intervals(target, target_intervals)

        coalesced_target_intervals = coalesce_adjacent_intervals(target_intervals)
        assert zenith_distance_intervals == coalesced_target_intervals

    def test_zd_intervals_zero_zd_non_sidereal_major_planet_target(self):
        target = self.major_planet_target

        target_intervals = self.v.get_target_intervals(target=target)
        zenith_distance_intervals = self.v.get_zenith_distance_intervals(target, target_intervals)

        coalesced_target_intervals = coalesce_adjacent_intervals(target_intervals)
        assert zenith_distance_intervals == coalesced_target_intervals

    def test_zd_intervals_zero_zd_non_sidereal_minor_planet_target(self):
        target = self.minor_planet_target

        target_intervals = self.v.get_target_intervals(target=target)
        zenith_distance_intervals = self.v.get_zenith_distance_intervals(target, target_intervals)

        coalesced_target_intervals = coalesce_adjacent_intervals(target_intervals)
        assert zenith_distance_intervals == coalesced_target_intervals

    def test_zd_intervals_zero_zd_non_sidereal_comet_target(self):
        target = self.comet_target

        target_intervals = self.v.get_target_intervals(target=target)
        zenith_distance_intervals = self.v.get_zenith_distance_intervals(target, target_intervals)

        coalesced_target_intervals = coalesce_adjacent_intervals(target_intervals)
        assert zenith_distance_intervals == coalesced_target_intervals


class TestZDIntervals180ZenithDistance(TestZenithDistanceCalculation):
    """Test that a giant zenith distance like 180 removes all intervals from a variety of targets.
    """

    def setup_method(self, method):
        super(TestZDIntervals180ZenithDistance, self).setup_method(method)
        start = datetime(2012, 2, 1)
        end = datetime(2012, 2, 2)
        self.v = Visibility(self.site, start, end, self.horizon, zenith_blind_spot=180.0)

    def test_zd_intervals_180_zd_sidereal_target(self):
        target = self.sidereal_target

        target_intervals = self.v.get_target_intervals(target=target)
        zenith_distance_intervals = self.v.get_zenith_distance_intervals(target, target_intervals)
        assert len(zenith_distance_intervals) == 0

    def test_zd_intervals_180_zd_non_sidereal_major_planet_target(self):
        target = self.major_planet_target

        target_intervals = self.v.get_target_intervals(target=target)
        zenith_distance_intervals = self.v.get_zenith_distance_intervals(target, target_intervals)

        assert len(zenith_distance_intervals) == 0

    def test_zd_intervals_180_zd_non_sidereal_minor_planet_target(self):
        target = self.minor_planet_target

        target_intervals = self.v.get_target_intervals(target=target)
        zenith_distance_intervals = self.v.get_zenith_distance_intervals(target, target_intervals)

        assert len(zenith_distance_intervals) == 0

    def test_zd_intervals_180_zd_non_sidereal_comet_target(self):
        target = self.comet_target

        target_intervals = self.v.get_target_intervals(target=target)
        zenith_distance_intervals = self.v.get_zenith_distance_intervals(target, target_intervals)

        assert len(zenith_distance_intervals) == 0


class TestMajorPlanetZenithDistanceIntervals(TestZenithDistanceCalculation):

    def setup_method(self, method):
        super(TestMajorPlanetZenithDistanceIntervals, self).setup_method(method)
        start = datetime(2012, 2, 1)
        end = datetime(2012, 2, 2)
        zenith_hole_radius = 9.0  # degrees; altitude > 81 should be excluded
        self.v = Visibility(self.site, start, end, self.horizon,
                            ha_limit_neg=-12.0, ha_limit_pos=12.0,
                            zenith_blind_spot=zenith_hole_radius)

    def test_jupiter(self):
        """Test interval removal for Jupiter > 81-degrees
        Excerpts from JPL Horizons, between the start and end datetimes given, Jupiter is over 81-degrees

        |-------------------+----------|
        | time              | altitude |
        |-------------------+----------|
        | 2012-Feb-01 03:09 |  80.9600 |
        | 2012-Feb-01 03:10 |  81.0298 |
        | 2012-Feb-01 03:11 |  81.0936 | enter zenith hole
        | 2012-Feb-01 03:12 |  81.1514 |
        | <snip>            |          |
        | 2012-Feb-01 03:30 |  81.0719 |
        | 2012-Feb-01 03:31 |  81.0060 |
        | 2012-Feb-01 03:32 |  80.9341 | exit zenith hole
        | 2012-Feb-01 03:33 |  80.8564 |
        |-------------------+----------|
        """
        target = self.major_planet_target  # Jupiter

        target_intervals = self.v.get_target_intervals(target=target)
        zenith_intervals = self.v.get_zenith_distance_intervals(target, target_intervals,
                                                                chunk_size=timedelta(minutes=1))

        # these intervals are tailored to the specific JPL Horizons output for this target
        expected_intervals = [
            (self.v.start_date, datetime(2012, 2, 1, 3, 10)),
            (datetime(2012, 2, 1, 3, 32), datetime(2012, 2, 1, 8, 30)),
            (datetime(2012, 2, 1, 22, 15), self.v.end_date)
        ]
        assert expected_intervals == zenith_intervals


class TestCometZenithDistanceIntervals(TestZenithDistanceCalculation):
    """
    Excerpts from JPL Horizons
    |-------------------+----------+-------------------|
    | time              | altitude | event             |
    |-------------------+----------+-------------------|
    | 2012-Jul-22 03:08 |  43.9532 |                   |
    | 2012-Jul-22 03:09 |  44.0039 | enter zenith hole |
    | 2012-Jul-22 03:10 |  44.0533 |                   |
    | 2012-Jul-22 04:25 |  44.0802 |                   |
    | 2012-Jul-22 04:26 |  44.0314 |                   |
    | 2012-Jul-22 04:27 |  43.9814 | exit zenith hole  |
    | 2012-Jul-22 04:28 |  43.9302 |                   |
    | 2012-Jul-22 04:29 |  43.8777 |                   |
    | 2012-Jul-22 07:54 |  15.0088 | comet set         |
    | 2012-Jul-22 07:55 |  14.8141 |                   |
    | 2012-Jul-22 23:37 |  14.8374 |                   |
    | 2012-Jul-22 23:38 |  15.0322 | comet rise        |
    |-------------------+----------+-------------------|
    """
    def setup_method(self, method):
        super(TestCometZenithDistanceIntervals, self).setup_method(method)
        start = datetime(2012, 7, 22)
        end = datetime(2012, 7, 23)
        zenith_hole_radius = 46.0  # degrees; altitude > 44 should be excluded
        self.v = Visibility(self.site, start, end, self.horizon,
                            ha_limit_neg=-12.0, ha_limit_pos=12.0,
                            zenith_blind_spot=zenith_hole_radius)

    def test_sidereal_target(self):
        target = self.comet_target

        target_intervals = self.v.get_target_intervals(target=target)
        zenith_intervals = self.v.get_zenith_distance_intervals(target, target_intervals,
                                                                chunk_size=timedelta(minutes=1))
        # first interval ends when the target enters the zenith hole (chunk_size = 1 min)
        # second interval begins when the target exits the zenith hole
        # second interval ends when the target goes below the horizon (15) (chunk_size=15 min default)
        # third interval begins when the target rises above the horizon (15) (chunk_size=15 min default)
        expected_intervals = [
            (self.v.start_date, datetime(2012, 7, 22, 3, 9)),
            (datetime(2012, 7, 22, 4, 27), datetime(2012, 7, 22, 7, 45)),
            (datetime(2012, 7, 22, 23, 45), self.v.end_date)
        ]
        assert expected_intervals == zenith_intervals


class TestZDvsAltitude(TestZenithDistanceCalculation):
    """Test that zd + alt = 90-degrees for a variety of targets.
    """

    def test_zd_vs_alt_degenerate(self):
        """
        degenerate case: test that zd + alt = 90-degrees
        """
        latitude = Angle(0.0)
        dec = Angle(0.0)
        local_hour_angle = Angle(0.0)
        zd = calculate_zenith_distance(latitude, dec, local_hour_angle)
        alt = calculate_altitude(latitude.in_degrees(),
                                 dec.in_degrees(),
                                 local_hour_angle.in_degrees())
        self.assertAlmostEqual((90.0 - zd.in_degrees()), alt.in_degrees())

    def test_zd_vs_alt_sidereal_target(self):
        """
        sidereal target: test that zd + alt = 90-degrees
        """
        target = self.sidereal_target
        start_time = datetime(2012, 2, 1)
        tdb = date_to_tdb(start_time)
        # get the apparent ra/dec for the target
        target_app_ra, target_app_dec = mean_to_apparent(target, tdb)  # for sidereal targets
        # get the local_hour_angle from
        latitude = self.site['latitude']

        local_hour_angle = calc_local_hour_angle(target_app_ra,
                                                 self.site['longitude'], start_time)

        zd = calculate_zenith_distance(latitude, target_app_dec, local_hour_angle)

        alt = calculate_altitude(latitude.in_degrees(),
                                 target_app_dec.in_degrees(),
                                 local_hour_angle.in_degrees())
        self.assertAlmostEqual((90.0 - zd.in_degrees()), alt.in_degrees(), places=4)

    def test_zd_vs_alt_non_sidereal_major_planet_target(self):
        """
        non-sidereal (major_planet): test that zd + alt = 90-degrees
        """
        target = self.major_planet_target
        start_time = datetime(2012, 2, 1)
        # get the apparent ra/dec for the target
        # target_app_ra, target_app_dec = mean_to_apparent(target, tdb)  # for sidereal targets
        # for non-sidereal targets
        target_app_ra, target_app_dec = elem_to_topocentric_apparent(start_time, target, self.site,
                                                                     target_to_jform(target))

        latitude = self.site['latitude']
        local_hour_angle = calc_local_hour_angle(target_app_ra, self.site['longitude'], start_time)
        zd = calculate_zenith_distance(latitude, target_app_dec, local_hour_angle)
        alt = calculate_altitude(latitude.in_degrees(),
                                 target_app_dec.in_degrees(),
                                 local_hour_angle.in_degrees())
        self.assertAlmostEqual((90.0 - zd.in_degrees()), alt.in_degrees(), places=4)

    def test_zd_vs_alt_non_sidereal_minor_planet_target(self):
        """
        non-sidereal (minor_planet): test that zd + alt = 90-degrees
        """
        target = self.minor_planet_target
        start_time = datetime(2012, 2, 1)
        # get the apparent ra/dec for the target
        # target_app_ra, target_app_dec = mean_to_apparent(target, tdb)  # for sidereal targets
        # for non-sidereal targets
        target_app_ra, target_app_dec = elem_to_topocentric_apparent(start_time, target, self.site,
                                                                     target_to_jform(target))

        latitude = self.site['latitude']
        local_hour_angle = calc_local_hour_angle(target_app_ra, self.site['longitude'], start_time)
        zd = calculate_zenith_distance(latitude, target_app_dec, local_hour_angle)
        alt = calculate_altitude(latitude.in_degrees(),
                                 target_app_dec.in_degrees(),
                                 local_hour_angle.in_degrees())
        self.assertAlmostEqual((90.0 - zd.in_degrees()), alt.in_degrees(), places=4)

    def test_zd_vs_alt_non_sidereal_comet_target(self):
        """
        non-sidereal (comet): test that zd + alt = 90-degrees
        """
        target = self.comet_target
        start_time = datetime(2012, 2, 1)
        # get the apparent ra/dec for the target
        # target_app_ra, target_app_dec = mean_to_apparent(target, tdb)  # for sidereal targets
        # for non-sidereal targets
        target_app_ra, target_app_dec = elem_to_topocentric_apparent(start_time, target, self.site,
                                                                     target_to_jform(target))

        latitude = self.site['latitude']
        local_hour_angle = calc_local_hour_angle(target_app_ra, self.site['longitude'], start_time)
        zd = calculate_zenith_distance(latitude, target_app_dec, local_hour_angle)
        alt = calculate_altitude(latitude.in_degrees(),
                                 target_app_dec.in_degrees(),
                                 local_hour_angle.in_degrees())
        self.assertAlmostEqual((90.0 - zd.in_degrees()), alt.in_degrees(), places=4)
