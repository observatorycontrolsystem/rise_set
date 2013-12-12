#!/usr/bin/env python

'''
moving_objects_test.py - Tests for moving object support

Author: Eric Saunders
December 2013
'''

from rise_set.angle import Angle
from rise_set.moving_objects import (Sites, initialise_sites,
                                     chunk_windows,
                                     elem_to_topocentric_apparent,
                                     calc_ephemerides,
                                     find_moving_object_up_intervals,
                                     find_moving_object_network_up_intervals)

from datetime import datetime, timedelta
from nose.tools import assert_equal, assert_almost_equal, nottest


def str_to_dt(dt_str):
   '''Simple helper to make the tests more readable.'''
   return datetime.strptime(dt_str, '%Y-%m-%d %H:%M:%S')


class TestSites(object):

    def setup(self):
        self.site_dict1 = {
                            'name'      : '1m0a.doma.cpt',
                            'latitude'  : -32.38059,
                            'longitude' :  20.8101083333,
                            'horizon'   : 30,
                          }
        self.site_dict2 = {
                            'name'      : '1m0a.doma.lsc',
                            'latitude'  : -30.1673472222,
                            'longitude' : -70.8046722222,
                            'horizon'   : 30,
                          }
        self.site_dict3 = {
                            'name'      : '1m0b.domb.cpt',
                            'latitude'  : -32.38059,
                            'longitude' :  20.8101083333,
                            'horizon'   : 30,
                          }

        self.sites = Sites()


    def test_can_add_new_site(self):
        expected = {
                     '1m0a.doma.cpt' : self.site_dict1,
                     '1m0a.doma.lsc' : self.site_dict2,
                   }

        self.sites.add_site(self.site_dict1)
        self.sites.add_site(self.site_dict2)

        assert_equal(self.sites.sites, expected)


    def test_can_add_new_sites(self):
        expected = {
                     '1m0a.doma.cpt' : self.site_dict1,
                     '1m0a.doma.lsc' : self.site_dict2,
                   }

        self.sites.add_sites([self.site_dict1, self.site_dict2])

        assert_equal(self.sites.sites, expected)


    def test_cant_add_same_site_coords_again(self):
        expected = {
                     '1m0a.doma.cpt' : self.site_dict1,
                   }

        self.sites.add_site(self.site_dict1)
        self.sites.add_site(self.site_dict3)

        assert_equal(self.sites.sites, expected)


    def test_can_get_site(self):
        self.sites.add_site(self.site_dict1)
        received = self.sites.get(self.site_dict1['name'])

        assert_equal(received, self.site_dict1)


    def test_none_if_no_site_found(self):
        received = self.sites.get(self.site_dict1['name'])

        assert_equal(received, None)


    def test_initialise_sites(self):
        expected_keys = ['1m0a.doma.elp', '1m0a.doma.coj',
                         '1m0a.doma.cpt', '1m0a.doma.lsc']
        sites = initialise_sites('test/telescopes.dat')

        assert_equal(len(sites.sites), 4)
        assert_equal(sites.sites.keys(), expected_keys)



class TestMovingObjects(object):

    def setup(self):
        # Asteroid (470) Kilia
        # From MPC webpages: http://minorplanetcenter.net/iau/MPEph/MPEph.html
        # Epoch 2013 Nov. 4.0 TT = JDT 2456600.5                  MPC
        # M  54.47380              (2000.0)            P               Q
        # n   0.26424476     Peri.   47.60658     -0.75565416     +0.65480399             T = 2456394.35096 JDT
        # a   2.4050925      Node   173.25052     -0.63179958     -0.72278517             q =     2.1781735
        # e   0.0943494      Incl.    7.22565     -0.17267333     -0.22093739
        # P   3.73           H   10.07          G   0.15           U   0
        # From 1313 observations at 37 oppositions, 1901-2013, mean residual 0".53.

        self.elements = {
                          'epoch'          : datetime(2013, 11, 4),
                          'inclination'    : Angle(degrees=7.22565),
                          'long_node'      : Angle(degrees=173.25052),
                          'arg_perihelion' : Angle(degrees=47.60658),
                          'semi_axis'      : 2.4050925,
                          'eccentricity'   : 0.0943494,
                          'mean_anomaly'   : Angle(degrees=54.47380),
                        }

        self.elp = {
                     'latitude'  : Angle(degrees=30.6801),
                     'longitude' : Angle(degrees=-104.015194444),
                   }


    def test_elem_to_topocentric_apparent_time1(self):
        dt = datetime(2013, 12, 9)
        tdb = dt - timedelta(seconds=67.184)

        ra, dec = elem_to_topocentric_apparent(tdb, self.elements, self.elp)


        # Coordinates of (470) Kilia from JPL Horizons: 284.70994 -18.19420
        assert_almost_equal(ra.in_degrees(), 284.70994, places=2)
        assert_almost_equal(dec.in_degrees(), -18.19420, places=4)


    def test_elem_to_topocentric_apparent_time2(self):
        dt = datetime(2013, 6, 9)
        tdb = dt - timedelta(seconds=67.184)

        ra, dec = elem_to_topocentric_apparent(tdb, self.elements, self.elp)


        # Coordinates from JPL Horizons: 225.22256  -5.29589
        assert_almost_equal(ra.in_degrees(), 225.22256, places=2)
        assert_almost_equal(dec.in_degrees(), -5.29589, places=2)


    def test_calc_ephemerides(self):
        window = {
                   'start' : datetime(2013, 12, 10),
                   'end'   : datetime(2013, 12, 11)
                 }
        site_dict1 = {
                        'name'      : '1m0a.doma.cpt',
                        'latitude'  : Angle(degrees=-32.38059),
                        'longitude' : Angle(degrees=20.8101083333),
                      }

        chunksize = timedelta(hours=6)

        expected = [
                     {
                       'ra_app'  : Angle(degrees=285.22681),
                       'dec_app' : Angle(degrees=-18.16350),
                       'start'   : window['start'],
                       'end'     : window['start'] + (1 * chunksize),
                     },
                     {
                       'ra_app'  : Angle(degrees=285.35677),
                       'dec_app' : Angle(degrees=-18.15597),
                       'start'   : window['start'] + (1 * chunksize),
                       'end'     : window['start'] + (2 * chunksize),
                     },
                     {
                       'ra_app'  : Angle(degrees=285.48538),
                       'dec_app' : Angle(degrees=-18.14841),
                       'start'   : window['start'] + (2 * chunksize),
                       'end'     : window['start'] + (3 * chunksize),
                     },
                     {
                       'ra_app'  : Angle(degrees=285.61369),
                       'dec_app' : Angle(degrees=-18.14033),
                       'start'   : window['start'] + (3 * chunksize),
                       'end'     : window['start'] + (4 * chunksize),
                     },
                     {
                       'ra_app'  : Angle(degrees=285.74337),
                       'dec_app' : Angle(degrees=-18.13206),
                       'start'   : window['start'] + (4 * chunksize),
                       'end'     : window['start'] + (4 * chunksize),
                     },
                   ]

        received = calc_ephemerides(window, self.elements, site_dict1, chunksize)

        for e, r in zip(expected, received):
            assert_almost_equal(e['ra_app'].in_degrees(), r['ra_app'].in_degrees(), places=3)
            assert_almost_equal(e['dec_app'].in_degrees(), r['dec_app'].in_degrees(), places=5)
            assert_equal(e['start'], r['start'])
            assert_equal(e['end'], r['end'])


    def test_find_moving_object_up_intervals1(self):
        window = {
                   'start' : datetime(2013, 12, 10),
                   'end'   : datetime(2013, 12, 11)
                 }
        site_dict1 = {
                        'name'      : '1m0a.doma.cpt',
                        'latitude'  : Angle(degrees=-32.38059),
                        'longitude' : Angle(degrees=20.8101083333),
                        'horizon'   : Angle(degrees=30),
                      }

        chunksize = timedelta(hours=3)

        expected_intervals = [
                                (
                                  datetime(2013, 12, 10, 9),
                                  datetime(2013, 12, 10, 12)
                                ),
                                (
                                  datetime(2013, 12, 10, 12),
                                  datetime(2013, 12, 10, 15)
                                )
                             ]
        expected_altitudes = [
                               Angle(degrees=42.7537),
                               Angle(degrees=74.9902)
                             ]

        received_ints, received_alts = find_moving_object_up_intervals(window,
                                                               self.elements,
                                                               site_dict1, chunksize)

        assert_equal(len(expected_intervals), len(received_ints))
        assert_equal(len(expected_altitudes), len(received_alts))

        for e, r in zip(expected_intervals, received_ints):
            assert_equal(e, r)

        for e, r in zip(expected_altitudes, received_alts):
            assert_almost_equal(e.in_degrees(), r.in_degrees(), places=3)


    def test_find_moving_object_up_intervals2(self):
        window = {
                   'start' : datetime(2013, 12, 10, 7, 30),
                   'end'   : datetime(2013, 12, 10, 8, 30)
                 }
        site_dict1 = {
                        'name'      : '1m0a.doma.cpt',
                        'latitude'  : Angle(degrees=-32.38059),
                        'longitude' : Angle(degrees=20.8101083333),
                        'horizon'   : Angle(degrees=30),
                      }

        chunksize = timedelta(minutes=15)

        expected_intervals = [
                               (
                                 datetime(2013, 12, 10, 8, 0),
                                 datetime(2013, 12, 10, 8, 15)
                               ),
                               (
                                 datetime(2013, 12, 10, 8, 15),
                                 datetime(2013, 12, 10, 8, 30)
                               ),
                             ]

        expected_altitudes = [
                               Angle(degrees=30.0814),
                               Angle(degrees=33.2494),
                             ]

        received_ints, received_alts = find_moving_object_up_intervals(window,
                                                               self.elements,
                                                               site_dict1, chunksize)

        assert_equal(len(expected_intervals), len(received_ints))
        assert_equal(len(expected_altitudes), len(received_alts))

        for e, r in zip(expected_intervals, received_ints):
            assert_equal(e, r)

        for e, r in zip(expected_altitudes, received_alts):
            assert_almost_equal(e.in_degrees(), r.in_degrees(), places=3)


    def test_find_moving_object_up_intervals3(self):
        window = {
                   'start' : datetime(2013, 12, 10, 7, 30),
                   'end'   : datetime(2013, 12, 10, 8, 31)
                 }
        site_dict1 = {
                        'name'      : '1m0a.doma.cpt',
                        'latitude'  : Angle(degrees=-32.38059),
                        'longitude' : Angle(degrees=20.8101083333),
                        'horizon'   : Angle(degrees=30),
                      }

        chunksize = timedelta(minutes=15)

        expected_intervals = [
                               (
                                 datetime(2013, 12, 10, 8, 0),
                                 datetime(2013, 12, 10, 8, 15)
                               ),
                               (
                                 datetime(2013, 12, 10, 8, 15),
                                 datetime(2013, 12, 10, 8, 30)
                               ),
                               (
                                 datetime(2013, 12, 10, 8, 30),
                                 datetime(2013, 12, 10, 8, 31)
                               ),
                             ]

        expected_altitudes = [
                               Angle(degrees=30.0814),
                               Angle(degrees=33.2494),
                               Angle(degrees=36.4200)
                             ]

        received_ints, received_alts = find_moving_object_up_intervals(window,
                                                               self.elements,
                                                               site_dict1, chunksize)

        assert_equal(len(expected_intervals), len(received_ints))
        assert_equal(len(expected_altitudes), len(received_alts))

        for e, r in zip(expected_intervals, received_ints):
            assert_equal(e, r)

        for e, r in zip(expected_altitudes, received_alts):
            assert_almost_equal(e.in_degrees(), r.in_degrees(), places=3)


    def test_find_moving_object_network_up_intervals(self):
        window = {
                   'start' : datetime(2013, 12, 10, 7, 30),
                   'end'   : datetime(2013, 12, 10, 8, 30)
                 }

        site_filename = 'test/telescopes.dat'

        chunksize = timedelta(minutes=15)

        expected = {
                     '1m0a.doma.elp' : [],
                     '1m0a.doma.coj' : [
                                         (
                                           datetime(2013, 12, 10, 7, 30),
                                           datetime(2013, 12, 10, 8, 0)
                                         )
                                       ],
                     '1m0a.doma.cpt' : [
                                         (
                                           datetime(2013, 12, 10, 8, 0),
                                           datetime(2013, 12, 10, 8, 30)
                                         )
                                       ],
                     '1m0a.doma.lsc' : [],
                   }

        received = find_moving_object_network_up_intervals(window, self.elements,
                                                           site_filename, chunksize)

        for site in expected.keys():
            print site
            assert_equal(expected[site], received[site])



    def test_chunk_windows(self):
        window = {
                   'start'    : str_to_dt('2013-12-05 01:38:53'),
                   'end'      : str_to_dt('2013-12-05 08:33:55'),
                 }

        chunksize = timedelta(hours=2)

        expected = [
                     (
                       str_to_dt('2013-12-05 01:38:53'),
                       str_to_dt('2013-12-05 03:38:53')
                     ),
                     (
                       str_to_dt('2013-12-05 03:38:53'),
                       str_to_dt('2013-12-05 05:38:53')
                     ),
                     (
                       str_to_dt('2013-12-05 05:38:53'),
                       str_to_dt('2013-12-05 07:38:53')
                     ),
                     (
                       str_to_dt('2013-12-05 07:38:53'),
                       str_to_dt('2013-12-05 08:33:55')
                     )
                   ]


        received = chunk_windows(window, chunksize)

        assert_equal(received, expected)


    def test_chunk_windows_span_ut_boundary(self):
        window = {
                   'start'    : str_to_dt('2013-12-04 23:38:53'),
                   'end'      : str_to_dt('2013-12-05 06:33:55'),
                 }

        chunksize = timedelta(hours=2)

        expected = [
                     (
                       str_to_dt('2013-12-04 23:38:53'),
                       str_to_dt('2013-12-05 01:38:53')
                     ),
                     (
                       str_to_dt('2013-12-05 01:38:53'),
                       str_to_dt('2013-12-05 03:38:53')
                     ),
                     (
                       str_to_dt('2013-12-05 03:38:53'),
                       str_to_dt('2013-12-05 05:38:53')
                     ),
                     (
                       str_to_dt('2013-12-05 05:38:53'),
                       str_to_dt('2013-12-05 06:33:55')
                     )
                   ]


        received = chunk_windows(window, chunksize)

        assert_equal(received, expected)
