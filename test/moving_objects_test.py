#!/usr/bin/env python

'''
moving_objects_test.py - Tests for moving object support

Authors: Eric Saunders
         Tim Lister
December 2013
'''
from __future__ import print_function
from __future__ import division
from builtins import zip
from builtins import str
from past.utils import old_div
from builtins import object


from rise_set.angle import Angle
from rise_set.moving_objects import (
                                     Sites, initialise_sites,
                                     chunk_windows,
                                     read_neocp_orbit,
                                     calc_ephemerides,
                                     hour_angle_within_limits,
                                     ephemeris_chunk_within_ha_limits,
                                     ephemeris_chunk_above_horizon,
                                     find_moving_object_up_intervals,
                                     find_moving_object_network_up_intervals,
                                     )
from rise_set.utils import MovingViolation, is_moving_object
from rise_set.astrometry import elem_to_topocentric_apparent

from datetime import datetime, timedelta
from nose.tools import (assert_equal, assert_almost_equal, assert_raises,
                        assert_false, assert_true)
from mock import patch


def str_to_dt(dt_str):
   '''Simple helper to make the tests more readable.'''
   return datetime.strptime(dt_str, '%Y-%m-%d %H:%M:%S')


class TestSites(object):

    def setup(self):
        self.site_dict1 = {
                             'name'         : '1m0a.doma.cpt',
                             'latitude'     : -32.38059,
                             'longitude'    :  20.8101083333,
                             'horizon'      : 30,
                             'ha_limit_neg' : -12.0,
                             'ha_limit_pos' : 12.0,
                          }
        self.site_dict2 = {
                             'name'         : '1m0a.doma.lsc',
                             'latitude'     : -30.1673472222,
                             'longitude'    : -70.8046722222,
                             'horizon'      : 30,
                             'ha_limit_neg' : -12.0,
                             'ha_limit_pos' : 12.0,
                          }
        self.site_dict3 = {
                             'name'         : '1m0b.domb.cpt',
                             'latitude'     : -32.38059,
                             'longitude'    :  20.8101083333,
                             'horizon'      : 30,
                             'ha_limit_neg' : -12.0,
                             'ha_limit_pos' : 12.0,
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
        expected_keys = set(['1m0a.doma.elp', '1m0a.doma.coj',
                         '1m0a.doma.cpt', '1m0a.doma.lsc'])
        sites = initialise_sites('test/telescopes.dat')

        assert_equal(len(sites.sites), 4)
        assert_equal(set(sites.sites.keys()), expected_keys)




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
                          'type'           : 'MPC_MINOR_PLANET',
                          'epoch'          : 56600.0,
                          'inclination'    : Angle(degrees=7.22565),
                          'long_node'      : Angle(degrees=173.25052),
                          'arg_perihelion' : Angle(degrees=47.60658),
                          'semi_axis'      : 2.4050925,
                          'eccentricity'   : 0.0943494,
                          'mean_anomaly'   : Angle(degrees=54.47380),
                        }

        # Element set for Comet C/2012 K1, formatted as an asteroid which
        # caused serious grief to the scheduler (Issue #7447)

        self.elements2 = {
                          'epoch'          : 56800.0,
                          'inclination'    : Angle(degrees=142.4282),
                          'long_node'      : Angle(degrees=317.7384),
                          'arg_perihelion' : Angle(degrees=203.1103),
                          'semi_axis'      : 2.9612496,
                          'eccentricity'   : 1.0,
                          'mean_anomaly'   : Angle(degrees=27.6567),
                        }

        # Fake element set to test MovingViolation handling in
        # elem_to_topocentric_apparent

        self.elements3 = {
                          'epoch'          : 56800.0,
                          'inclination'    : Angle(degrees=142.4282),
                          'long_node'      : Angle(degrees=317.7384),
                          'arg_perihelion' : Angle(degrees=203.1103),
                          'semi_axis'      : -2.9612496,
                          'eccentricity'   : 0.99,
                          'mean_anomaly'   : Angle(degrees=27.6567),
                        }
        # Element set for Comet C/2012 K1 from JPL HORIZONS on 2014-06-11

        self.comet_elements1 = {
                          'epochofperih'   : 56896.65636138478,
                          'inclination'    : Angle(degrees=142.4281063052002),
                          'long_node'      : Angle(degrees=317.731588784467),
                          'arg_perihelion' : Angle(degrees=203.0929610210846),
                          'perihdist'      : 1.054777188675267,
                          'eccentricity'   : 1.000209136966309,
                        }


        # Element set for 2013 UQ4 from JPL HORIZONS on 2014-06-11

        self.comet_elements2 = {
                          'epochofperih'   : 56843.96309797978,
                          'inclination'    : Angle(degrees=145.25830684677),
                          'long_node'      : Angle(degrees=317.6619568986719),
                          'arg_perihelion' : Angle(degrees=23.31494833774585),
                          'perihdist'      : 1.080940992734334,
                          'eccentricity'   : 0.9820934471344616,
                        }

        self.elp = {
                     'latitude'  : Angle(degrees=30.6801),
                     'longitude' : Angle(degrees=-104.015194444),
                   }

        self.HOURS_TO_DEGREES = 15

        self.cpt = {
                        'name'      : '1m0a.doma.cpt',
                        'latitude'  : Angle(degrees=-32.38059),
                        'longitude' : Angle(degrees=20.8101083333),
                        'horizon'      : Angle(degrees=30),
                        'ha_limit_neg' : Angle(degrees=-12.0*self.HOURS_TO_DEGREES),
                        'ha_limit_pos' : Angle(degrees=12.0*self.HOURS_TO_DEGREES),
                    }


    def test_is_moving_object_defaults_to_extra_solar_objects(self):
        target = {}
        assert_false(is_moving_object(target))

    def test_is_moving_object_understands_planets(self):
        target = {'type' : 'MPC_MINOR_PLANET'}
        assert_true(is_moving_object(target))

    def test_is_moving_object_understands_comets(self):
        target = {'type' : 'MPC_MINOR_PLANET'}
        assert_true(is_moving_object(target))

    def test_is_moving_object_is_case_insensitive(self):
        target = {'type' : 'MpC_MiNoR_PlAnEt'}
        assert_true(is_moving_object(target))

    def test_is_moving_object_other_types_dont_move(self):
        target = {'type' : 'TIMS_CAT'}
        assert_false(is_moving_object(target))


    def test_read_neocp_orbit1(self):

        # Element dictionary produced by TAL's code from NEOCP orbit file, migrated to
        # rise_set format.
        expected = {
                     'name': 'P109rXK',
                     'epoch': 56620.0,
                     'long_node':  Angle(degrees=224.42273),
                     'eccentricity': 0.5131088,
                     'semi_axis': 2.0773493,
                     'mean_anomaly':  Angle(degrees=322.22271),
                     'arg_perihelion': Angle(degrees=307.96804),
                     'inclination':  Angle(degrees=13.72592),
                     'MDM':  Angle(degrees=0.3291848),
                     'H': 17.5,
                     'G': 0.15,
                     'n_obs': 4,
                     'n_oppos' : 1,
                     'n_nights': 0,
                   }

        recieved = read_neocp_orbit('test/P109rXK.neocp')

        for e in list(expected.keys()):
            assert_equal(expected[e], recieved[e])

    def test_read_neocp_orbit2(self):

        # Element dictionary produced by TAL's code from NEOCP orbit file, migrated to
        # rise_set format.
        expected = {
                     'name'           : 'K13TB7L',
                     'epoch'          : 56600.0,
                     'long_node'      : Angle(degrees=3.36296),
                     'eccentricity'   : 0.6898696,
                     'semi_axis'      : 3.6038493,
                     'mean_anomaly'   : Angle(degrees=344.69702),
                     'arg_perihelion' : Angle(degrees=112.20127),
                     'inclination'    : Angle(degrees=9.36592),
                     'MDM'            : Angle(degrees=0.14406356),
                     'H'              : 17.7,
                     'G'              : 0.15,
                     'n_obs'          : 146,
                     'n_oppos'        : 1,
                     'n_nights'       : 67,
                     'reference'      : 'E2013-X58',
                     'uncertainty'    : '3',
                   }

        recieved = read_neocp_orbit('test/2013TL117.neocp')

        for e in list(expected.keys()):
            assert_equal(expected[e], recieved[e])


    def test_read_neocp_orbit3(self):

        # Element dictionary produced by TAL's code from NEOCP orbit file, migrated to
        # rise_set format.
        # This test tests multi-opposition orbits
        expected = {
                     'name'           : '00001',
                     'epoch'          : 57200.0,
                     'long_node'      : Angle(degrees=80.3272),
                     'eccentricity'   : 0.0757825,
                     'semi_axis'      : 2.7679724,
                     'mean_anomaly'   : Angle(degrees=138.66222),
                     'arg_perihelion' : Angle(degrees= 72.65410),
                     'inclination'    : Angle(degrees= 10.59230),
                     'MDM'            : Angle(degrees=0.21402349),
                     'H'              :  3.34,
                     'G'              : 0.12,
                     'n_obs'          : 6541,
                     'n_oppos'        : 107,
                     'n_nights'       : '1802-2014',
                     'reference'      : 'MPO328644',
                     'uncertainty'    : '0',
                   }

        recieved = read_neocp_orbit('test/00001_Ceres.neocp')

        for e in list(expected.keys()):
            assert_equal(expected[e], recieved[e])



    def test_read_neocp_orbit4(self):

        # Element dictionary produced by TAL's code from NEOCP orbit file, migrated to
        # rise_set format.
        # This test tests very short orbits from find_orb with an arc length in
        # minutes
        expected = {
                     'name'           : 'LSCTLF3',
                     'epoch'          : 57109.0,
                     'semi_axis'      : 2.2601939,
                     'mean_anomaly'   : Angle(degrees= 23.58746),
                     'arg_perihelion' : Angle(degrees= 75.79622),
                     'long_node'      : Angle(degrees= 69.92181),
                     'inclination'    : Angle(degrees=  9.43941),
                     'eccentricity'   : 0.3369232,
                     'MDM'            : Angle(degrees=0.29005846),
                     'H'              : 18.90,
                     'G'              : 0.15,
                     'n_obs'          : 6,
                     'n_oppos'        : 1,
                     'n_nights'       : old_div(18.8,1440.0),    # 18.8 mins->days
                     'reference'      : 'FO 150328',
                     'residual'       : 0.08,
                     'uncertainty'    : 'U',
                   }

        recieved = read_neocp_orbit('test/LSCTLF3.neocp')

        for e in list(expected.keys()):
            assert_equal(expected[e], recieved[e])



    def test_read_neocp_orbit5(self):

        # Element dictionary produced by TAL's code from NEOCP orbit file, migrated to
        # rise_set format.
        # This test tests very short orbits from find_orb with an arc length in
        # hours
        expected = {
                     'name'           : 'LSCTLF3',
                     'epoch'          : 57109.0,
                     'mean_anomaly'   : Angle(degrees= 21.73763),
                     'arg_perihelion' : Angle(degrees= 80.02788),
                     'long_node'      : Angle(degrees= 64.35739),
                     'inclination'    : Angle(degrees= 12.17788),
                     'eccentricity'   : 0.3895454,
                     'MDM'            : Angle(degrees=0.21752314),
                     'semi_axis'      : 2.7382037,
                     'H'              : 18.23,
                     'G'              : 0.15,
                     'n_obs'          : 10,
                     'n_oppos'        : 1,
                     'n_nights'       : old_div(13.0,24.0),    # 13.0 hrs->days
                     'reference'      : 'FO 150328',
                     'residual'       : 0.09,
                     'uncertainty'    : 'U',
                   }

        recieved = read_neocp_orbit('test/LSCTLF3_V2.neocp')

        for e in list(expected.keys()):
            assert_equal(expected[e], recieved[e])



    def test_read_neocp_orbit6(self):

        # Element dictionary produced by TAL's code from NEOCP orbit file, migrated to
        # rise_set format.
        # This test tests NEOCP orbits without a rms residual.
        expected = {
                     'name'           : 'N007rdx',
                     'epoch'          : 57120.0,
                     'mean_anomaly'   : Angle(degrees=119.11683),
                     'arg_perihelion' : Angle(degrees= 75.91136),
                     'long_node'      : Angle(degrees=351.90874),
                     'inclination'    : Angle(degrees=  3.64505),
                     'eccentricity'   : 0.1054741,
                     'MDM'            : Angle(degrees=1.07358714),
                     'semi_axis'      : 0.9445925,
                     'H'              : 21.8,
                     'G'              : 0.15,
                     'n_obs'          : 5,
                     'n_oppos'        : 1,
                     'n_nights'       : 0,
                     'residual'       : 9.99,   # No value, assume default
                     'reference'      : '',
                   }

        recieved = read_neocp_orbit('test/N007rdx.neocp')

        for e in list(expected.keys()):
            assert_equal(expected[e], recieved[e])


    def test_read_neocp_orbit7(self):

        # Element dictionary produced by TAL's code from NEOCP orbit file, migrated to
        # rise_set format.
        # This test tests NEOCP orbits without an uncertainty.
        expected = {
                     'name'           : 'N007rdx',
                     'epoch'          : 57120.0,
                     'mean_anomaly'   : Angle(degrees=119.11683),
                     'arg_perihelion' : Angle(degrees= 75.91136),
                     'long_node'      : Angle(degrees=351.90874),
                     'inclination'    : Angle(degrees=  3.64505),
                     'eccentricity'   : 0.1054741,
                     'MDM'            : Angle(degrees=1.07358714),
                     'semi_axis'      : 0.9445925,
                     'H'              : 21.8,
                     'G'              : 0.15,
                     'n_obs'          : 5,
                     'n_oppos'        : 1,
                     'n_nights'       : 0,
                     'uncertainty'    : 'U',   # No value, assume default
                     'reference'      : '',
                   }

        recieved = read_neocp_orbit('test/N007rdx.neocp')

        for e in list(expected.keys()):
            assert_equal(expected[e], recieved[e])


    def test_read_neocp_comet_orbit1(self):

        # Hand-crafted element dictionary from MPC MPES output on 2014-09-26.
        # This test tests long-period comet parsing
        expected = {
                     'name'           : 'CK13A010',
                     'epoch'          : 57000.0,
                     'long_node'      : Angle(degrees=300.9764),
                     'eccentricity'   : 1.000434,
                     'semi_axis'      : None,
                     'mean_anomaly'   : None,
                     'arg_perihelion' : Angle(degrees=2.4226),
                     'inclination'    :  Angle(degrees=129.0428),
                     'MDM'            :  None,
                     'H'              : 8.2,
                     'G'              : 2.4,
                     'n_obs'          : None,
                     'n_nights'       : None,
                     'type'           : 'MPC_COMET'
                   }

        recieved = read_neocp_orbit('test/Comet_SidingSpring.neocp')

        for e in list(expected.keys()):
            assert_equal(expected[e], recieved[e])

    def test_read_neocp_comet_orbit2(self):

        # Hand-crafted element dictionary from MPC MPES output on 2014-09-26.
        # This test tests periodic comet parsing
        expected = {
                     'name'           : '0299P',
                     'epoch'          : 57000.0,
                     'long_node'      : Angle(degrees=271.6800),
                     'eccentricity'   : 0.282122,
                     'semi_axis'      : None,
                     'mean_anomaly'   : None,
                     'arg_perihelion' : Angle(degrees=323.5091),
                     'inclination'    : Angle(degrees=10.4798),
                     'MDM'            : None,
                     'H'              : 11.5,
                     'G'              : 4.00,
                     'n_obs'          : None,
                     'n_nights'       : None,
                     'type'           : 'MPC_COMET'
                   }

        recieved = read_neocp_orbit('test/Comet_299P.neocp')

        for e in list(expected.keys()):
            assert_equal(expected[e], recieved[e])

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


    def test_elem_to_topocentric_apparent_bad_eccen(self):
        dt = datetime(2014,  6, 14)
        tdb = dt - timedelta(seconds=67.184)

        assert_raises(MovingViolation, elem_to_topocentric_apparent, tdb, self.elements2, self.elp)
        try:
            ra, dec = elem_to_topocentric_apparent(tdb, self.elements2, self.elp)
        except MovingViolation as e:
            assert_equal(str(e), 'Error: -2 (illegal eccentricity)')
        else:
            assert False, "Didn't raise expected MovingViolation error"


    def test_elem_to_topocentric_apparent_bad_meandist(self):
        dt = datetime(2014,  6, 14)
        tdb = dt - timedelta(seconds=67.184)

        assert_raises(MovingViolation, elem_to_topocentric_apparent, tdb, self.elements3, self.elp)
        try:
            ra, dec = elem_to_topocentric_apparent(tdb, self.elements3, self.elp)
        except MovingViolation as e:
            assert_equal(str(e), 'Error: -3 (illegal mean distance)')
        else:
            assert False, "Didn't raise expected MovingViolation error"


    def test_elem_to_topocentric_apparent_bad_jform(self):
        dt = datetime(2014,  6, 14)
        tdb = dt - timedelta(seconds=67.184)

        assert_raises(MovingViolation, elem_to_topocentric_apparent, tdb, self.comet_elements1, self.elp, 42)
        try:
            ra, dec = elem_to_topocentric_apparent(tdb, self.comet_elements1, self.elp, 42)
        except MovingViolation as e:
            assert_equal(str(e), 'Error: -1 (illegal JFORM)')
        else:
            assert False, "Didn't raise expected MovingViolation error"


    def test_elem_to_topocentric_apparent_comet1(self):
        dt = datetime(2014, 6, 14)
        tdb = dt - timedelta(seconds=67.184)

        ra, dec = elem_to_topocentric_apparent(tdb, self.comet_elements1, self.elp, 3)


        # Coordinates from JPL Horizons: 153.36306  36.46146
        assert_almost_equal(ra.in_degrees(), 153.36306, places=2)
        assert_almost_equal(dec.in_degrees(), 36.46146, places=2)


    def test_elem_to_topocentric_apparent_comet2(self):
        dt = datetime(2013, 12, 12, 6, 10)
        tdb = dt - timedelta(seconds=67.184)

        ra, dec = elem_to_topocentric_apparent(tdb, self.comet_elements2, self.cpt, 3)


        # Coordinates from JPL Horizons:  37.78519 -31.04102
        assert_almost_equal(ra.in_degrees(),  37.78519, places=4)
        assert_almost_equal(dec.in_degrees(),-31.04102, places=4)


    def test_calc_ephemerides(self):
        window = {
                   'start' : datetime(2013, 12, 10),
                   'end'   : datetime(2013, 12, 11)
                 }

        chunk_size = timedelta(hours=6)

        expected = [
                     {
                       'ra_app'  : Angle(degrees=285.22681),
                       'dec_app' : Angle(degrees=-18.16350),
                       'start'   : window['start'],
                       'end'     : window['start'] + (1 * chunk_size),
                     },
                     {
                       'ra_app'  : Angle(degrees=285.35677),
                       'dec_app' : Angle(degrees=-18.15597),
                       'start'   : window['start'] + (1 * chunk_size),
                       'end'     : window['start'] + (2 * chunk_size),
                     },
                     {
                       'ra_app'  : Angle(degrees=285.48538),
                       'dec_app' : Angle(degrees=-18.14841),
                       'start'   : window['start'] + (2 * chunk_size),
                       'end'     : window['start'] + (3 * chunk_size),
                     },
                     {
                       'ra_app'  : Angle(degrees=285.61369),
                       'dec_app' : Angle(degrees=-18.14033),
                       'start'   : window['start'] + (3 * chunk_size),
                       'end'     : window['start'] + (4 * chunk_size),
                     },
                     {
                       'ra_app'  : Angle(degrees=285.74337),
                       'dec_app' : Angle(degrees=-18.13206),
                       'start'   : window['start'] + (4 * chunk_size),
                       'end'     : window['start'] + (4 * chunk_size),
                     },
                   ]

        received = calc_ephemerides(window, self.elements, self.cpt, chunk_size)

        for e, r in zip(expected, received):
            assert_almost_equal(e['ra_app'].in_degrees(), r['ra_app'].in_degrees(), places=3)
            assert_almost_equal(e['dec_app'].in_degrees(), r['dec_app'].in_degrees(), places=4)
            assert_equal(e['start'], r['start'])
            assert_equal(e['end'], r['end'])


    def test_find_moving_object_up_intervals1(self):
        window = {
                   'start' : datetime(2013, 12, 10),
                   'end'   : datetime(2013, 12, 11)
                 }

        chunk_size = timedelta(hours=3)

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
                                                               self.cpt, chunk_size)

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

        chunk_size = timedelta(minutes=15)

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
                                                               self.cpt, chunk_size)

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

        chunk_size = timedelta(minutes=15)

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
                                                               self.cpt, chunk_size)

        assert_equal(len(expected_intervals), len(received_ints))
        assert_equal(len(expected_altitudes), len(received_alts))

        for e, r in zip(expected_intervals, received_ints):
            assert_equal(e, r)

        for e, r in zip(expected_altitudes, received_alts):
            assert_almost_equal(e.in_degrees(), r.in_degrees(), places=3)


    @patch('rise_set.moving_objects.calc_ephemerides')
    @patch('rise_set.moving_objects.calc_local_hour_angle')
    @patch('rise_set.moving_objects.calculate_altitude')
    def test_find_moving_object_sets_asteroid_jform(self, alt_mock, ha_mock,
                                                    ephem_mock):
        window   = None
        elements = {'type' : 'MPC_MINOR_PLANET'}
        site     = None
        find_moving_object_up_intervals(window, elements, site)
        assert_equal(ephem_mock.call_args[0][4], 2)


    @patch('rise_set.moving_objects.calc_ephemerides')
    @patch('rise_set.moving_objects.calc_local_hour_angle')
    @patch('rise_set.moving_objects.calculate_altitude')
    def test_find_moving_object_sets_comet_jform(self, alt_mock, ha_mock,
                                                    ephem_mock):
        window   = None
        elements = {'type' : 'MPC_COMET'}
        site     = None
        find_moving_object_up_intervals(window, elements, site)
        assert 3 in ephem_mock.call_args[0], 'expected jform=3 not passed to calc_ephemerides'


    def test_find_moving_object_invalid_type_raises_exception(self):
        window   = None
        elements = {'type' : 'TIMS_CAT'}
        site     = None

        try:
            find_moving_object_up_intervals(window, elements, site)
        except MovingViolation as e:
            assert_equal(str(e), "Unsupported target type: 'tims_cat'")
        else:
            assert False, "Didn't raise expected MovingViolation error"


    def test_hour_angle_beyond_neg_limit(self):
        ha        = Angle(degrees=-8*self.HOURS_TO_DEGREES)
        neg_limit = Angle(degrees=-4.6*self.HOURS_TO_DEGREES)
        pos_limit = Angle(degrees=4.6*self.HOURS_TO_DEGREES)

        assert_equal(hour_angle_within_limits(ha, neg_limit, pos_limit), False)


    def test_hour_angle_beyond_pos_limit(self):
        ha        = Angle(degrees=8*self.HOURS_TO_DEGREES)
        neg_limit = Angle(degrees=-4.6*self.HOURS_TO_DEGREES)
        pos_limit = Angle(degrees=4.6*self.HOURS_TO_DEGREES)

        assert_equal(hour_angle_within_limits(ha, neg_limit, pos_limit), False)


    def test_hour_angle_within_limits(self):
        ha        = Angle(degrees=3*self.HOURS_TO_DEGREES)
        neg_limit = Angle(degrees=-4.6*self.HOURS_TO_DEGREES)
        pos_limit = Angle(degrees=4.6*self.HOURS_TO_DEGREES)

        assert_equal(hour_angle_within_limits(ha, neg_limit, pos_limit), True)


    def test_ephemeris_chunk_within_ha_limits(self):
        ha1        = Angle(degrees=3*self.HOURS_TO_DEGREES)
        ha2        = Angle(degrees=3.5*self.HOURS_TO_DEGREES)
        neg_limit  = Angle(degrees=-4.6*self.HOURS_TO_DEGREES)
        pos_limit  = Angle(degrees=4.6*self.HOURS_TO_DEGREES)

        assert_equal(ephemeris_chunk_within_ha_limits(ha1, ha2, neg_limit, pos_limit), True)


    def test_ephemeris_chunk_partially_outside_positive_ha_limit(self):
        ha1        = Angle(degrees=4.5*self.HOURS_TO_DEGREES)
        ha2        = Angle(degrees=4.7*self.HOURS_TO_DEGREES)
        neg_limit  = Angle(degrees=-4.6*self.HOURS_TO_DEGREES)
        pos_limit  = Angle(degrees=4.6*self.HOURS_TO_DEGREES)

        assert_equal(ephemeris_chunk_within_ha_limits(ha1, ha2, neg_limit, pos_limit), False)


    def test_ephemeris_chunk_fully_outside_positive_ha_limit(self):
        ha1        = Angle(degrees=4.7*self.HOURS_TO_DEGREES)
        ha2        = Angle(degrees=4.8*self.HOURS_TO_DEGREES)
        neg_limit  = Angle(degrees=-4.6*self.HOURS_TO_DEGREES)
        pos_limit  = Angle(degrees=4.6*self.HOURS_TO_DEGREES)

        assert_equal(ephemeris_chunk_within_ha_limits(ha1, ha2, neg_limit, pos_limit), False)


    def test_ephemeris_chunk_partially_outside_negative_ha_limit(self):
        ha1        = Angle(degrees=-4.5*self.HOURS_TO_DEGREES)
        ha2        = Angle(degrees=-4.7*self.HOURS_TO_DEGREES)
        neg_limit  = Angle(degrees=-4.6*self.HOURS_TO_DEGREES)
        pos_limit  = Angle(degrees=4.6*self.HOURS_TO_DEGREES)

        assert_equal(ephemeris_chunk_within_ha_limits(ha1, ha2, neg_limit, pos_limit), False)


    def test_ephemeris_chunk_fully_outside_negative_ha_limit(self):
        ha1        = Angle(degrees=-4.7*self.HOURS_TO_DEGREES)
        ha2        = Angle(degrees=-4.8*self.HOURS_TO_DEGREES)
        neg_limit  = Angle(degrees=-4.6*self.HOURS_TO_DEGREES)
        pos_limit  = Angle(degrees=4.6*self.HOURS_TO_DEGREES)

        assert_equal(ephemeris_chunk_within_ha_limits(ha1, ha2, neg_limit, pos_limit), False)


    def test_ephemeris_chunk_fully_above_horizon(self):
        alt1    = Angle(degrees=20)
        alt2    = Angle(degrees=21)
        horizon = Angle(degrees=15)

        assert_equal(ephemeris_chunk_above_horizon(alt1, alt2, horizon), True)


    def test_ephemeris_chunk_partially_above_horizon(self):
        alt1    = Angle(degrees=14)
        alt2    = Angle(degrees=16)
        horizon = Angle(degrees=15)

        assert_equal(ephemeris_chunk_above_horizon(alt1, alt2, horizon), False)


    def test_ephemeris_chunk_below_horizon(self):
        alt1    = Angle(degrees=10)
        alt2    = Angle(degrees=11)
        horizon = Angle(degrees=15)

        assert_equal(ephemeris_chunk_above_horizon(alt1, alt2, horizon), False)


    @patch('rise_set.moving_objects.calc_ephemerides')
    @patch('rise_set.moving_objects.calc_local_hour_angle')
    @patch('rise_set.moving_objects.calculate_altitude')
    def test_moving_objects_respect_negative_hour_angle_limit(self, alt_mock,
                                                              ha_mock, ephem_mock):
        window   = None
        elements = {'type' : 'mpc_minor_planet'}
        site_dict = {
                        'name'         : '1m0a.doma.cpt',
                        'latitude'     : Angle(degrees=-32.38059),
                        'longitude'    : Angle(degrees=20.8101083333),
                        'horizon'      : Angle(degrees=15),
                        'ha_limit_neg' : Angle(degrees=-5.0*15),
                        'ha_limit_pos' : Angle(degrees=5.0*15),
                    }

        ephem = [
                    {
                      'start' : 1,
                      'end'   : 2,
                      'ra_app' : Angle(3),
                      'dec_app' : Angle(4),
                    },
                    {
                      'start' : 5,
                      'end'   : 6,
                      'ra_app' : Angle(7),
                      'dec_app' : Angle(8),
                    },
               ]

        chunk_size = timedelta(minutes=15)

        alt_mock.return_value   = Angle(degrees=20)
        ha_mock.return_value    = Angle(degrees=-7*15)
        ephem_mock.return_value = ephem

        received_ints, received_alts = find_moving_object_up_intervals(window, elements,
                                                                  site_dict, chunk_size)
        assert_equal(received_ints, [])


    def test_find_moving_object_network_up_intervals(self):
        window = {
                   'start' : datetime(2013, 12, 10, 7, 30),
                   'end'   : datetime(2013, 12, 10, 8, 30)
                 }

        site_filename = 'test/thirty_degree_telescopes.dat'

        chunk_size = timedelta(minutes=15)

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
                                                           site_filename, chunk_size)

        for site in list(expected.keys()):
            print(site)
            assert_equal(expected[site], received[site])


    def test_chunk_windows(self):
        window = {
                   'start'    : str_to_dt('2013-12-05 01:38:53'),
                   'end'      : str_to_dt('2013-12-05 08:33:55'),
                 }

        chunk_size = timedelta(hours=2)

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


        received = chunk_windows(window, chunk_size)

        assert_equal(received, expected)


    def test_chunk_windows_span_ut_boundary(self):
        window = {
                   'start'    : str_to_dt('2013-12-04 23:38:53'),
                   'end'      : str_to_dt('2013-12-05 06:33:55'),
                 }

        chunk_size = timedelta(hours=2)

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


        received = chunk_windows(window, chunk_size)

        assert_equal(received, expected)


    def test_special_neo_target(self):
        window = {
                   'start'    : str_to_dt('2014-07-17 23:00:00'),
                   'end'      : str_to_dt('2014-07-18 11:00:00'),
                 }

        elements = read_neocp_orbit('test/special_neo_target.neocp')

        for x, y in list(elements.items()):
            print(x, y)

        site_filename = 'test/telescopes.dat'
        chunk_size = timedelta(minutes=10)

        received = find_moving_object_network_up_intervals(window, elements,
                                                           site_filename, chunk_size)

        expected = {
                     '1m0a.doma.elp' : [
                                         (
                                           datetime(2014, 7, 18, 3, 50),
                                           datetime(2014, 7, 18, 9, 10)
                                         )
                                       ],
                     '1m0a.doma.coj' : [
                                         (
                                           datetime(2014, 7, 18, 8, 50),
                                           datetime(2014, 7, 18, 11, 0)
                                         )
                                       ],
                     '1m0a.doma.cpt' : [
                                         (
                                           datetime(2014, 7, 17, 23, 0),
                                           datetime(2014, 7, 18, 3, 10)
                                         )
                                       ],
                     '1m0a.doma.lsc' : [
                                         (
                                           datetime(2014, 7, 17, 23, 30),
                                           datetime(2014, 7, 18, 9, 20)
                                         )
                                       ],
                   }


        for site in list(expected.keys()):
            print(site)
            assert_equal(expected[site], received[site])
