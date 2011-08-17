#!/usr/bin/python
from __future__ import division

# Standard library imports
import datetime
from nose.tools import assert_equal

# Module imports
from rise_set.angle import Angle
from rise_set.sky_coordinates import RightAscension, Declination
from rise_set.rates import ProperMotion

#Import the module to test
from rise_set.telescope import Telescope


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
                         'ra_proper_motion'  : ProperMotion(RightAscension('00 00 00.01999')),
                         'dec_proper_motion' : ProperMotion(Declination('00 00 00.02367')),
                         'parallax'          : 0.01043,   # Units: arcsec
                         'rad_vel'           : 20.5,      # Units: km/s (-ve approaches)
                         'epoch'             : 2000,
                       }

        # Site (East +ve longitude)
        self.siding_spring = Telescope(latitude  = Angle(degrees = -31.273),
                                       longitude = Angle(degrees = 149.070593))

        # Date
        self.date = datetime.datetime(year = 2010, month = 3, day = 12)

        self.siding_spring.set_date(self.date)


    def test_rise_set(self):
        (transit, rise, set) = self.siding_spring.calc_rise_set(self.canopus)

        # Second accuracy isn't meaningful, so just check we would round up or
        # down as apppropriate.
        exp_secs = 30

        # Expected times taken from http://aa.usno.navy.mil/data/docs/mrst.php

        # Transit
        exp_t_hr  = 9
        exp_t_min = 8

        transit_time = transit + self.date

        assert_equal(transit_time.hour, exp_t_hr)
        assert_equal(transit_time.minute, exp_t_min)
        assert transit_time.second <= exp_secs, '%r !> %r' % (transit_time.second, exp_secs)


        # Rise
        exp_r_hr  = 23
        exp_r_min = 27

        rise_time = rise + self.date

        assert_equal(rise_time.hour, exp_r_hr)
        assert_equal(rise_time.minute, exp_r_min)
        assert rise_time.second <= exp_secs, '%r !> %r' % (rise_time.second, exp_secs)


        # Set
        exp_s_hr  = 18
        exp_s_min = 45    # USNO actually gives 18.46, but let's not worry...

        set_time = set + self.date

        assert_equal(set_time.hour, exp_s_hr)
        assert_equal(set_time.minute, exp_s_min)
        assert set_time.second <= exp_secs, '%r !< %r' % (set_time.second, exp_secs)
