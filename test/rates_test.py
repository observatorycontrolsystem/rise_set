#!/usr/bin/env python

from __future__ import division
from builtins import object

from nose.tools import assert_equal, assert_almost_equal, raises
from math import pi

#Import the module to test
from rise_set.rates import ProperMotion, RatesConfigError
from rise_set.sky_coordinates import RightAscension, Declination


class TestProperMotion(object):
    '''Unit tests for the ProperMotion object.'''

    #Test that conifg errors are raised if wrong key words are inputted
    @raises(RatesConfigError)
    def test_configuration(self):
        ra = RightAscension(degrees = 90)
        proper_motion = ProperMotion(ra, time = 'millennia')



    # Generate proper motion for right ascension of degrees in arc
    # Return in same time units as entered
    def test_in_degrees_per_year_ra_deg_milliarcseconds(self):
        ra = RightAscension(degrees = '0 0 0.001')
        proper_motion = ProperMotion(ra)
        assert_almost_equal(proper_motion.in_degrees_per_year(),  2.77778e-07, 5)

    def test_in_radians_per_year_ra_deg_arc(self):
        ra  = RightAscension(degrees = 90)
        proper_motion = ProperMotion(ra)
        assert_equal(proper_motion.in_radians_per_year(), pi/2)

    def test_in_degrees_per_year_ra_deg_arc(self):
        ra  = RightAscension(degrees = 90)
        proper_motion = ProperMotion(ra)
        assert_equal(proper_motion.in_degrees_per_year(), 90)

    def test_in_radians_per_century_ra_deg_arc(self):
        ra  = RightAscension(degrees = 90)
        proper_motion = ProperMotion(ra, time = 'century')
        assert_equal(proper_motion.in_radians_per_century(), pi/2)

    def test_in_degrees_per_century_ra_deg_arc(self):
        ra  = RightAscension(degrees = 90)
        proper_motion = ProperMotion(ra, time = 'century')
        assert_equal(proper_motion.in_degrees_per_century(), 90)



    # Generate proper motion for right ascension of degrees in arc
    # Return in different time units as entered
    def test_in_radians_per_year_from_century_ra_deg_arc(self):
        ra  = RightAscension(degrees = 90)
        proper_motion = ProperMotion(ra, time = 'century')
        assert_equal(proper_motion.in_radians_per_year(), pi/2 / 100 )

    def test_in_degrees_per_year_from_century_ra_deg_arc(self):
        ra  = RightAscension(degrees = 90)
        proper_motion = ProperMotion(ra, time = 'century')
        assert_equal(proper_motion.in_degrees_per_year(), 90 / 100)

    def test_in_radians_per_century_from_year_ra_deg_arc(self):
        ra  = RightAscension(degrees = 90)
        proper_motion = ProperMotion(ra)
        assert_equal(proper_motion.in_radians_per_century(), pi/2 * 100)

    def test_in_degrees_per_century_from_year_ra_deg_arc(self):
        ra  = RightAscension(degrees = 90)
        proper_motion = ProperMotion(ra)
        assert_equal(proper_motion.in_degrees_per_century(), 90 * 100)



    # Generate proper motion for right ascension of degrees in time
    # Return in same time units as entered
    def test_in_radians_per_year_ra_time(self):
        ra  = RightAscension('12:00:00')
        proper_motion = ProperMotion(ra)
        assert_equal(proper_motion.in_radians_per_year(), pi)

    def test_in_degrees_per_year_ra_time(self):
        ra  = RightAscension('12:00:00')
        proper_motion = ProperMotion(ra)
        assert_equal(proper_motion.in_degrees_per_year(), 180)

    def test_in_radians_per_century_ra_time(self):
        ra  = RightAscension('12:00:00')
        proper_motion = ProperMotion(ra, time = 'century')
        assert_equal(proper_motion.in_radians_per_century(), pi)

    def test_in_degrees_per_century_ra_time(self):
        ra  = RightAscension('12:00:00')
        proper_motion = ProperMotion(ra, time = 'century')
        assert_equal(proper_motion.in_degrees_per_century(), 180)



    # Generate proper motion for right ascension of degrees in time
    # Return in different time units as entered
    def test_in_radians_per_year_from_century_ra_time(self):
        ra  = RightAscension('12:00:00')
        proper_motion = ProperMotion(ra, time = 'century')
        assert_equal(proper_motion.in_radians_per_year(), pi / 100 )

    def test_in_degrees_per_year_from_century_ra_time(self):
        ra  = RightAscension('12:00:00')
        proper_motion = ProperMotion(ra, time = 'century')
        assert_equal(proper_motion.in_degrees_per_year(), 180 / 100)

    def test_in_radians_per_century_from_year_ra_time(self):
        ra  = RightAscension('12:00:00')
        proper_motion = ProperMotion(ra)
        assert_equal(proper_motion.in_radians_per_century(), pi * 100)

    def test_in_degrees_per_century_from_year_ra_time(self):
        ra  = RightAscension('12:00:00')
        proper_motion = ProperMotion(ra)
        assert_equal(proper_motion.in_degrees_per_century(), 180 * 100)



    # Generate proper motion for right ascension of degrees in radians
    # Return in same time units as entered
    def test_in_radians_per_year_ra_radians(self):
        ra  = RightAscension(radians = pi/2)
        proper_motion = ProperMotion(ra)
        assert_equal(proper_motion.in_radians_per_year(), pi/2)

    def test_in_degrees_per_year_ra_radians(self):
        ra  = RightAscension(radians = pi/2)
        proper_motion = ProperMotion(ra)
        assert_equal(proper_motion.in_degrees_per_year(), 90)

    def test_in_radians_per_century_ra_radians(self):
        ra  = RightAscension(radians = pi/2)
        proper_motion = ProperMotion(ra, time = 'century')
        assert_equal(proper_motion.in_radians_per_century(), pi/2)

    def test_in_degrees_per_century_ra_radians(self):
        ra  = RightAscension(radians = pi/2)
        proper_motion = ProperMotion(ra, time = 'century')
        assert_equal(proper_motion.in_degrees_per_century(), 90)



    # Generate proper motion for right ascension of degrees in radians
    # Return in different time units as entered
    def test_in_radians_per_year_from_century_ra_radians(self):
        ra  = RightAscension(radians = pi/2)
        proper_motion = ProperMotion(ra, time = 'century')
        assert_equal(proper_motion.in_radians_per_year(), pi/2 / 100 )

    def test_in_degrees_per_year_from_century_ra_radians(self):
        ra  = RightAscension(radians = pi/2)
        proper_motion = ProperMotion(ra, time = 'century')
        assert_equal(proper_motion.in_degrees_per_year(), 90 / 100)

    def test_in_radians_per_century_from_year_ra_radians(self):
        ra  = RightAscension(radians = pi/2)
        proper_motion = ProperMotion(ra)
        assert_equal(proper_motion.in_radians_per_century(), pi/2 * 100)

    def test_in_degrees_per_century_from_year_ra_radians(self):
        ra  = RightAscension(radians = pi/2)
        proper_motion = ProperMotion(ra)
        assert_equal(proper_motion.in_degrees_per_century(), 90 * 100)

#_______________________________________________________________________

    # Generate proper motion for declination of degrees in arc
    # Return in same time units as entered
    def test_in_radians_per_year_dec_deg_arc(self):
        dec  = Declination(degrees = 90)
        proper_motion = ProperMotion(dec)
        assert_equal(proper_motion.in_radians_per_year(), pi/2)

    def test_in_degrees_per_year_dec_deg_arc(self):
        dec  = Declination(degrees = 90)
        proper_motion = ProperMotion(dec)
        assert_equal(proper_motion.in_degrees_per_year(), 90)

    def test_in_radians_per_century_dec_deg_arc(self):
        dec  = Declination(degrees = 90)
        proper_motion = ProperMotion(dec, time = 'century')
        assert_equal(proper_motion.in_radians_per_century(), pi/2)

    def test_in_degrees_per_century_dec_deg_arc(self):
        dec  = Declination(degrees = 90)
        proper_motion = ProperMotion(dec, time = 'century')
        assert_equal(proper_motion.in_degrees_per_century(), 90)



    # Generate proper motion for declination of degrees in arc
    # Return in different time units as entered
    def test_in_radians_per_year_from_century_dec_deg_arc(self):
        dec  = Declination(degrees = 90)
        proper_motion = ProperMotion(dec, time = 'century')
        assert_equal(proper_motion.in_radians_per_year(), pi/2 / 100 )

    def test_in_degrees_per_year_from_century_dec_deg_arc(self):
        dec  = Declination(degrees = 90)
        proper_motion = ProperMotion(dec, time = 'century')
        assert_equal(proper_motion.in_degrees_per_year(), 90 / 100)

    def test_in_radians_per_century_from_year_dec_deg_arc(self):
        dec  = Declination(degrees = 90)
        proper_motion = ProperMotion(dec)
        assert_equal(proper_motion.in_radians_per_century(), pi/2 * 100)

    def test_in_degrees_per_century_from_year_dec_deg_arc(self):
        dec  = Declination(degrees = 90)
        proper_motion = ProperMotion(dec)
        assert_equal(proper_motion.in_degrees_per_century(), 90 * 100)



    # Generate proper motion for declination of degrees in time
    # Return in same time units as entered
    def test_in_radians_per_year_dec_time(self):
        dec  = Declination(degrees = 6, units = 'time')
        proper_motion = ProperMotion(dec)
        assert_equal(proper_motion.in_radians_per_year(), pi/2)

    def test_in_degrees_per_year_dec_time(self):
        dec  = Declination(degrees = 6, units = 'time')
        proper_motion = ProperMotion(dec)
        assert_equal(proper_motion.in_degrees_per_year(), 90)

    def test_in_radians_per_century_dec_time(self):
        dec  = Declination(degrees = 6, units = 'time')
        proper_motion = ProperMotion(dec, time = 'century')
        assert_equal(proper_motion.in_radians_per_century(), pi/2)

    def test_in_degrees_per_century_dec_time(self):
        dec  = Declination(degrees = 6, units = 'time')
        proper_motion = ProperMotion(dec, time = 'century')
        assert_equal(proper_motion.in_degrees_per_century(), 90)



    # Generate proper motion for declination of degrees in time
    # Return in different time units as entered
    def test_in_radians_per_year_from_century_dec_time(self):
        dec  = Declination(degrees = 6, units = 'time')
        proper_motion = ProperMotion(dec, time = 'century')
        assert_equal(proper_motion.in_radians_per_year(), pi/2 / 100 )

    def test_in_degrees_per_year_from_century_dec_time(self):
        dec  = Declination(degrees = 6, units = 'time')
        proper_motion = ProperMotion(dec, time = 'century')
        assert_equal(proper_motion.in_degrees_per_year(), 90 / 100)

    def test_in_radians_per_century_from_year_dec_time(self):
        dec  = Declination(degrees = 6, units = 'time')
        proper_motion = ProperMotion(dec)
        assert_equal(proper_motion.in_radians_per_century(), pi/2 * 100)

    def test_in_degrees_per_century_from_year_dec_time(self):
        dec  = Declination(degrees = 6, units = 'time')
        proper_motion = ProperMotion(dec)
        assert_equal(proper_motion.in_degrees_per_century(), 90 * 100)



    # Generate proper motion for declination of degrees in radians
    # Return in same time units as entered
    def test_in_radians_per_year_dec_radians(self):
        dec  = Declination(radians = pi/4)
        proper_motion = ProperMotion(dec)
        assert_equal(proper_motion.in_radians_per_year(), pi/4)

    def test_in_degrees_per_year_dec_radians(self):
        dec  = Declination(radians = pi/4)
        proper_motion = ProperMotion(dec)
        assert_equal(proper_motion.in_degrees_per_year(), 45)

    def test_in_radians_per_century_dec_radians(self):
        dec  = Declination(radians = pi/4)
        proper_motion = ProperMotion(dec, time = 'century')
        assert_equal(proper_motion.in_radians_per_century(), pi/4)

    def test_in_degrees_per_century_dec_radians(self):
        dec  = Declination(radians = pi/4)
        proper_motion = ProperMotion(dec, time = 'century')
        assert_equal(proper_motion.in_degrees_per_century(), 45)



    # Generate proper motion for declination of degrees in radians
    # Return in different time units as entered
    def test_in_radians_per_year_from_century_dec_radians(self):
        dec  = Declination(radians = pi/4)
        proper_motion = ProperMotion(dec, time = 'century')
        assert_equal(proper_motion.in_radians_per_year(), pi/4 / 100 )

    def test_in_degrees_per_year_from_century_dec_radians(self):
        dec  = Declination(radians = pi/4)
        proper_motion = ProperMotion(dec, time = 'century')
        assert_equal(proper_motion.in_degrees_per_year(), 45 / 100)

    def test_in_radians_per_century_from_year_dec_radians(self):
        dec  = Declination(radians = pi/4)
        proper_motion = ProperMotion(dec)
        assert_equal(proper_motion.in_radians_per_century(), pi/4 * 100)

    def test_in_degrees_per_century_from_year_dec_radians(self):
        dec  = Declination(radians = pi/4)
        proper_motion = ProperMotion(dec)
        assert_equal(proper_motion.in_degrees_per_century(), 45 * 100)
