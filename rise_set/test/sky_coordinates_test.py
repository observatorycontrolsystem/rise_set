#!/usr/bin/python
from __future__ import division

from nose.tools import raises

from math import pi

#Import the module to test
from rise_set.sky_coordinates import RightAscension, Declination
from rise_set.angle import InvalidAngleError


class TestCoordinateSystem(object):
    '''Unit tests for the sky_coordinates.RightAscension and 
    sky_coordinates.Declination classes. Note: These classes inherit
    from the angle.Angle class'''
    
    # Test invalid inputs for RA
    @raises(InvalidAngleError)
    def test_validate_ra_invalid_hr_too_small(self):
        self.angle = RightAscension('-12:00:00')

    @raises(InvalidAngleError)
    def test_validate_ra_invalid_hr_too_big(self):
        self.angle = RightAscension('24:00:00')
    
    @raises(InvalidAngleError)
    def test_validate_ra_invalid_hr_too_big_rad(self):
        self.angle = RightAscension(radians = 2*pi)

    @raises(InvalidAngleError)
    def test_validate_ra_invalid_min_too_big(self):
        self.angle = RightAscension('23:61:00')



    # Test invalid inputs for Dec
    @raises(InvalidAngleError)
    def test_validate_dec_invalid_deg_too_small(self):
        self.angle = Declination(-91)

    @raises(InvalidAngleError)
    def test_validate_dec_invalid_deg_too_big(self):
        self.angle = Declination(degrees = 91)
    
    @raises(InvalidAngleError)
    def test_validate_dec_invalid_deg_too_big_rad(self):
        self.angle = Declination(radians = pi)

    @raises(InvalidAngleError)
    def test_validate_dec_invalid_min_too_big(self):
        self.angle = Declination(degrees = '89:61:00')



    # Test valid inputs for RA
    def test_validate_ra_valid_at_zero_degrees(self):
        self.angle = RightAscension(degrees = 0)
        assert(self.angle.validate_ra())
        
    def test_validate_ra_valid_at_zero_radians(self):
        self.angle = RightAscension(radians = 0)
        assert(self.angle.validate_ra())   
        
    def test_validate_ra_valid_at_middle2(self):
        self.angle = RightAscension(degrees = 24)
        assert(self.angle.validate_ra())   
        
    def test_validate_ra_valid_at_middle(self):
        self.angle = RightAscension(degrees = 180)
        assert(self.angle.validate_ra())

    def test_validate_ra_valid_at_middle_rad(self):
        self.angle = RightAscension(radians = pi)
        assert(self.angle.validate_ra())



    # Test valid inputs for Dec
    def test_validate_dec_valid_at_zero(self):
        self.angle = Declination(degrees = 0)
        assert(self.angle.validate_dec())
    
    def test_validate_dec_valid_at_middle(self):
        self.angle = Declination(degrees = 45)
        assert(self.angle.validate_dec())

    def test_validate_dec_valid_at_middle_rad(self):
        self.angle = Declination(radians = pi/4)
        assert(self.angle.validate_dec())


