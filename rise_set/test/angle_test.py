#!/usr/bin/env python

from __future__ import division

from nose.tools import assert_equal, assert_almost_equal, raises

from math import pi

#Import the module to test
from rise_set.angle import Angle, AngleConfigError, InvalidAngleError


class TestAngle(object):
    '''Unit tests for the angle.Angle class.'''

    # Test constructor errors
    @raises(AngleConfigError)
    def test_invalid_angle_none(self):
        self.angle = Angle()

    @raises(AngleConfigError)
    def test_invalid_angle_type(self):
        self.angle = Angle(units = 'time')

    @raises(AngleConfigError)
    def test_invalid_angle_units(self):
        self.angle = Angle(degrees = 45, units = 'kilos')

    @raises(AngleConfigError)
    def test_invalid_multiple_angles(self):
        self.angle = Angle(degrees = 45, radians = pi/4)


    # Test degree functionality
    def test_in_degrees(self):
        self.angle = Angle(degrees=37)
        assert_equal(self.angle.in_degrees(), 37)

    def test_in_degrees_rads_provided(self):
        self.angle = Angle(radians=pi)
        assert_equal(self.angle.in_degrees(), 180)

    def test_in_degrees_negative(self):
        self.angle = Angle(degrees = -37)
        assert_equal(self.angle.in_degrees(), -37)



    # Test radian functionality
    def test_in_radians_rads_provided(self):
        self.angle = Angle(radians=2*pi)
        assert_equal(self.angle.in_radians(), 2*pi)

    def test_in_radians_degrees_provided(self):
        self.angle = Angle(degrees=180)
        assert_equal(self.angle.in_radians(), pi)

    def test_in_radians_degrees_provided_time(self):
        self.angle = Angle(degrees=12, units = 'time')
        assert_equal(self.angle.in_radians(), pi)


    def test_in_radians_negative(self):
        self.angle = Angle(radians = -pi)
        assert_equal(self.angle.in_radians(), -pi)



    # Test valid sexegesimal->degrees conversion given units = time
    def test_from_sexegesimal_hrs_time(self):
        self.angle = Angle(degrees='12:00:00', units= 'time')
        assert_equal(self.angle.in_degrees(), 180)

    def test_from_sexegesimal_more_hrs_time(self):
        self.angle = Angle(degrees='120:00:00', units= 'time')
        assert_equal(self.angle.in_degrees(), 1800)

    def test_from_sexegesimal_hrs_mins_time(self):
        self.angle = Angle(degrees='12:30:00', units = 'time')
        assert_equal(self.angle.in_degrees(), 187.5)

    def test_from_sexegesimal_hrs_secs_time(self):
        self.angle = Angle(degrees='12:30:30', units = 'time')
        assert_equal(self.angle.in_degrees(), 187.625)

    def test_from_sexegesimal_fractional_secs_time(self):
        self.angle = Angle(degrees='12:30:30.1', units = 'time')
        assert_equal(self.angle.in_degrees(), (187.625 + (360/24/36000)))

    def test_from_sexegesimal_negative_time(self):
        self.angle = Angle(degrees='-12:00:00', units = 'time')
        assert_equal(self.angle.in_degrees(), -180)

    def test_from_sexegesimal_positive_time(self):
        self.angle = Angle(degrees='+12:00:00', units = 'time')
        assert_equal(self.angle.in_degrees(), 180)

    def test_from_sexegesimal_zero_time(self):
        self.angle = Angle(degrees='0 0 0', units = 'time')
        assert_equal(self.angle.in_degrees(), 0.0)



    # Test valid sexegesimal->degrees conversion given units = arc
    def test_from_sexegesimal_hrs_arc(self):
        self.angle = Angle(degrees='12:00:00')
        assert_equal(self.angle.in_degrees(), 12)

    def test_from_sexegesimal_more_hrs_arc(self):
        self.angle = Angle(degrees='120:00:00')
        assert_equal(self.angle.in_degrees(), 120)

    def test_from_sexegesimal_hrs_mins_arc(self):
        self.angle = Angle(degrees='12:30:00')
        assert_equal(self.angle.in_degrees(), 12.5)

    def test_from_sexegesimal_hrs_secs_arc(self):
        self.angle = Angle(degrees='12:30:30')
        assert_almost_equal(self.angle.in_degrees(), 12.50833, 5)

    def test_from_sexegesimal_fractional_secs_arc(self):
        self.angle = Angle(degrees='12:30:30.1')
        assert_almost_equal(self.angle.in_degrees(), 12.50836, 5)

    def test_from_sexegesimal_negative_arc(self):
        self.angle = Angle(degrees='-12:00:00')
        assert_equal(self.angle.in_degrees(), -12.0)

    def test_from_sexegesimal_positive_arc(self):
        self.angle = Angle(degrees='+12:00:00')
        assert_equal(self.angle.in_degrees(), 12.0)

    def test_from_sexegesimal_zero_arc(self):
        self.angle = Angle(degrees='0 0 0')
        assert_equal(self.angle.in_degrees(), 0.0)



    # Test various valid input forms for sexegesimal
    def test_from_sexegesimal_valid_format_colons(self):
        self.angle = Angle(degrees='12:30:30', units = 'time')
        assert_equal(self.angle.in_degrees(), 187.625)

    def test_from_sexegesimal_valid_format_one_space(self):
        self.angle = Angle(degrees='12 30 30', units = 'time')
        assert_equal(self.angle.in_degrees(), 187.625)

    def test_from_sexegesimal_valid_format_many_spaces(self):
        self.angle = Angle(degrees='12        30   30', units = 'time')
        assert_equal(self.angle.in_degrees(), 187.625)

    def test_from_sexegesimal_valid_format_weird_delims(self):
        self.angle = Angle(degrees='12$$30$$30', units = 'time')
        assert_equal(self.angle.in_degrees(), 187.625)



    # Test sexegesimal validation on various illegal inputs
    @raises(InvalidAngleError)
    def test_from_sexegesimal_invalid_format_number_mins_too_long(self):
        self.angle = Angle(degrees='12:0120:00')

    @raises(InvalidAngleError)
    def test_from_sexegesimal_invalid_format_secs_too_long(self):
        self.angle = Angle(degrees='12:00:0032')

    @raises(InvalidAngleError)
    def test_from_sexegesimal_invalid_format_min_too_small(self):
        self.angle = Angle(degrees='12:-1:00')

    @raises(InvalidAngleError)
    def test_from_sexegesimal_invalid_format_sec_too_small(self):
        self.angle = Angle(degrees='12:00:-1.0')

    @raises(InvalidAngleError)
    def test_from_sexegesimal_invalid_format_no_delimiters(self):
        self.angle = Angle(degrees = '123030')

    @raises(InvalidAngleError)
    def test_from_sexegesimal_invalid_format_minuses(self):
        self.angle = Angle(degrees = '12-30-30', units = 'time')



    # Test returning degrees in sexegesimal
    def test_in_sexegesimal_degrees_str_arc(self):
        self.angle = Angle(degrees = '12 30 30')
        assert_equal(self.angle.in_sexegesimal(), '12 30 30.0')

    def test_in_sexegesimal_degrees_num_arc(self):
        self.angle = Angle(degrees = 45)
        assert_equal(self.angle.in_sexegesimal(), '45 0 0')

    def test_in_sexegesimal_degrees_negative_num_arc(self):
        self.angle = Angle(degrees = -90)
        assert_equal(self.angle.in_sexegesimal(), '-90 0 0')

    def test_in_sexegesimal_degrees_negative_str_arc(self):
        self.angle = Angle(degrees = '-12 00 00')
        assert_equal(self.angle.in_sexegesimal(), '-12 0 0.0')

    def test_in_sexegesimal_degrees_str_time(self):
        self.angle = Angle(degrees = '12 30 30', units = 'time')
        assert_equal(self.angle.in_sexegesimal(), '12 30 30.0')

    def test_in_sexegesimal_degrees_num_time(self):
        self.angle = Angle(degrees = 12.5, units = 'time')
        assert_equal(self.angle.in_sexegesimal(), '12 30 0.0')

    def test_in_sexegesimal_degrees_negative_num_time(self):
        self.angle = Angle(degrees = -12.5, units = 'time')
        assert_equal(self.angle.in_sexegesimal(), '-12 30 0.0')

    def test_in_sexegesimal_degrees_negative_str_time(self):
        self.angle = Angle(degrees = '-12 00 00', units = 'time')
        assert_equal(self.angle.in_sexegesimal(), '-12 0 0.0')



    # Test returning radians in sexegesimal
    def test_in_sexegesimal_radians_str_arc(self):
        self.angle = Angle(radians = '12 30 30')
        assert_equal(self.angle.in_sexegesimal(radians = True), '12 30 30.0')

    def test_in_sexegesimal_radians_negative_str_arc(self):
        self.angle = Angle(radians = '-2 00 00')
        assert_equal(self.angle.in_sexegesimal(radians = True), '-2 0 0.0')

    def test_in_sexegesimal_radians_str_time(self):
        self.angle = Angle(radians = '12 30 30', units = 'time')
        assert_equal(self.angle.in_sexegesimal(radians = True), '12 30 30.0')

    def test_in_sexegesimal_radians_negative_str_time(self):
        self.angle = Angle(radians = '-3 00 00', units = 'time')
        assert_equal(self.angle.in_sexegesimal(radians = True), '-3 0 0.0')



    # Test converting degrees to radians sexegesimal
    def test_in_sexegesimal_degrees_to_radians_arc(self):
        self.angle = Angle(degrees = 180)
        assert_equal(self.angle.in_sexegesimal(radians = True), '3 8 29.7335529233')

    def test_in_sexegesimal__negative_degrees_to_radians_arc(self):
        self.angle = Angle(degrees = -180)
        assert_equal(self.angle.in_sexegesimal(radians = True), '-3 8 29.7335529233')

    def test_in_sexegesimal_degrees_to_radians_time(self):
        self.angle = Angle(degrees = 12, units = 'time')
        assert_equal(self.angle.in_sexegesimal(radians = True), '3 8 29.7335529233')

    def test_in_sexegesimal_negative_degrees_to_radians_time(self):
        self.angle = Angle(degrees = -12, units = 'time')
        assert_equal(self.angle.in_sexegesimal(radians = True), '-3 8 29.7335529233')

    def test_in_sexegesimal_radians_to_degrees_arc(self):
        self.angle = Angle(radians = pi)
        assert_equal(self.angle.in_sexegesimal(), '180 0 0.0')

    def test_in_sexegesimal_negative_radians_to_degrees_arc(self):
        self.angle = Angle(radians = -2*pi)
        assert_equal(self.angle.in_sexegesimal(), '-360 0 0.0')

    def test_in_sexegesimal_radians_to_degrees_time(self):
        self.angle = Angle(radians = pi, units = 'time')
        assert_equal(self.angle.in_sexegesimal(), '180 0 0.0')

    def test_in_sexegesimal_negative_radians_to_degrees_time(self):
        self.angle = Angle(radians = -2*pi, units = 'time')
        assert_equal(self.angle.in_sexegesimal(), '-360 0 0.0')

