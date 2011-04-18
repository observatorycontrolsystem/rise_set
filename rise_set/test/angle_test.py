#!/usr/bin/python
from __future__ import division

from nose.tools import eq_, assert_equal, assert_almost_equal, raises

from math import pi

#Import the module to test
from rise_set.angle import Angle, AngleConfigError, InvalidAngleError


class TestAngle(object):
    '''Unit tests for the angle.Angle class.'''

    def setUp(self):
        pass

    def tearDown(self):
        pass

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



    # Test degree functionality
    def test_in_degrees(self):
        self.angle = Angle(degrees=37)
        assert_equal(self.angle.in_degrees(), 37)
    
    def test_in_degrees_rads_provided(self):
        self.angle = Angle(radians=pi)
        assert_equal(self.angle.in_degrees(), 180)
    

    # Test radian functionality
    def test_in_radians_rads_provided(self):
        self.angle = Angle(radians=2*pi)
        assert_equal(self.angle.in_radians(), 2*pi)

    def test_in_radians_degrees_provided(self):
        self.angle = Angle(degrees=180)
        assert_equal(self.angle.in_radians(), pi)



    # Test valid sexegesimal->degrees conversion given units = time
    def test_from_sexegesimal_hrs_time(self):
        self.angle = Angle(degrees='12:00:00', units= 'time')
        assert_equal(self.angle.in_degrees(), 180)

    def test_from_sexegesimal_hrs_mins_time(self):
        self.angle = Angle(degrees='12:30:00', units = 'time')
        assert_equal(self.angle.in_degrees(), 187.5)

    def test_from_sexegesimal_hrs_secs_time(self):
        self.angle = Angle(degrees='12:30:30', units = 'time')
        assert_equal(self.angle.in_degrees(), 187.625)

    def test_from_sexegesimal_fractional_secs_time(self):
        self.angle = Angle(degrees='12:30:30.1', units = 'time')
        assert_equal(self.angle.in_degrees(), (187.625 + (360/24/36000)))



    # Test valid sexegesimal->degrees conversion given units = arc
    def test_from_sexegesimal_hrs_arc(self):
        self.angle = Angle(degrees='12:00:00')
        assert_equal(self.angle.in_degrees(), 12)

    def test_from_sexegesimal_hrs_mins_arc(self):
        self.angle = Angle(degrees='12:30:00')
        assert_equal(self.angle.in_degrees(), 12.5)

    def test_from_sexegesimal_hrs_secs_arc(self):
        self.angle = Angle(degrees='12:30:30')
        assert_almost_equal(self.angle.in_degrees(), 12.50833, 5)

    def test_from_sexegesimal_fractional_secs_arc(self):
        self.angle = Angle(degrees='12:30:30.1')
        assert_almost_equal(self.angle.in_degrees(), (12.50836), 5)
        


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
    def test_from_sexegesimal_invalid_no_delimiters(self):###
        try:
            self.angle = Angle(degrees = '123030')
            assert False
        except InvalidAngleError as e:
            expected_message = ("Invalid sexegesimal format: Try colon or space"
                                " delimiters instead (e.g. -12:34:56)")
            assert_equal(e.value, expected_message)

    def test_from_sexegesimal_invalid_format_minuses(self):
        try:
            self.angle = Angle(degrees = '12-30-30', units = 'time')
            assert False
        except InvalidAngleError as e:
            expected_message = ("Invalid sexegesimal format: Try colon or space"
                                " delimiters instead (e.g. -12:34:56)")
            assert_equal(e.value, expected_message)


    # Test returning degrees in sexegesimal
    def test_in_sexegesimal_degree_string_provided_arc(self):
        self.angle = Angle(degrees = '12 30 30')
        assert_equal(self.angle.in_sexegesimal(), '12 30 30.0')
    
    def test_in_sexegesimal_degree_number_provided_arc(self):
        self.angle = Angle(degrees = 45)
        assert_equal(self.angle.in_sexegesimal(), '45 00 00')
        
    def test_in_sexegesimal_degree_string_provided_time(self):
        self.angle = Angle(degrees = '12 30 30', units = 'time')
        assert_equal(self.angle.in_sexegesimal(), '12 30 30.0')
        
    def test_in_sexegesimal_degree_number_provided_time(self):
        self.angle = Angle(degrees = 45, units = 'time')
        assert_equal(self.angle.in_sexegesimal(), '03 00 00.0')
        
    
    # Test returning radians in sexegesimal
    def test_in_sexegesimal_radians_string_provided_arc(self):
        self.angle = Angle(radians = '12 30 30')
        assert_equal(self.angle.in_sexegesimal(), '12 30 30.0')
        
    def test_in_sexegesimal_radians_number_provided_arc(self):
        self.angle = Angle(radians = pi)
        assert_equal(self.angle.in_sexegesimal(), '180 00 00.0')            
        
    def test_in_sexegesimal_radians_string_provided_time(self):
        self.angle = Angle(radians = '12 30 30', units = 'time')
        assert_equal(self.angle.in_sexegesimal(), '12 30 30.0')
        
    def test_in_sexegesimal_radians_number_provided_time(self):
        self.angle = Angle(radians = pi, units = 'time')
        assert_equal(self.angle.in_sexegesimal(), '12 00 00.0') 
        

    '''# Test sexegesimal RA validation
    @raises(InvalidAngleError)
    def test_from_sexegesimal_invalid_format_number_hrs_too_long(self):
        self.angle = Angle(ra='1222:12:00')

    @raises(InvalidAngleError)
    def test_from_sexegesimal_invalid_format_number_mins_too_long(self):
        self.angle = Angle(ra='12:0120:00')

    @raises(InvalidAngleError)
    def test_from_sexegesimal_invalid_format_secs_too_long(self):
        self.angle = Angle(ra='12:00:0032')

    @raises(InvalidAngleError)
    def test_from_sexegesimal_invalid_hr_too_small(self):
        self.angle = Angle(ra='-12:00:00')

    @raises(InvalidAngleError)
    def test_from_sexegesimal_invalid_hr_too_big(self):
        self.angle = Angle(ra='24:00:00')

    @raises(InvalidAngleError)
    def test_from_sexegesimal_invalid_min_too_big(self):
        self.angle = Angle(ra='12:75:00')

    @raises(InvalidAngleError)
    def test_from_sexegesimal_invalid_min_too_small(self):
        self.angle = Angle(ra='12:-10:00')

    @raises(InvalidAngleError)
    def test_from_sexegesimal_invalid_sec_too_big(self):
        self.angle = Angle(ra='12:00:99')

    @raises(InvalidAngleError)
    def test_from_sexegesimal_invalid_sec_too_small(self):
        self.angle = Angle(ra='12:00:-19')'''


    '''# Test various valid declination inputs
    def test_from_sexegesimal_dec_north_pole(self):
        self.angle = Angle(dec='+90:00:00')
        assert_equal(self.angle.in_degrees(), 90)

    def test_from_sexegesimal_dec_south_pole(self):
        self.angle = Angle(dec='-90:00:00')
        assert_equal(self.angle.in_degrees(), -90)

    def test_from_sexegesimal_dec_positive_number(self):
        self.angle = Angle(dec='30:30:01')
        assert_equal(self.angle.in_degrees(), (30.5 + 1/3600))

    def test_from_sexegesimal_dec_negative_number(self):
        self.angle = Angle(dec='-30:30:01')
        assert_equal(self.angle.in_degrees(), (-30.5 - 1/3600))

    def test_from_sexegesimal_dec_small_positive_number(self):
        self.angle = Angle(dec='00:06:01')
        assert_equal(self.angle.in_degrees(), (0.1 + 1/3600))

    def test_from_sexegesimal_dec_small_negative_number(self):
        self.angle = Angle(dec='-00:06:01')
        assert_equal(self.angle.in_degrees(), (-0.1 - 1/3600))


    # Test sexegesimal dec validation
    def test_from_sexegesimal_invalid_format_number_too_big_mins(self):
        try:
            self.angle = Angle(dec='90:10:00')
            assert False
        except InvalidAngleError as e:
            expected_message = ("'90 10 0' exceeds maximum allowable"
                                " declination")
            assert_equal(e.value, expected_message)


    @raises(InvalidAngleError)
    def test_from_sexegesimal_invalid_format_number_too_big_secs(self):
        self.angle = Angle(dec='90:00:01')'''
