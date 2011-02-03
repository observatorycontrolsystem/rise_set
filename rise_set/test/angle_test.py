#!/usr/bin/python
from __future__ import division

from nose.tools import eq_, assert_equal, raises

from math import pi

#Import the module to test
from rise_set.angle import Angle, AngleConfigError, InvalidAngleError


class TestAngle(object):
    '''Unit tests for the angle.Angle class.'''

    def setUp(self):
        pass

    def tearDown(self):
        pass

    @raises(AngleConfigError)
    def test_invalid_angle_type(self):
        self.angle = Angle(flibflobs=45)


    def test_in_degrees(self):
        self.angle = Angle(degrees=37)
        assert_equal(self.angle.in_degrees(), 37)

    # Test radian functionality
    def test_in_radians_rads_provided(self):
        self.angle = Angle(radians=2*pi)
        assert_equal(self.angle.in_radians(), 2*pi)

    def test_in_radians_rads_degrees_provided(self):
        self.angle = Angle(degrees=180)
        assert_equal(self.angle.in_radians(), pi)



    # Test valid sexegesimal->degrees conversion
    def test_from_sexegesimal_hrs(self):
        self.angle = Angle(ra='12:00:00')
        assert_equal(self.angle.in_degrees(), 180)

    def test_from_sexegesimal_hrs_mins(self):
        self.angle = Angle(ra='12:30:00')
        assert_equal(self.angle.in_degrees(), 187.5)

    def test_from_sexegesimal_hrs_mins(self):
        self.angle = Angle(ra='12:30:30')
        assert_equal(self.angle.in_degrees(), 187.625)

    def test_from_sexegesimal_fractional_secs(self):
        self.angle = Angle(ra='12:30:30.1')
        assert_equal(self.angle.in_degrees(), (187.625 + (360/24/36000)))



    # Test various valid input forms for sexegesimal
    def test_from_sexegesimal_valid_format_no_delimiters(self):
        self.angle = Angle(ra='123030')
        assert_equal(self.angle.in_degrees(), 187.625)

    def test_from_sexegesimal_valid_format_colons(self):
        self.angle = Angle(ra='12:30:30')
        assert_equal(self.angle.in_degrees(), 187.625)

    def test_from_sexegesimal_valid_format_one_space(self):
        self.angle = Angle(ra='12 30 30')
        assert_equal(self.angle.in_degrees(), 187.625)

    def test_from_sexegesimal_valid_format_many_spaces(self):
        self.angle = Angle(ra='12        30   30')
        assert_equal(self.angle.in_degrees(), 187.625)

    def test_from_sexegesimal_valid_format_weird_delims(self):
        self.angle = Angle(ra='12$$30$$30')
        assert_equal(self.angle.in_degrees(), 187.625)



    # Test sexegesimal validation on various illegal inputs
    def test_from_sexegesimal_invalid_no_ra_or_dec_specified(self):
        try:
            self.angle = Angle()
            self.angle.from_sexegesimal('12:30:30')
            assert False
        except InvalidAngleError as e:
            expected_message = ("Sexegesimal number must be qualified with "
                                "either 'ra' or 'dec'")
            assert_equal(e.__str__(), expected_message)
            assert_equal(e.value, expected_message)

    def test_from_sexegesimal_invalid_format_minuses(self):
        try:
            self.angle = Angle()
            self.angle.from_sexegesimal('12-30-30', ra=True)
            assert False
        except InvalidAngleError as e:
            expected_message = ("Invalid sexagesimal format: Try colon or space"
                                " delimiters instead (e.g. -12:34:56)")
            assert_equal(e.__str__(), expected_message)
            assert_equal(e.value, expected_message)


    # Test sexegesimal RA validation
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
        self.angle = Angle(ra='12:00:-19')


    # Test various valid declination inputs
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
        self.angle = Angle(dec='90:00:01')
