#!/usr/bin/python

import rise_set.angle as angle
from math import pi

from nose.tools import eq_

def test_from_degrees():
    a = angle.Angle()
    a.from_degrees(37)
    eq_(a.in_degrees(), 37)



class some_test_class_i_just_made_up(object):
    def setup(self):
        self.angle = angle.Angle()
        print "Hello"

    def test_from_degrees(self):
        self.angle.from_degrees(37)
        eq_(self.angle.in_degrees(), 37)


    def test_from_radians(self):
        self.angle.from_radians(pi)
        eq_(self.angle.in_degrees(), 180)
