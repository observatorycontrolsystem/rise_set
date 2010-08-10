#!/usr/bin/python

from setuptools import setup, find_packages

setup(
    name = 'rise_set',
    version = '0.2',
    description = 'Routines for accurate rise/set/transit calculations',
    author = 'Eric Saunders',
    author_email = 'esaunders@lcogt.net',
    packages = ['rise_set', 'rise_set.test'],
    dependency_links = [
       "http://pluto/pyslalib-1.0.tar.gz"
    ],
    install_requires = [
      "pyslalib"
    ],
    
    package_data = {
        '': ['*.conf'],    
    },
    
    
    tests_require = [
        "nose"
    ],
    test_suite = 'nose.collector',
)
