#!/usr/bin/env python

'''This file controls the packaging and deployment of the rise_set module.

    Preferred Installation methods are found in the INSTALLATION file in this directory.

Eric Saunders
August 2012
'''

from setuptools import setup, find_packages

setup(
    name = 'rise_set',
    version = '0.4.1',
    description = 'Routines for accurate rise/set/transit calculations',
    author = 'Eric Saunders',
    author_email = 'esaunders@lcogt.net',
    packages = ['rise_set'],
    dependency_links = [
        "http://github.com/LCOGT/pyslalib/tarball/master#egg=pyslalib",
    ],
    setup_requires = [
        "numpy"
    ],
    install_requires = [
        "pySLALIB",
        "numpy",
        "future",
    ],
    package_data = {
        '': ['*.conf'],
    },

    tests_require = [
        "mock",
        "pylint",
        "nose",
        "nosexcover"
    ],
    test_suite = 'nose.collector',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.3',
    ]

)
