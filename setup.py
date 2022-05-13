#!/usr/bin/env python

'''This file controls the packaging and deployment of the rise_set module.

    Preferred Installation methods are found in the INSTALLATION file in this directory.

Eric Saunders
August 2012
'''

from setuptools import setup

setup(
    name = 'ocs-rise-set',
    version = '0.5.2',
    description = 'Routines for accurate rise/set/transit calculations',
    author = 'Eric Saunders',
    author_email = 'esaunders@lcogt.net',
    packages = ['rise_set'],
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
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Topic :: Scientific/Engineering :: Astronomy',
    ]

)
