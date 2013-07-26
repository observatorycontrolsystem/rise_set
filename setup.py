#!/usr/bin/env python

'''This file controls the packaging and deployment of the rise_set module.

    1) RUN TESTS: Unit tests may be run with

       python setup.py test

    (ignore any logger errors)

    2) CREATE SOURCE: Do
       python setup.py sdist

       to create a source tarball suitable for distributing.

    3) INSTALL:
       python setup.py install

       This will download and compile any needed dependencies, then install to
       your Python package directory.

Eric Saunders
August 2012
'''


from setuptools import setup, find_packages

setup(
    name = 'rise_set',
    version = '0.2.9',
    description = 'Routines for accurate rise/set/transit calculations',
    author = 'Eric Saunders',
    author_email = 'esaunders@lcogt.net',
#    packages = ['rise_set', 'rise_set.test'],
    packages = ['rise_set'],
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
