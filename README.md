# Rise-set Library

![Build](https://github.com/observatorycontrolsystem/rise_set/workflows/Build/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/observatorycontrolsystem/rise_set/badge.svg?branch=master)](https://coveralls.io/github/observatorycontrolsystem/rise_set?branch=master)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/8b05ac4107534c7297dfff464360af57)](https://www.codacy.com/gh/observatorycontrolsystem/rise_set/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=observatorycontrolsystem/rise_set&amp;utm_campaign=Badge_Grade)
[![PyPI](https://img.shields.io/pypi/v/ocs-rise-set?style=flat)](https://pypi.org/project/ocs-rise-set/)

***WARNING: The API of this library will be updated in the near future***

This library provides Python routines for finding the positions of astronomical bodies to reasonable precision. It is primarily used to calculate target uptime and sunrise, sunset and twilight times, accurate to 30s-1m. It also supports calculating target airmass over time. The library depends on Fortran SLALIB to perform the computations. The rise/set/transit algorithms are implementations of Astronomical Algorithms, Ch. 14 (Jean Meeus). This library was initially authored by [Eric Saunders](https://github.com/ire-and-curses).

## Prerequisites

-   Python >= 3.7
-   Ability to compile fortran (gfortran installed)

## Installation

It is highly recommended that you install and run your python code inside a dedicated python
[virtual environment](https://docs.python.org/3/tutorial/venv.html).

Add the `ocs_rise_set` package to your python environment (you may need to install numpy first):

```bash
(venv) $ pip install numpy; pip install ocs_rise_set
```

## For Developers

### Running the Tests

The unit tests are currently using nosetests, so there are a few tests dependencies to install. After cloning this project, from the project root and inside your virtual environment (using Poetry):

```bash
$ poetry env use python3.10
$ poetry install
$ poetry run pytest
```
