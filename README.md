# Rise-set Library

***WARNING: The API of this library will be updated in the near future***

This library provides Python routines for finding the positions of astronomical bodies to reasonable precision. It is primarily used to calculate target uptime and sunrise, sunset and twilight times, accurate to 30s-1m. It also supports calculating target airmass over time. The library depends on Fortran SLALIB to perform the computations. The rise/set/transit algorithms are implementations of Astronomical Algorithms, Ch. 14 (Jean Meeus). This library was initially authored by Eric Saunders (eric.saunders@gmail.com).

## Prerequisites

-   Python >= 2.7
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

The unit tests are currently using nosetests, so there are a few tests dependencies to install. After cloning this project, from the project root and inside your virtual environment:

```bash
(venv) $ pip install -r test_requirements.pip
(venv) $ nosetests
```
