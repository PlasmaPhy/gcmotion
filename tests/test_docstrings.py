import doctest
from gcmotion.entities import (
    tokamak,
    physical_parameters,
    initial_conditions,
    profile,
)

V = False  # Verbosity


def test_tokamak_docstring():
    doctest_results = doctest.testmod(m=tokamak, verbose=V)
    assert doctest_results.failed == 0


def test_physical_parameters_docstring():
    doctest_results = doctest.testmod(m=physical_parameters, verbose=V)
    assert doctest_results.failed == 0


def test_initial_conditions_docstring():
    doctest_results = doctest.testmod(m=initial_conditions, verbose=V)
    assert doctest_results.failed == 0


def test_profile_docstring():
    doctest_results = doctest.testmod(m=profile, verbose=V)
    assert doctest_results.failed == 0