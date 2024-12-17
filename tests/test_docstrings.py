import doctest
from gcmotion.entities import (
    tokamak,
    initial_conditions,
    profile,
    particle,
)

V = False  # Verbosity


def test_tokamak_docstring():
    doctest_results = doctest.testmod(m=tokamak, verbose=V)
    assert doctest_results.failed == 0


def test_initial_conditions_docstring():
    doctest_results = doctest.testmod(m=initial_conditions, verbose=V)
    assert doctest_results.failed == 0


def test_profile_docstring():
    doctest_results = doctest.testmod(m=profile, verbose=V)
    assert doctest_results.failed == 0


def test_particle_docstring():
    doctest_results = doctest.testmod(m=particle, verbose=V)
    assert doctest_results.failed == 0
