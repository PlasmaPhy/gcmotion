import doctest

path = "../docs/source/reference/"

V = False  # Verbosity


def test_quantity_reference():
    doctest_results = doctest.testfile(
        filename=f"{path}quantity.rst", verbose=V
    )
    assert doctest_results[0] == 0


def test_qfactor_reference():
    doctest_results = doctest.testfile(
        filename=f"{path}qfactor.rst", verbose=V
    )
    assert doctest_results[0] == 0


def test_bfield_reference():
    doctest_results = doctest.testfile(filename=f"{path}bfield.rst", verbose=V)
    assert doctest_results[0] == 0


def test_efield_reference():
    doctest_results = doctest.testfile(filename=f"{path}efield.rst", verbose=V)
    assert doctest_results[0] == 0
