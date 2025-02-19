import pytest
import gcmotion as gcm

initializers = [
    gcm.SmartPositiveInit,
    gcm.SmartNegativeInit,
    gcm.SmartNegative2Init,
    gcm.DTTPositiveInit,
    gcm.DTTPositiveInit,
]

species = "p"

qfactors = [
    gcm.qfactor.SmartPositive,
    gcm.qfactor.SmartNegative,
    gcm.qfactor.SmartNegative2,
    gcm.qfactor.DTTPositive,
    gcm.qfactor.DTTNegative,
]

bfields = [
    gcm.bfield.SmartPositive,
    gcm.bfield.SmartNegative,
    gcm.bfield.SmartNegative2,
    gcm.bfield.DTTPositive,
    gcm.bfield.DTTNegative,
]


@pytest.mark.filterwarnings("ignore:numpy.ndarray size changed")
def test_initializers_initialization():
    r"""Tests that all initializers are initialized."""
    for init in initializers:
        init = init(species)
        init.QuantityConstructor()
        assert hasattr(init, "R")
        assert hasattr(init, "B0")
        assert hasattr(init, "psi_wallNU")
        assert hasattr(init, "a")


@pytest.mark.filterwarnings("ignore:numpy.ndarray size changed")
def test_numerical_qfactor_initialization():
    r"""Tests that all numerical qfactors are initialized."""
    for qfactor in qfactors:
        qfactor = qfactor()
        assert hasattr(qfactor, "is_numerical")


@pytest.mark.filterwarnings("ignore:numpy.ndarray size changed")
def test_numerical_bfield_initialization():
    r"""Tests that all numerical bfields are initialized."""
    for bfield in bfields:
        bfield = bfield()
        assert hasattr(bfield, "is_numerical")
