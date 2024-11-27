#############################
Tokamak Configuration Objects
#############################

The qfactor, bfield and efield modules contain Classes for creating instances of q-factors, magnetic fields and electric fields respectively. Those instaces contain all the information about the object, like their parameters, and implement a method for every needed calculation.

About Tokamak Configuration Objects
===================================

In general, a configuration object looks like this:

.. code-block:: python

    class MyElectricField(ElectricField):

        def __init__(self, *parameters):
            """Parameter setup. Quantites and intermediate
            values are setup here"""

        def solverPhiderNU(self, psi, theta):
            """Return the two Φ derivatives in NU as needed 
            by the solver"""
            return [Phi_der_psip, Phi_der_theta]

        def PhiNU(self, psi, theta):
            """Return the Φ value in NU"""
            return Phi

        def Er(self, psi, theta):
            """Return the radial field's magnitude in NU"""
            return Er

        def __repr__():
            """Object representation"""
            return str

The instantiation parameters must be Quantities in either SI or NU. The `__init__()` method makes sure they have the correct dimensionality and units, calculates intermediate values and stores their magnitudes in different fields. The rest of the methods deal only with purely numerical values and *not* with Quantities, since Quantities greatly impact performance..

All electric field classes inherit the abstract class `ElectricField` to ensure functionality.

Contents
========

.. toctree::
   :maxdepth: 2

   qfactor
   bfield
   efield
