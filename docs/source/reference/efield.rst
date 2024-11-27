############################
Electric field configuration
############################

.. currentmodule:: gcmotion.efield

Here is a list of the availiable q-factor configurations:

======================     ===================
No Electric field          :py:class:`Nofield`
Raidial Electric field     :py:class:`Radial`
======================     ===================

Their parameters are documented below.

.. rubric:: Example

>>> import gcmotion as gcm
>>>
>>> # Quantity Constructor
>>> Rnum = 1.6
>>> anum = 0.5
>>> B0num = 1
>>> species = "p"
>>> Q = gcm.QuantityConstructor(R=Rnum, a=anum, B0=B0num, species=species)
>>>
>>> # Intermediate values
>>> a = Q(anum, "meters")
>>> B0 = Q(B0num, "Tesla")
>>> Ea = Q(75, "kiloVolts/meter")
>>>
>>> # Some qfactors
>>> efield1 = gcm.efield.Nofield()
>>> efield2 = gcm.efield.Radial(a=a, B0=B0, Ea=Ea, peak=0.98, rw=1/50)

.. note::

    The values of `a` and `B0` are necessary whenever we define the electric
    field with reference to its value at the wall.

*******
Methods
*******

The functions `solverqNU`, `psipNU` and `Er` work identically in every class, 
so I list their methods here as to not repeat myself:

.. autofunction:: gcmotion.efield.ElectricField.solverPhiderNU

.. autofunction:: gcmotion.efield.ElectricField.PhiNU

.. autofunction:: gcmotion.efield.ElectricField.Er

.. _available_efields:

*************************
Available electric fields
*************************

.. rubric:: Nofield

.. autoclass:: Nofield
   :class-doc-from: class

.. rubric:: Radial

.. autoclass:: Radial
   :class-doc-from: class
