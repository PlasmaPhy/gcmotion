.. _bfield_configuration:

###############
gcmotion.bfield
###############

.. currentmodule:: gcmotion.bfield

Here is a list of the availiable q-factor configurations:

==================               ===============
Large Aspect Ratio               :py:class:`LAR`
==================               ===============

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
>>> B0 = Q(B0num, "Tesla")
>>> i = Q(0, "NUPlasma_current")
>>> g = Q(1, "NUPlasma_current")
>>>
>>> # A magnetic field
>>> bfield1 = gcm.bfield.LAR(B0=B0, i=i, g=g)

*******
Methods
*******

The functions `bigNU` and `solverbNU` work identically in every class, so I
list their methods here as to not repeat myself:

.. autofunction:: gcmotion.bfield.MagneticField.bigNU
   :no-index:

.. autofunction:: gcmotion.bfield.MagneticField.solverbNU
   :no-index:

.. _available_bfields:

*************************
Available magnetic fields
*************************

.. rubric:: LAR

.. autoclass:: LAR
   :class-doc-from: class

.. rubric:: SmartPositive

.. autoclass:: SmartPositive
   :class-doc-from: class

.. rubric:: SmartNegative

.. autoclass:: SmartNegative
   :class-doc-from: class
