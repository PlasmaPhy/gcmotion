.. _qfactor_configuration:

################
gcmotion.qfactor
################

.. currentmodule:: gcmotion.qfactor

Here is a list of the availiable q-factor configurations:

=========================      ==========================
Unity q-factor                 :py:class:`Unity`
Parabolic q-factor             :py:class:`Parabolic`
Hypergeometric q-factor        :py:class:`Hypergeometric`
(Numerical) SmartPositive      :py:class:`SmartPositive`
(Numerical) SmartNegative      :py:class:`SmartNegative`
=========================      ==========================

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
>>>
>>> # Some qfactors
>>> qfactor1 = gcm.qfactor.Unity()
>>> qfactor3 = gcm.qfactor.Parabolic(a, B0, q0=1.1, q_wall=3.8)
>>> qfactor3 = gcm.qfactor.Hypergeometric(a, B0, q0=1.1, q_wall=3.8, n=2)

.. note::

    The values of `a` and `B0` are necessary whenever we define the qfactor
    with reference to its value at the wall.

*******
Methods
*******

The functions `solverqNU` and `psipNU` work identically in every class, so I
list their methods here as to not repeat myself:

.. autofunction:: gcmotion.qfactor.QFactor.solverqNU

.. autofunction:: gcmotion.qfactor.QFactor.psipNU

.. _available_qfactors:

*******************
Available q-factors
*******************

.. rubric:: Unity

.. autoclass:: Unity
   :class-doc-from: class

.. rubric:: Parabolic

.. autoclass:: Parabolic
   :class-doc-from: class

.. rubric:: Hypergeometric

.. autoclass:: Hypergeometric
   :class-doc-from: class

.. rubric:: SmartPositive

.. autoclass:: SmartPositive
   :class-doc-from: class

.. rubric:: SmartNegative

.. autoclass:: SmartNegative
   :class-doc-from: class
