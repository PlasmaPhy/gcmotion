.. _efields:

Electric field profiles
=================================

Note that the methods :func:`Phi_der()`, :func:`Er` and 
:func:`PhiNU()` as implemented
in :py:class:`~gcmotion.tokamak.efield.ElectricField` are implied.

.. currentmodule:: gcmotion.tokamak.efield

.. Show clickable efield classes

.. rubric:: Currently supported electric fields

.. autosummary::
   :template: _templates/autosummary/class
   :recursive:

   Nofield
   Radial

.. rubric:: Class methods (besides Phi_der(), Er() and PhiNU()):

.. autoclass:: Nofield

.. autoclass:: Radial
    :members: __init__