Electric field profiles
=================================

Note that the methods :func:`Phi_der()`, :func:`Er_of_psi` and 
:func:`Phi_of_psi()` as implemented
in :py:class:`~gcmotion.tokamak.efield.ElectricField` are implied.

.. currentmodule:: gcmotion.tokamak.efield

.. Show clickable q factor classes

.. rubric:: Currently supported electric fields

.. autosummary::
   :template: _templates/autosummary/class
   :recursive:

   Nofield
   Parabolic
   Radial

.. rubric:: Class methods (besides Phi_der(), Er_of_psi() and Phi_of_psi()):

.. autoclass:: Nofield

.. autoclass:: Parabolic
    :members: __init__

.. autoclass:: Radial
    :members: __init__