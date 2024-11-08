.. _qfactors:

q-factor profiles
===========================

Note that the methods :func:`q_of_psi()` and :func:`psipNU()` as implemented
in :py:class:`~gcmotion.tokamak.qfactor.QFactor` are implied.

.. currentmodule:: gcmotion.tokamak.qfactor

.. Show clickable q factor classes

.. rubric:: Currently supported q-factors
   
.. autosummary::
   :template: _templates/autosummary/class
   :recursive:

   Unity
   Parabolic
   Hypergeometric

.. rubric:: Class methods (besides q_of_psi() and psipNU()):

.. autoclass:: Unity

.. autoclass:: Parabolic
   :members: __init__

.. autoclass:: Hypergeometric
   :members: __init__